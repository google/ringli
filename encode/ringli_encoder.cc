// Copyright 2024 The Ringli Authors. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "encode/ringli_encoder.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include "absl/log/check.h"
#include "common/block_predictor.h"
#include "common/convolve.h"
#include "common/covariance_lattice.h"
#include "common/data_defs/constants.h"
#include "common/data_defs/data_vector.h"
#include "common/dct.h"
#include "common/entropy_coding.h"
#include "common/fast_online_predictor.h"
#include "common/log2floor.h"
#include "common/online_predictor.h"
#include "common/predictor.h"
#include "common/ringli_header.h"
#include "common/wav_header.h"
#include "common/wav_reader.h"
#include "encode/entropy_encode.h"
#include "encode/noise_shaping.h"

namespace ringli {
namespace {

int QuantizationBySignal(double mean_signal,
                         const RingliEncoderConfig& config) {
  return round(config.quantization_curve.GetValue(mean_signal));
}

void CalculateQuantizationMean(const AudioBlock& channels,
                               const RingliEncoderConfig& config,
                               RingliDCTHeader& header) {
  double mean_abs_value = 0;
  for (const auto& input_channel : channels.GetChannels()) {
    for (auto value : input_channel) {
      mean_abs_value += std::abs(value);
    }
  }
  mean_abs_value /= (channels.GetChannels().size() * kRingliBlockSize);

  const int quantization_coef = QuantizationBySignal(mean_abs_value, config);
  for (int i = 0; i < kNumDctBands; ++i) {
    header.quant[i] = quantization_coef;
  }
}

void CalculateQuantizationMeanBand(const AudioBlock& channels,
                                   const RingliEncoderConfig& config,
                                   RingliDCTHeader& header) {
  double mean_abs_value[kNumDctBands] = {0};
  for (const auto& input_channel : channels.GetChannels()) {
    // Calculate mean values for dct bands among all DCT blocks.
    for (int i = 0; i < kRingliBlockSize; i++) {
      const int k = i % kDctLength;
      for (int b = 0; b < kNumDctBands; ++b) {
        if (k < kDctBandEnd[b]) {
          mean_abs_value[b] += std::abs(input_channel[i]);
          break;
        }
      }
    }
  }
  for (int b = 0; b < kNumDctBands; ++b) {
    int bw = kDctBandEnd[b] - (b > 0 ? kDctBandEnd[b - 1] : 0);
    mean_abs_value[b] /= bw * kDctNumber * channels.GetChannels().size();
    header.quant[b] = QuantizationBySignal(mean_abs_value[b], config);
  }
}

void CalculateQuantization(const AudioBlock& channels,
                           const RingliEncoderConfig& config,
                           RingliDCTHeader& header) {
  switch (config.quantization_type) {
    case QuantizationType::CONST: {
      for (int i = 0; i < kNumDctBands; ++i) {
        header.quant[i] = config.quantization_curve.GetValue(0);
      }
      break;
    }
    case QuantizationType::ADAPTIVE_PER_PACKET:
      CalculateQuantizationMean(channels, config, header);
      break;
    case QuantizationType::ADAPTIVE_PER_BAND:
      CalculateQuantizationMeanBand(channels, config, header);
      break;
  }
}

RingliBlock EncodePredictive(const RingliEncoderConfig& config,
                             const AudioBlock& block) {
  const size_t num_channels = block.GetChannels().size();
  const int quant = config.dconfig.pred_quant;
  const int order_min = config.pred_order_min;
  const int order_max = config.pred_order_max;
  CHECK_GE(order_min, 2);
  CHECK_LE(order_max, kMaxPredictorOrder);
  RingliBlock encoded_block(num_channels);
  double sample;
  const auto& num_bits = [&](int residual) {
    return residual == 0 ? 0 : Log2FloorNonZero(std::abs(residual) + 1);
  };
  for (size_t c = 0; c < num_channels; ++c) {
    if (config.dconfig.use_online_predictive_coding) {
      const bool fast_mode = config.dconfig.predictor_fast_mode();
      std::unique_ptr<Predictor> predictor;
      if (fast_mode) {
        predictor = std::make_unique<FastOnlinePredictor>();
      } else {
        predictor =
            std::make_unique<OnlinePredictor>(kOnlinePredictorRegulariser);
      }
      NoiseShaper noise_shaper;
      const float iquant = 1.0 / quant;
      for (int i = 0; i < kRingliBlockSize; i++) {
        float prediction = predictor->Predict();
        if (quant == 1) prediction = std::round(prediction);
        sample = block[c][i];
        if (config.use_noise_shaping) {
          float filtered_noise = noise_shaper.GetFilteredNoise();
          sample += filtered_noise;
        }
        const float error = sample - prediction;
        encoded_block.channels[c][i] = std::round(error * iquant);
        const float sample_deq =
            prediction + quant * encoded_block.channels[c][i];
        predictor->AddNewSample(sample_deq);
        if (config.use_noise_shaping) {
          noise_shaper.AddNewSample(sample_deq - sample);
        }
      }
    } else {
      RingliVector best_residuals;
      RingliPredictiveHeader best_header;
      bool last_is_best = true;
      double best_score = 0;
      double regulariser = (quant * quant) * (1.0 / 12.0);
      CovarianceLattice<int32_t> covlattice_orig(
          block[c].Data(), kRingliBlockSize, kMaxPredictorOrder, regulariser);
      for (int order = order_min; order <= order_max; order += 2) {
        RingliPredictiveHeader& header = encoded_block.header.pred[c];
        int total_num_bits = 0;
        float iquant = 1.0 / quant;
        BlockPredictor<kRingliBlockSize> block_predictor =
            BlockPredictor<kRingliBlockSize>::CreateForEncoder(order, &header,
                                                               covlattice_orig);
        if (quant == 1) {
          for (int i = 0; i < kRingliBlockSize; i++) {
            const float prediction = std::round(block_predictor.Predict());
            block_predictor.AddNewSample(block[c][i]);
            encoded_block.channels[c][i] = block[c][i] - prediction;
            total_num_bits += num_bits(encoded_block.channels[c][i]);
          }
        } else {
          for (int i = 0; i < kRingliBlockSize; i++) {
            const float prediction = block_predictor.Predict();
            const float error = block[c][i] - prediction;
            encoded_block.channels[c][i] = std::round(error * iquant);
            const float sample_deq =
                prediction + quant * encoded_block.channels[c][i];
            block_predictor.AddNewSample(sample_deq);
            total_num_bits += num_bits(encoded_block.channels[c][i]);
          }
        }
        double score = total_num_bits + order * 1.5;
        if (best_score == 0 || score < best_score) {
          best_score = score;
          if (order < order_max) {
            best_header = header;
            best_residuals = encoded_block.channels[c];
          }
          last_is_best = true;
        } else {
          last_is_best = false;
        }
      }
      if (!last_is_best) {
        encoded_block.channels[c] = best_residuals;
        encoded_block.header.pred[c] = best_header;
      }
    }
  }
  return encoded_block;
}

RingliBlock EncodeWithDCT(const RingliEncoderConfig& config,
                          const AudioBlock& prev, const AudioBlock& current,
                          const AudioBlock& next) {
  const size_t num_channels = current.GetChannels().size();
  RingliBlock encoded_block(num_channels);

  RingliDCTHeader prev_header;
  RingliDCTHeader curr_header;
  RingliDCTHeader next_header;
  CalculateQuantization(prev, config, prev_header);
  CalculateQuantization(current, config, curr_header);
  CalculateQuantization(next, config, next_header);

  for (size_t c = 0; c < num_channels; ++c) {
    DataVector<float, kACPredictionWindowSize> coeff_window;
    for (int i = 0; i < kACPredictionWindowSize; ++i) {
      if (i < kACPredictionBorder) {
        coeff_window[i] = prev[c][kRingliBlockSize - kACPredictionBorder + i];
      } else if (i < kRingliBlockSize + kACPredictionBorder) {
        coeff_window[i] = current[c][i - kACPredictionBorder];
      } else {
        coeff_window[i] = next[c][i - kRingliBlockSize - kACPredictionBorder];
      }
    }
    for (int i = 0; i < kACPredictionWindowSize; i += kDctLength) {
      DataVector<float, kDctLength> signal_data;
      for (int k = 0; k < kDctLength; ++k) {
        signal_data[k] = coeff_window[i + k];
      }
      const DataVector<float, kDctLength> dct_data = ForwardDCT(signal_data);
      for (int k = 0; k < kDctLength; ++k) {
        coeff_window[i + k] = dct_data[k];
      }
    }
    DataVector<float, kACPredictionWindowSize> predictor_window;
    DataVector<float, kACPredictionWindowSize> output_window;
    int prev_k_limit = 0;
    for (int step = 0; step <= kNumACPredictionSteps; ++step) {
      int k_limit = step == kNumACPredictionSteps ? kDctLength
                                                  : kACPredictionStart << step;
      for (int i = 0; i < kACPredictionWindowSize; i += kDctLength) {
        for (int k = prev_k_limit; k < k_limit; ++k) {
          float quant = i < kACPredictionBorder
                            ? prev_header.GetQuantizationCoef(k)
                        : i < kRingliBlockSize + kACPredictionBorder
                            ? curr_header.GetQuantizationCoef(k)
                            : next_header.GetQuantizationCoef(k);
          float qcoef = coeff_window[i + k] / quant;
          int32_t icoef = std::round(qcoef);
          if (i >= kACPredictionBorder &&
              i < kRingliBlockSize + kACPredictionBorder) {
            encoded_block.channels[c][i - kACPredictionBorder + k] = icoef;
          }
          if (step < kNumACPredictionSteps) {
            coeff_window[i + k] = icoef * quant + predictor_window[i + k];
          }
        }
        if (step < kNumACPredictionSteps) {
          DataVector<float, kDctLength> dct_data;
          for (int k = 0; k < k_limit; ++k) {
            dct_data[k] = coeff_window[i + k];
          }
          const DataVector<float, kDctLength> signal_data =
              InverseDCT(dct_data);
          for (int k = 0; k < kDctLength; ++k) {
            output_window[i + k] = signal_data[k];
          }
        }
      }
      if (step == kNumACPredictionSteps) {
        break;
      }
      output_window =
          Convolve(GaussianKernel<kDctLength>(kACPredictionSigma / k_limit),
                   output_window);
      for (int i = 0; i < kACPredictionWindowSize; i += kDctLength) {
        DataVector<float, kDctLength> signal_data;
        for (int k = 0; k < kDctLength; ++k) {
          signal_data[k] = output_window[i + k];
        }
        const DataVector<float, kDctLength> dct_data = ForwardDCT(signal_data);
        for (int k = k_limit; k < 2 * k_limit; ++k) {
          predictor_window[i + k] += dct_data[k];
          coeff_window[i + k] -= dct_data[k];
        }
      }
      prev_k_limit = k_limit;
    }
  }
  encoded_block.header.dct = curr_header;
  return encoded_block;
}

}  // namespace

StreamingRingliEncoder::StreamingRingliEncoder(
    const RingliEncoderConfig& config)
    : config_(config) {
  format_.format_chunk_size = 0;
  wav_reader_.RegisterCallback("fmt ", this, ParseFormatCb,
                               sizeof(FormatChunk));
}

void StreamingRingliEncoder::Reset() {
  ringli_data_.clear();
  wav_reader_.Reset();
  format_.format_chunk_size = 0;
  output_pos_ = 0;
  ringli_blocks_.clear();
  ringli_headers_.clear();
}

bool StreamingRingliEncoder::ParseFormatCb(void* opaque, const uint8_t* data,
                                           size_t len, size_t chunk_pos,
                                           size_t chunk_size) {
  auto enc = reinterpret_cast<StreamingRingliEncoder*>(opaque);
  if (chunk_size != len) {
    fprintf(stderr, "Invalid format chunk size %zu\n", chunk_size);
    return false;
  }
  enc->format_.format_chunk_size = chunk_size;
  size_t pos = 0;
  if (!ParseFormatChunk(data, len, &pos, &enc->format_)) {
    return false;
  }
  enc->InitForFormat();
  return true;
}

void StreamingRingliEncoder::InitForFormat() {
  const bool fully_streaming = config_.dconfig.use_online_predictive_coding &&
                               config_.dconfig.ecparams.arithmetic_only;
  const size_t num_channels = format_.number_of_channels;
  const size_t bytes_per_sample = format_.bits_per_sample / 8;
  const size_t block_size = (fully_streaming ? 1 : kRingliBlockSize) *
                            num_channels * bytes_per_sample;
  wav_reader_.RegisterCallback("data", this, ProcessDataCb, block_size);
  if (!config_.dconfig.use_predictive_coding) {
    prev_ = std::make_unique<AudioBlock>(num_channels);
    current_ = std::make_unique<AudioBlock>(num_channels);
    next_ = std::make_unique<AudioBlock>(num_channels);
  } else {
    entropy_coder_ = std::make_unique<EntropyCoder>(
        config_.dconfig.ecparams, format_.sampling_frequency, num_channels,
        config_.dconfig.use_predictive_coding,
        config_.dconfig.use_online_predictive_coding);
    if (fully_streaming) {
      const bool fast_mode = config_.dconfig.predictor_fast_mode();
      idx_ = 0;
      for (int i = 0; i < num_channels; ++i) {
        if (fast_mode) {
          predictors_.emplace_back(std::make_unique<FastOnlinePredictor>());
        } else {
          predictors_.emplace_back(
              std::make_unique<OnlinePredictor>(kOnlinePredictorRegulariser));
        }
      }

      noise_shapers_.resize(num_channels);
      adaptive_quantizers_.resize(num_channels);

      for (size_t c = 0; c < num_channels; ++c) {
        predictors_[c]->Reset();
        noise_shapers_[c].Reset();
        adaptive_quantizers_[c].Reset();
      }
    }
  }
}

bool StreamingRingliEncoder::ProcessDataCb(void* opaque, const uint8_t* data,
                                           size_t len, size_t chunk_pos,
                                           size_t chunk_size) {
  auto enc = reinterpret_cast<StreamingRingliEncoder*>(opaque);
  if (enc->format_.format_chunk_size == 0) {
    fprintf(stderr, "Missing fmt chunk.\n");
    return false;
  }
  return enc->ProcessData(data, len, chunk_pos, chunk_size);
}

bool StreamingRingliEncoder::ProcessData(const uint8_t* data, size_t len,
                                         size_t chunk_pos, size_t chunk_size) {
  if (chunk_pos == 0) {
    WriteHeader(chunk_size);
  }
  if (config_.dconfig.use_predictive_coding) {
    if (config_.dconfig.use_online_predictive_coding &&
        config_.dconfig.ecparams.arithmetic_only) {
      const size_t num_channels = format_.number_of_channels;
      const size_t bytes_per_sample = format_.bits_per_sample / 8;
      std::vector<int> encoded(num_channels);
      for (int ci = 0; ci < num_channels; ++ci) {
        int16_t value;
        memcpy(&value, &data[ci * bytes_per_sample], bytes_per_sample);
        float sample = static_cast<float>(value);
        if (config_.use_noise_shaping) {
          const float filtered_noise = noise_shapers_[ci].GetFilteredNoise();
          sample += filtered_noise;
        }

        float quant;
        if (config_.dconfig.use_adaptive_quantization) {
          quant = adaptive_quantizers_[ci].QuantStep();
        } else {
          quant = config_.dconfig.pred_quant;
        }
        const float iquant = 1.0f / quant;
        // printf("idx %d, c %d, quant: %f\n", int(idx_), ci, quant);

        const float prediction = predictors_[ci]->Predict();
        const float error = sample - prediction;
        encoded[ci] = std::round(error * iquant);
        const float sample_deq = prediction + quant * encoded[ci];
        predictors_[ci]->AddNewSample(sample_deq);
        if (config_.dconfig.use_adaptive_quantization) {
          adaptive_quantizers_[ci].ProcessSample(sample_deq);
        }
        if (config_.use_noise_shaping) {
          noise_shapers_[ci].AddNewSample(sample_deq - sample);
        }
      }
      entropy_coder_->ProcessSamples(encoded.data(), &ringli_data_);
      ++idx_;
      if (idx_ == kRingliBlockSize) {
        idx_ = 0;
        for (int ci = 0; ci < num_channels; ++ci) {
          predictors_[ci]->Reset();
          adaptive_quantizers_[ci].Reset();
        }
      }
    } else {
      AudioBlock block(format_.number_of_channels);
      CopyBlock(data, len, &block);
      return entropy_coder_->ProcessBlock(EncodePredictive(config_, block),
                                          &ringli_data_);
    }
  } else if (chunk_pos == 0) {
    CopyBlock(data, len, current_.get());
  } else {
    CopyBlock(data, len, next_.get());
    ProcessBlock(EncodeWithDCT(config_, *prev_, *current_, *next_));
    *prev_ = *current_;
    *current_ = *next_;
  }
  return true;
}

void StreamingRingliEncoder::WriteHeader(size_t chunk_size) {
  data_len_ = chunk_size;
  RingliHeader ringli_header;
  memcpy(ringli_header.ringli_id, kRingliId, sizeof(kRingliId));
  ringli_header.header_length =
      sizeof(RingliHeader) -
      (sizeof(ringli_header.ringli_id) + sizeof(ringli_header.header_length));
  ringli_header.number_of_channels = format_.number_of_channels;
  ringli_header.sampling_frequency = format_.sampling_frequency;
  ringli_header.bits_per_sample = format_.bits_per_sample;
  ringli_header.data_length = chunk_size;
  // Save encoder config fields that are needed for decoding in ringli header.
  ringli_header.config = config_.dconfig;
  // Write header to ringli bitstream.
  // TODO(szabadka): Make it work for big-endian machines.
  ringli_data_.append(reinterpret_cast<char*>(&ringli_header),
                      sizeof(RingliHeader));
}

void StreamingRingliEncoder::CopyBlock(const uint8_t* data, size_t len,
                                       AudioBlock* block) {
  const size_t n_channels = format_.number_of_channels;
  std::vector<int16_t> block_data(kRingliBlockSize * n_channels);
  if (data) {
    memcpy(block_data.data(), data, len);
  }
  for (int j = 0; j < kRingliBlockSize; j++) {
    for (int channel = 0; channel < n_channels; channel++) {
      (*block)[channel][j] = block_data[j * n_channels + channel];
    }
  }
}

bool StreamingRingliEncoder::ProcessBlock(const RingliBlock& ringli_block) {
  if (config_.dconfig.use_predictive_coding) {
    return entropy_coder_->ProcessBlock(ringli_block, &ringli_data_);
  }
  ringli_headers_.push_back(ringli_block.header);
  ringli_blocks_.push_back(ringli_block);
  return true;
}

bool StreamingRingliEncoder::ProcessInput(const uint8_t* data, size_t len) {
  return wav_reader_.ProcessInput(data, len);
}

bool StreamingRingliEncoder::Flush() {
  if (config_.dconfig.use_predictive_coding) {
    return entropy_coder_->Flush(&ringli_data_);
  }
  CopyBlock(nullptr, 0, next_.get());
  ProcessBlock(EncodeWithDCT(config_, *prev_, *current_, *next_));
  CompressCoefficients(ringli_blocks_, format_.number_of_channels,
                       config_.dconfig.ecparams, &ringli_data_);
  return true;
}

size_t StreamingRingliEncoder::OutputSize() const {
  return ringli_data_.size() - output_pos_;
}

size_t StreamingRingliEncoder::CopyOutput(uint8_t* buffer, size_t len) {
  size_t nbytes = std::min(len, OutputSize());
  memcpy(buffer, reinterpret_cast<const uint8_t*>(&ringli_data_[output_pos_]),
         nbytes);
  output_pos_ += nbytes;
  return nbytes;
}

void RingliCompress(const std::string& wav_data,
                    const RingliEncoderConfig& config,
                    std::string* ringli_data) {
  StreamingRingliEncoder encoder(config);
  encoder.ProcessInput(reinterpret_cast<const uint8_t*>(wav_data.data()),
                       wav_data.size());
  encoder.Flush();
  ringli_data->resize(encoder.OutputSize());
  encoder.CopyOutput(const_cast<uint8_t*>(
                         reinterpret_cast<const uint8_t*>(ringli_data->data())),
                     ringli_data->size());
}

}  // namespace ringli
