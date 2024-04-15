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

#include "decode/ringli_decoder.h"

#include <stdio.h>
#include <string.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include "absl/log/check.h"
#include "common/block_predictor.h"
#include "common/convolve.h"
#include "common/data_defs/constants.h"
#include "common/data_defs/data_vector.h"
#include "common/dct.h"
#include "common/entropy_coding.h"
#include "common/fast_online_predictor.h"
#include "common/online_predictor.h"
#include "common/predictor.h"
#include "common/ringli_header.h"
#include "common/wav_header.h"
#include "common/wav_writer.h"
#include "decode/entropy_decode.h"

namespace ringli {
namespace {

std::unique_ptr<Predictor> GetPredictor(const RingliDecoderConfig& config,
                                        const RingliPredictiveHeader& header) {
  if (config.use_online_predictive_coding) {
    const bool fast_mode = config.predictor_fast_mode();
    if (fast_mode) {
      return std::make_unique<FastOnlinePredictor>();
    } else {
      return std::make_unique<OnlinePredictor>(kOnlinePredictorRegulariser);
    }
  } else {
    return BlockPredictor<kRingliBlockSize>::CreateForDecoder(&header);
  }
}

AudioBlock DecodePredictive(const RingliDecoderConfig& config,
                            const RingliBlock& encoded_block) {
  const size_t num_channels = encoded_block.channels.GetChannels().size();
  AudioBlock decoded_block(num_channels);
  for (size_t c = 0; c < num_channels; ++c) {
    const RingliPredictiveHeader& header = encoded_block.header.pred[c];
    const int quant = config.pred_quant;
    std::unique_ptr<Predictor> predictor = GetPredictor(config, header);
    for (int i = 0; i < kRingliBlockSize; i++) {
      float prediction = predictor->Predict();
      if (quant == 1) prediction = std::round(prediction);
      const float residual = quant * encoded_block.channels[c][i];
      const float sample_deq = prediction + residual;
      predictor->AddNewSample(sample_deq);
      decoded_block[c][i] = std::round(sample_deq);
    }
  }
  return decoded_block;
}

AudioBlock DecodeWithDCT(const RingliDecoderConfig& config,
                         const RingliBlock& prev, const RingliBlock& current,
                         const RingliBlock& next) {
  const size_t num_channels = current.channels.GetChannels().size();
  AudioBlock decoded_result(num_channels);
  for (size_t c = 0; c < num_channels; ++c) {
    DataVector<float, kACPredictionWindowSize> coeff_window;
    for (int i = 0; i < kACPredictionWindowSize; ++i) {
      const int k = i % kDctLength;
      if (i < kACPredictionBorder) {
        coeff_window[i] =
            prev.channels[c][kRingliBlockSize - kACPredictionBorder + i] *
            prev.header.dct.GetQuantizationCoef(k);
      } else if (i < kRingliBlockSize + kACPredictionBorder) {
        coeff_window[i] = current.channels[c][i - kACPredictionBorder] *
                          current.header.dct.GetQuantizationCoef(k);
      } else {
        coeff_window[i] =
            next.channels[c][i - kRingliBlockSize - kACPredictionBorder] *
            next.header.dct.GetQuantizationCoef(k);
      }
    }
    DataVector<float, kACPredictionWindowSize> output_window;
    for (int step = 0; step <= kNumACPredictionSteps; ++step) {
      int k_limit = step == kNumACPredictionSteps ? kDctLength
                                                  : kACPredictionStart << step;
      for (int i = 0; i < kACPredictionWindowSize; i += kDctLength) {
        DataVector<float, kDctLength> dct_data;
        for (int k = 0; k < k_limit; ++k) {
          dct_data[k] = coeff_window[i + k];
        }
        const DataVector<float, kDctLength> signal_data = InverseDCT(dct_data);
        for (int k = 0; k < kDctLength; ++k) {
          output_window[i + k] = signal_data[k];
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
        DataVector<float, kDctLength> dct_data = ForwardDCT(signal_data);
        for (int k = k_limit; k < 2 * k_limit; ++k) {
          coeff_window[i + k] += dct_data[k];
        }
      }
    }
    for (int i = 0; i < kRingliBlockSize; i++) {
      decoded_result[c][i] = std::round(output_window[i + kACPredictionBorder]);
    }
  }
  return decoded_result;
}

}  // namespace

void StreamingRingliDecoder::Reset() {
  ringli_data_.clear();
  wav_data_.clear();
  idx_ = 0;
  samples_written_ = 0;
  num_blocks_ = 0;
  input_pos_ = 0;
  output_pos_ = 0;
}

bool StreamingRingliDecoder::ProcessInput(const uint8_t* data, size_t len) {
  const size_t header_size = sizeof(RingliHeader);
  if (input_pos_ == 0 && len >= header_size) {
    if (!ProcessHeader(data, header_size)) {
      return false;
    }
    input_pos_ = header_size;
    data += header_size;
    len -= header_size;
  }
  if (input_pos_ >= header_size && len > 0) {
    if (ringli_header_.config.use_online_predictive_coding &&
        ringli_header_.config.ecparams.arithmetic_only) {
      entropy_decoder_->ProcessInput(data, len);
    } else {
      ringli_data_.append(reinterpret_cast<const char*>(data), len);
    }
    input_pos_ += len;
  }
  return true;
}

bool StreamingRingliDecoder::ProcessHeader(const uint8_t* data, size_t len) {
  CHECK_EQ(len, sizeof(ringli_header_));
  memcpy(&ringli_header_, data, len);
  if (memcmp(ringli_header_.ringli_id, kRingliId, sizeof(kRingliId)) != 0) {
    fprintf(stderr, "Unsupported ringli stream version\n");
    return false;
  }
  const size_t num_channels = ringli_header_.number_of_channels;
  const size_t bytes_per_sample = ringli_header_.bits_per_sample / 8;
  // Generate wav header based on ringli header.
  WavHeader wav_header;
  // RIFF chunk
  memcpy(wav_header.riff_chunk.riff_chunk_id, "RIFF", 4);
  wav_header.riff_chunk.riff_chunk_size = ringli_header_.data_length + 36;
  memcpy(wav_header.riff_chunk.wave_format, "WAVE", 4);
  // Format chunk
  memcpy(wav_header.format_chunk.format_chunk_id, "fmt ", 4);
  wav_header.format_chunk.format_chunk_size = 16;
  wav_header.format_chunk.audio_format = 1;
  wav_header.format_chunk.number_of_channels = num_channels;
  wav_header.format_chunk.sampling_frequency =
      ringli_header_.sampling_frequency;
  wav_header.format_chunk.byte_rate = num_channels *
                                      ringli_header_.sampling_frequency *
                                      ringli_header_.bits_per_sample / 8;
  wav_header.format_chunk.block_align = num_channels * 2;
  wav_header.format_chunk.bits_per_sample = ringli_header_.bits_per_sample;
  // Data chunk
  memcpy(wav_header.channel_id, "data", 4);
  wav_header.channel_data_length = ringli_header_.data_length;
  WriteWavHeader(wav_header, &wav_data_);
  wav_data_.reserve(wav_data_.size() + wav_header.channel_data_length);
  if (!ringli_header_.config.use_predictive_coding) {
    prev_ = std::make_unique<RingliBlock>(num_channels);
    current_ = std::make_unique<RingliBlock>(num_channels);
    next_ = std::make_unique<RingliBlock>(num_channels);
  } else if (ringli_header_.config.ecparams.arithmetic_only &&
             ringli_header_.config.use_online_predictive_coding) {
    entropy_decoder_ =
        std::make_unique<EntropyDecoder>(num_channels, this, ProcessSamplesCb);
    const bool fast_mode = ringli_header_.config.predictor_fast_mode();
    for (int i = 0; i < num_channels; ++i) {
      if (fast_mode) {
        predictors_.emplace_back(std::make_unique<FastOnlinePredictor>());
      } else {
        predictors_.emplace_back(
            std::make_unique<OnlinePredictor>(kOnlinePredictorRegulariser));
      }
    }

    noise_filters_.resize(num_channels);
    adaptive_quantizers_.resize(num_channels);
    for (size_t c = 0; c < num_channels; ++c) {
      predictors_[c]->Reset();
      noise_filters_[c].Reset();
      adaptive_quantizers_[c].Reset();
    }
  }
  remaining_samples_ = ringli_header_.data_length / bytes_per_sample;
  return true;
}

bool StreamingRingliDecoder::ProcessSamples(const int* samples) {
  // samples is a array of ints of size num_channels
  const size_t num_channels = ringli_header_.number_of_channels;
  std::vector<int32_t> decoded(num_channels);
  for (size_t c = 0; c < num_channels; ++c) {
    float quant;
    if (ringli_header_.config.use_adaptive_quantization) {
      quant = adaptive_quantizers_[c].QuantStep();
    } else {
      quant = ringli_header_.config.pred_quant;
    }
    const float prediction = predictors_[c]->Predict();
    const float residual = quant * samples[c];
    const float sample_deq = prediction + residual;
    predictors_[c]->AddNewSample(sample_deq);
    if (ringli_header_.config.use_adaptive_quantization) {
      adaptive_quantizers_[c].ProcessSample(sample_deq);
    }
    if (ringli_header_.config.use_noise_filter) {
      noise_filters_[c].AddNewSample(sample_deq);
      decoded[c] = std::round(noise_filters_[c].GetFilteredSample());
    } else {
      decoded[c] = std::round(sample_deq);
    }
  }

  if (!ringli_header_.config.use_noise_filter ||
      idx_ >= noise_filters_[0].get_delay()) {
    CHECK_LT(samples_written_, remaining_samples_);
    WriteSamples(decoded.data(), num_channels, &wav_data_);
    samples_written_ += num_channels;
  }
  ++idx_;
  if (idx_ % kRingliBlockSize == 0) {
    for (size_t ci = 0; ci < num_channels; ++ci) {
      predictors_[ci]->Reset();
      adaptive_quantizers_[ci].Reset();
    }
  }
  return true;
}

bool StreamingRingliDecoder::ProcessBlock(const RingliBlock& block) {
  if (ringli_header_.config.use_predictive_coding) {
    WriteBlock(DecodePredictive(ringli_header_.config, block));
  } else {
    if (num_blocks_ == 0) {
      *current_ = block;
    } else {
      *next_ = block;
      WriteBlock(
          DecodeWithDCT(ringli_header_.config, *prev_, *current_, *next_));
      *prev_ = *current_;
      *current_ = *next_;
    }
  }
  ++num_blocks_;
  return true;
}

void StreamingRingliDecoder::WriteBlock(const AudioBlock& block) {
  const size_t num_channels = ringli_header_.number_of_channels;
  const size_t samples_per_block = kRingliBlockSize * num_channels;
  const size_t num_samples = std::min(remaining_samples_, samples_per_block);
  WriteWavBlock(block, num_samples, &wav_data_);
  remaining_samples_ -= num_samples;
}

bool StreamingRingliDecoder::Flush() {
  if (ringli_header_.config.use_online_predictive_coding &&
      ringli_header_.config.ecparams.arithmetic_only) {
    if (ringli_header_.config.use_noise_filter) {
      const std::vector<int> flush_samples(ringli_header_.number_of_channels);
      for (int i = 0; i < noise_filters_[0].get_delay(); ++i) {
        ProcessSamples(flush_samples.data());
      }
    }
    CHECK_EQ(samples_written_, remaining_samples_);
    return true;
  }
  // Stream dimensions derived from header
  const size_t num_channels = ringli_header_.number_of_channels;
  const size_t bytes_per_sample = ringli_header_.bits_per_sample / 8;
  const size_t vec_size = kRingliBlockSize * bytes_per_sample;
  const size_t block_size = vec_size * num_channels;
  const size_t num_blocks =
      (ringli_header_.data_length + block_size - 1) / block_size;
  std::vector<RingliBlock> ringli_blocks;
  ringli_blocks.reserve(num_blocks);
  size_t ringli_pos = 0;
  const bool predictive_coding = ringli_header_.config.use_predictive_coding;
  if (predictive_coding) {
    if (!DecompressPredictiveRingliBlocks(
            &ringli_data_[ringli_pos], ringli_data_.size() - ringli_pos,
            num_channels, num_blocks, ringli_header_.config, &ringli_blocks)) {
      return false;
    }
  } else {
    if (!DecompressCoefficients(
            &ringli_data_[ringli_pos], ringli_data_.size() - ringli_pos,
            num_channels, num_blocks, ringli_header_.config, &ringli_blocks)) {
      return false;
    }
  }
  for (size_t i = 0; i < num_blocks; ++i) {
    ProcessBlock(ringli_blocks[i]);
  }
  if (!predictive_coding) {
    ProcessBlock(RingliBlock(num_channels));
  }
  return true;
}

size_t StreamingRingliDecoder::OutputSize() const {
  return wav_data_.size() - output_pos_;
}

size_t StreamingRingliDecoder::CopyOutput(uint8_t* buffer, size_t len) {
  size_t nbytes = std::min(len, OutputSize());
  memcpy(buffer, reinterpret_cast<const uint8_t*>(&wav_data_[output_pos_]),
         nbytes);
  output_pos_ += nbytes;
  return nbytes;
}

bool RingliDecompress(const std::string& ringli_data, std::string* wav_data) {
  StreamingRingliDecoder decoder;
  if (!decoder.ProcessInput(
          reinterpret_cast<const uint8_t*>(ringli_data.data()),
          ringli_data.size())) {
    return false;
  }
  if (!decoder.Flush()) {
    return false;
  }
  wav_data->resize(decoder.OutputSize());
  if (wav_data->size() !=
      decoder.CopyOutput(reinterpret_cast<uint8_t*>(wav_data->data()),
                         wav_data->size())) {
    return false;
  }
  return true;
}

}  // namespace ringli
