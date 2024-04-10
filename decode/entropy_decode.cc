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

#include "decode/entropy_decode.h"

#include <stddef.h>
#include <stdint.h>

#include <utility>
#include <vector>

#include "common/context.h"
#include "common/data_defs/constants.h"
#include "common/distributions.h"
#include "common/entropy_coding.h"
#include "common/ringli_header.h"
#include "decode/ans_decode.h"
#include "decode/arith_decode.h"
#include "decode/bit_reader.h"
#include "decode/context_map_decode.h"
#include "decode/ringli_input.h"

namespace ringli {

bool DecodeBase128(const uint8_t* data, const size_t len, size_t* pos,
                   size_t* val) {
  int shift = 0;
  uint64_t b;
  *val = 0;
  do {
    if (*pos >= len || shift > 57) {
      return false;
    }
    b = data[(*pos)++];
    *val |= (b & 0x7f) << shift;
    shift += 7;
  } while (b & 0x80);
  return true;
}

int16_t ConvertToSigned(uint16_t val) {
  int ux = val;
  if (ux & 1) {
    return -1 * ((ux + 1) / 2);
  } else {
    return ux / 2;
  }
}

bool DecodeDataLength(const uint8_t* data, const size_t len, size_t* pos,
                      size_t* data_len) {
  if (!DecodeBase128(data, len, pos, data_len)) {
    return false;
  }
  return *data_len <= len && *pos <= len - *data_len;
}

int DecodeSymbol(int alphabet_size, Prob* p, BinaryArithmeticDecoder* ac,
                 RingliInput* in) {
  int val0 = 0;
  int val1 = alphabet_size;
  while (val0 + 1 < val1) {
    int mid = (val0 + val1) >> 1;
    const int bit = ac->ReadBit(p[mid - 1].get_proba(), in);
    p[mid - 1].Add(bit);
    if (bit) {
      val0 = mid;
    } else {
      val1 = mid;
    }
  }
  return val0;
}

int DecodeValue(int symbol, int ndirect, RingliInput* in,
                BinaryArithmeticDecoder* ac = nullptr) {
  const int ndirect_symbols = 2 * ndirect - 1;
  if (symbol < ndirect_symbols) {
    return ConvertToSigned(symbol);
  } else {
    int sym = symbol - ndirect_symbols;
    const int sign = sym & 1;
    sym >>= 1;
    const int msb = sym & 1;
    const int nbits = sym >> 1;
    const int absval = ndirect - 2 + ((2 + msb) << nbits) +
                       (ac ? ac->ReadBits(nbits, in) : in->ReadBits(nbits));
    return (1 - 2 * sign) * absval;
  }
}

bool DecompressCoefficients(const char* input, size_t input_size,
                            size_t num_channels, size_t num_blocks,
                            const RingliDecoderConfig& config,
                            std::vector<RingliBlock>* ringli_blocks) {
  const uint8_t* data = reinterpret_cast<const uint8_t*>(input);
  const size_t num_contexts = 2 + kNumZeroDensityContexts;
  size_t pos = 0;
  std::vector<uint8_t> context_map;
  std::vector<ANSDecodingData> entropy_codes;
  if (!config.ecparams.arithmetic_only) {
    size_t histograms_size;
    if (!DecodeDataLength(data, input_size, &pos, &histograms_size) ||
        histograms_size == 0) {
      return false;
    }
    RingliBitReader br;
    RingliBitReaderInit(&br, &data[pos], histograms_size);
    int num_histograms;
    context_map.resize(num_contexts);
    if (!DecodeContextMap(num_contexts, &context_map[0], &num_histograms,
                          &br)) {
      return false;
    }
    entropy_codes.resize(num_histograms);
    for (int i = 0; i < num_histograms; ++i) {
      if (!entropy_codes[i].ReadFromBitStream(&br)) {
        return false;
      }
    }
    pos += histograms_size;
  }
  size_t coeff_data_size;
  if (!DecodeDataLength(data, input_size, &pos, &coeff_data_size) ||
      coeff_data_size == 0) {
    return false;
  }
  RingliInput in(&data[pos], coeff_data_size);
  BinaryArithmeticDecoder ac;
  ANSDecoder ans;
  if (!config.ecparams.arithmetic_only) {
    ans.Init(&in);
  }
  in.InitBitReader();
  ac.Init(&in);
  std::vector<Prob> last_nz_prob(kDctLength - 1);
  std::vector<Prob> is_zero_prob(kNumZeronessContexts);
  std::vector<Prob> sign_prob(kDctLength);
  std::vector<Prob> symbol_prob(num_contexts * (MAX_SYMBOLS - 1));

  size_t total_num_zeros = 0;
  size_t total_extra_bits = 0;
  for (size_t i = 0; i < num_blocks; ++i) {
    RingliBlock ringli_block(num_channels);
    for (int band = 0; band < kNumDctBands; ++band) {
      int quant_msb, quant_lsb;
      if (config.ecparams.arithmetic_only) {
        quant_msb = DecodeSymbol(MAX_SYMBOLS, &symbol_prob[0], &ac, &in);
        quant_lsb =
            DecodeSymbol(MAX_SYMBOLS, &symbol_prob[MAX_SYMBOLS - 1], &ac, &in);
      } else {
        quant_msb = ans.ReadSymbol(entropy_codes[context_map[0]], &in);
        quant_lsb = ans.ReadSymbol(entropy_codes[context_map[1]], &in);
      }
      ringli_block.header.dct.quant[band] = (quant_msb << 8) + quant_lsb;
    }
    for (size_t c = 0; c < num_channels; ++c) {
      for (size_t b = 0; b < kDctNumber; ++b) {
        auto& block = ringli_block.channels[c];
        int last_nz = DecodeSymbol(kDctLength, &last_nz_prob[0], &ac, &in);
        int num_nzeros = 0;
        for (int k = last_nz; k >= 0; --k) {
          int is_zero = 0;
          if (k == 0 || k < last_nz) {
            const int is_zero_ctx = ZeronessContext(num_nzeros, k);
            Prob* const p = &is_zero_prob[is_zero_ctx];
            is_zero = ac.ReadBit(p->get_proba(), &in);
            p->Add(is_zero);
          }
          total_num_zeros += is_zero;
          if (!is_zero) {
            int absval = 1;
            const int sign_ctx = k;
            Prob* const sign_p = &sign_prob[sign_ctx];
            int sign = ac.ReadBit(sign_p->get_proba(), &in);
            sign_p->Add(sign);
            const int absval_ctx = 2 + ZeroDensityContext(num_nzeros, k);
            int code = 0;
            if (config.ecparams.arithmetic_only) {
              code = DecodeSymbol(MAX_SYMBOLS,
                                  &symbol_prob[absval_ctx * (MAX_SYMBOLS - 1)],
                                  &ac, &in);
            } else {
              const int entropy_ix = context_map[absval_ctx];
              code = ans.ReadSymbol(entropy_codes[entropy_ix], &in);
            }
            if (code < NUM_DIRECT_CODES) {
              absval = code + 1;
            } else {
              int nbits = code - NUM_DIRECT_CODES + 1;
              int extra_bits_val = in.ReadBits(nbits);
              absval = NUM_DIRECT_CODES - 1 + (1 << nbits) + extra_bits_val;
              total_extra_bits += nbits;
            }
            block[b * kDctLength + k] = (1 - 2 * sign) * absval;
            ++num_nzeros;
          }
        }
      }
    }
    if (!in.ok()) {
      return false;
    }
    ringli_blocks->emplace_back(std::move(ringli_block));
  }
  if (!config.ecparams.arithmetic_only && !ans.CheckCRC()) {
    return false;
  }
  return true;
}

bool DecompressPredictiveRingliBlocks(const char* input, size_t input_size,
                                      size_t num_channels, size_t num_blocks,
                                      const RingliDecoderConfig& config,
                                      std::vector<RingliBlock>* ringli_blocks) {
  const uint8_t* data = reinterpret_cast<const uint8_t*>(input);
  size_t pos = 0;

  PredictiveContextModel context_model;
  const size_t num_contexts = 3 + kNumLSFContexts + context_model.NumContexts();
  std::vector<Prob> symbol_prob;
  if (config.ecparams.arithmetic_only) {
    symbol_prob.resize(num_contexts * (MAX_SYMBOLS - 1));
  }
  std::vector<uint8_t> context_map;
  std::vector<ANSDecodingData> entropy_codes;
  if (!config.ecparams.arithmetic_only) {
    size_t histograms_size;
    if (!DecodeDataLength(data, input_size, &pos, &histograms_size) ||
        histograms_size == 0) {
      return false;
    }
    RingliBitReader br;
    RingliBitReaderInit(&br, &data[pos], histograms_size);
    int num_histograms;
    context_map.resize(num_contexts);
    if (!DecodeContextMap(num_contexts, &context_map[0], &num_histograms,
                          &br)) {
      return false;
    }
    entropy_codes.resize(num_histograms);
    for (int i = 0; i < num_histograms; ++i) {
      if (!entropy_codes[i].ReadFromBitStream(&br)) {
        return false;
      }
    }
    pos += histograms_size;
  }

  RingliInput in(&data[pos], input_size - pos);
  BinaryArithmeticDecoder ac;
  ANSDecoder ans;
  if (!config.ecparams.arithmetic_only) {
    ans.Init(&in);
    in.InitBitReader();
  }
  ac.Init(&in);

  for (size_t bi = 0; bi < num_blocks; ++bi) {
    RingliBlock ringli_block(num_channels);
    for (size_t ci = 0; ci < num_channels; ++ci) {
      auto& header = ringli_block.header.pred[ci];
      if (!config.use_online_predictive_coding) {
        uint8_t order;
        if (config.ecparams.arithmetic_only) {
          order = DecodeSymbol(MAX_SYMBOLS, &symbol_prob[2 * (MAX_SYMBOLS - 1)],
                               &ac, &in);
        } else {
          order = ans.ReadSymbol(entropy_codes[context_map[2]], &in);
        }
        header.quant_lsf.resize(order);
        for (int p = 0; p < order; ++p) {
          const int pred_lsf = p * (kLSFQuant[p] / order);
          const int ctx = 3 + LSFContext(p, order);
          if (config.ecparams.arithmetic_only) {
            const int symbol = DecodeSymbol(
                MAX_SYMBOLS, &symbol_prob[ctx * (MAX_SYMBOLS - 1)], &ac, &in);
            header.quant_lsf[p] = pred_lsf + DecodeValue(symbol, 16, &in, &ac);
          } else {
            const int symbol =
                ans.ReadSymbol(entropy_codes[context_map[ctx]], &in);
            header.quant_lsf[p] = pred_lsf + DecodeValue(symbol, 16, &in);
          }
        }
      }
      auto& block = ringli_block.channels[ci];
      context_model.Reset();
      for (int i = 0; i < kRingliBlockSize; ++i) {
        const int ctx = 3 + kNumLSFContexts + context_model.Context();
        const int entropy_ix = context_map[ctx];
        int val;
        if (config.ecparams.arithmetic_only) {
          const int symbol = DecodeSymbol(
              MAX_SYMBOLS, &symbol_prob[ctx * (MAX_SYMBOLS - 1)], &ac, &in);
          val = DecodeValue(symbol, kPredNumDirectAbsval, &in, &ac);
        } else {
          const int symbol = ans.ReadSymbol(entropy_codes[entropy_ix], &in);
          val = DecodeValue(symbol, kPredNumDirectAbsval, &in);
        }
        context_model.Add(val);
        block[i] = val;
      }
    }
    ringli_blocks->emplace_back(std::move(ringli_block));
  }
  if (!config.ecparams.arithmetic_only && !ans.CheckCRC()) {
    return false;
  }
  return true;
}

IntegerArithmeticDecoder::IntegerArithmeticDecoder(int ndirect, int max_sym,
                                                   void* opaque,
                                                   ProcessOutput output_cb)
    : ndirect_absval_(ndirect),
      ndirect_symbols_(2 * ndirect - 1),
      max_symbols_(max_sym),
      opaque_(opaque),
      output_cb_(output_cb),
      state_(SYMBOL_DECODING),
      val0_(0),
      val1_(max_symbols_),
      distribution_(nullptr) {}

bool IntegerArithmeticDecoder::ProcessInput(uint16_t next_word) {
  ac_.Fill(next_word);
  while (ac_.HasBit()) {
    if (state_ == SYMBOL_DECODING) {
      if (!distribution_) {
        return false;
      }
      const int mid = (val0_ + val1_) >> 1;
      Prob& p = distribution_[mid - 1];
      const uint8_t prob = p.get_proba();
      const int bit = ac_.ReadBitNoFill(prob);
      p.Add(bit);
      if (bit) {
        val0_ = mid;
      } else {
        val1_ = mid;
      }
      if (val0_ + 1 == val1_) {
        if (val0_ < ndirect_symbols_) {
          if (!Output(ConvertToSigned(val0_))) {
            return false;
          }
        } else {
          val0_ -= ndirect_symbols_;
          sign_ = 1 - 2 * (val0_ & 1);
          val0_ >>= 1;
          msb_ = val0_ & 1;
          nbits_ = val0_ >> 1;
          if (nbits_ == 0) {
            if (!Output(sign_ * (ndirect_absval_ + msb_))) {
              return false;
            }
          } else {
            bitpos_ = 0;
            extra_bits_val_ = 0;
            state_ = EXTRA_BITS_DECODING;
          }
        }
      }
    } else if (state_ == EXTRA_BITS_DECODING) {
      extra_bits_val_ += ac_.ReadBitNoFill(128) << bitpos_;
      ++bitpos_;
      if (bitpos_ == nbits_) {
        const int absval =
            ndirect_absval_ - 2 + ((2 + msb_) << nbits_) + extra_bits_val_;
        if (!Output(sign_ * absval)) {
          return false;
        }
      }
    }
  }
  return true;
}

bool IntegerArithmeticDecoder::Output(int value) {
  state_ = SYMBOL_DECODING;
  val0_ = 0;
  val1_ = max_symbols_;
  return output_cb_(opaque_, value);
}

EntropyDecoder::EntropyDecoder(size_t num_channels, void* opaque,
                               ProcessSamples process_samples)
    : num_channels_(num_channels),
      opaque_(opaque),
      process_samples_(process_samples),
      int_decoder_(kPredNumDirectAbsval, MAX_SYMBOLS, this, ProcessOutputCb),
      context_model_(num_channels_),
      next_word_(0),
      input_shift_(0),
      samples_(num_channels),
      channel_idx_(0),
      idx_(0) {
  symbol_prob_.resize(context_model_[0].NumContexts() * (MAX_SYMBOLS - 1));
  SetContext();
}

bool EntropyDecoder::ProcessInput(const uint8_t* data, size_t len) {
  for (size_t i = 0; i < len; ++i) {
    next_word_ |= data[i] << input_shift_;
    input_shift_ += 8;
    if (input_shift_ == 16) {
      if (!int_decoder_.ProcessInput(next_word_)) {
        return false;
      }
      input_shift_ = 0;
      next_word_ = 0;
    }
  }
  return true;
}

bool EntropyDecoder::ProcessOutput(int value) {
  samples_[channel_idx_] = value;
  context_model_[channel_idx_].Add(value);
  ++channel_idx_;
  if (channel_idx_ == num_channels_) {
    channel_idx_ = 0;
    if (!process_samples_(opaque_, samples_.data())) {
      return false;
    }
    ++idx_;
    if (idx_ == kRingliBlockSize) {
      idx_ = 0;
      for (size_t ci = 0; ci < num_channels_; ++ci) {
        context_model_[ci].Reset();
      }
    }
  }
  SetContext();
  return true;
}

void EntropyDecoder::SetContext() {
  const int ctx = context_model_[channel_idx_].Context();
  int_decoder_.set_distribution(&symbol_prob_[ctx * (MAX_SYMBOLS - 1)]);
}

}  // namespace ringli
