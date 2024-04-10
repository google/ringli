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

#ifndef COMMON_RINGLI_HEADER_H_
#define COMMON_RINGLI_HEADER_H_

#include <stddef.h>

#include <cstdint>
#include <vector>

#include "common/data_defs/constants.h"
#include "common/data_defs/data_vector.h"
#include "common/entropy_coding.h"

namespace ringli {

struct RingliDecoderConfig {
  EntropyCodingParams ecparams;
  uint8_t use_predictive_coding = false;
  uint8_t use_online_predictive_coding = false;
  uint16_t pred_quant = 1;
  uint8_t effort = 5;

  bool use_noise_filter = false;
  bool use_adaptive_quantization = false;

  bool predictor_fast_mode() const { return effort <= 5; }
} __attribute__((packed));

struct RingliHeader {
  char ringli_id[8];
  uint32_t header_length;
  uint16_t number_of_channels;
  uint32_t sampling_frequency;
  uint16_t bits_per_sample;
  uint32_t data_length;
  RingliDecoderConfig config;
} __attribute__((packed));

struct RingliDCTHeader {
  uint16_t quant[kNumDctBands] = {};

  uint16_t GetQuantizationCoef(int index) const {
    for (int b = 0; b < kNumDctBands; ++b) {
      if (index < kDctBandEnd[b]) {
        return quant[b];
      }
    }
    return 1;  // should not happen
  }
} __attribute__((packed));

struct RingliPredictiveHeader {
  // Line spectral frequencies normalized in the [0, 1] interval and quantized
  // with a variable-precision quantizer: the pth coefficient is quantized to
  // round(kLSFQuant[p]) + 1 levels.
  std::vector<uint16_t> quant_lsf;
};

struct RingliBlockHeader {
  explicit RingliBlockHeader(size_t num_channesl) : pred(num_channesl) {}
  RingliDCTHeader dct;
  std::vector<RingliPredictiveHeader> pred;
};

struct RingliBlock {
  explicit RingliBlock(size_t num_channels)
      : header(num_channels), channels(num_channels) {}

  RingliBlockHeader header;
  DataVectorPack<int32_t, kRingliBlockSize> channels;
};

}  // namespace ringli

#endif  // COMMON_RINGLI_HEADER_H_
