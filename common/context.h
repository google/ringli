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

#ifndef COMMON_CONTEXT_H_
#define COMMON_CONTEXT_H_

#include <stdint.h>

#include <algorithm>
#include <cstdlib>

#include "absl/log/check.h"
#include "common/data_defs/constants.h"
#include "common/log2floor.h"

namespace ringli {

static const uint8_t kNonzeroBuckets[kDctLength] = {
    0,  1,  2,  3,  4,  4,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,
    7,  7,  7,  7,  7,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,
    9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  10, 10, 10,
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
};
static const int kNumNonzeroBuckets = 11;

static const int kNumZeronessContexts = kNumNonzeroBuckets * kDctLength;

inline int ZeronessContext(int nonzeros_left, int k) {
  DCHECK_LT(nonzeros_left, kDctLength);
  return kNonzeroBuckets[nonzeros_left] * kDctLength + k;
}

static const uint8_t kFreqContext[kDctLength] = {
    0,  1,  2,  3,  4,  4,  5,  5,  6,  6,  7,  7,  8,  8,  8,  8,
    9,  9,  9,  9,  10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12,
    13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
};

static const uint16_t kNumNonzeroContext[kDctLength] = {
    0,   16,  16,  32,  32,  32,  48,  48,  48,  48,  64,  64,  64,
    64,  64,  64,  80,  80,  80,  80,  80,  80,  80,  80,  95,  95,
    95,  95,  95,  95,  95,  95,  95,  95,  95,  95,  95,  95,  95,
    95,  109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109,
    109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109};

static const int kNumZeroDensityContexts = 120;

inline int ZeroDensityContext(int nonzeros_left, int k) {
  DCHECK_LT(nonzeros_left, kDctLength);
  DCHECK_LT(k, kDctLength);
  return kNumNonzeroContext[nonzeros_left] + kFreqContext[k];
}

static const size_t kNumLSFContexts =
    ((kMaxPredictorOrder + 2) * kMaxPredictorOrder) / 4;

inline int LSFContext(int p, int order) {
  return ((order - 2) * order) / 4 + p;
}

class PredictiveContextModel {
 public:
  PredictiveContextModel() { Init(); }

  void Init() {
    context_mul_[kOrder - 1] = 1;
    for (int i = kOrder - 2; i >= 0; --i) {
      const int next_ncontexts = 2 * nbits_context_map_[i + 1][kMaxNumBits] + 1;
      context_mul_[i] = context_mul_[i + 1] * next_ncontexts;
    }
    Reset();
  }

  void Reset() {
    pos_ = 0;
    for (int i = 0; i < kOrder; ++i) {
      nbits_[i] = 0;
      sign_[i] = 0;
    }
  }

  void Add(int value) {
    nbits_[pos_ % kOrder] =
        value == 0
            ? 0
            : std::min(kMaxNumBits, 1 + Log2FloorNonZero(std::abs(value)));
    sign_[pos_ % kOrder] = value >= 0 ? 0 : 1;
    ++pos_;
  }

  int Context() const {
    int ctx = 0;
    for (int i = 0; i < kOrder; ++i) {
      const int idx = (pos_ + kOrder - 1 - i) % kOrder;
      const int nbits_ctx = nbits_context_map_[i][nbits_[idx]];
      if (nbits_ctx > 0) {
        ctx += (2 * nbits_ctx - sign_[idx]) * context_mul_[i];
      }
    }
    return ctx;
  }

  int NumContexts() const {
    int num = 1;
    for (int i = 0; i < kOrder; ++i) {
      num *= (2 * nbits_context_map_[i][kMaxNumBits] + 1);
    }
    return num;
  }

 private:
  static constexpr int kOrder = 3;
  static constexpr int kMaxNumBits = 11;
  const int nbits_context_map_[kOrder][kMaxNumBits + 1] = {
      {0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6},
      {0, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3},
      {0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2},
  };
  int context_mul_[kOrder];

  uint32_t pos_;
  int nbits_[kOrder];
  int sign_[kOrder];
};

}  // namespace ringli

#endif  // COMMON_CONTEXT_H_
