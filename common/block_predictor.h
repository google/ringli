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

#ifndef COMMON_BLOCK_PREDICTOR_H_
#define COMMON_BLOCK_PREDICTOR_H_

#include <stddef.h>
#include <stdint.h>

#include <utility>
#include <vector>

#include "absl/log/check.h"
#include "common/covariance_lattice.h"
#include "common/predictor.h"
#include "common/ringli_header.h"

namespace ringli {
// Reconstructs the linear predictor coefficients from the given line spectral
// frequencies.
void ComputeLinearPredictorCoeffs(const uint16_t* quant_lsf, float* pcoefs,
                                  int order);

void ComputePredictorParams(RingliPredictiveHeader* header, float* pcoefs,
                            int order);

void DefaultLineSpectralFrequencies(uint16_t* quant_lsf, int order);

template <typename T>
class BlockPredictor : public Predictor {
 public:
  static BlockPredictor<T> CreateForEncoder(
      const T* block, int order, RingliPredictiveHeader* header,
      CovarianceLattice<int32_t>& covlattice_orig) {
    std::vector<float> pcoefs(order);
    covlattice_orig.FitPredictorCoeffs(&pcoefs[0], order);
    ComputePredictorParams(header, &pcoefs[0], order);
    return BlockPredictor(block, order, std::move(pcoefs));
  }

  static BlockPredictor<T> CreateForDecoder(
      const T* block, const RingliPredictiveHeader* header) {
    const int order = header->quant_lsf.size();
    std::vector<float> pcoefs(order);
    ComputeLinearPredictorCoeffs(&header->quant_lsf[0], &pcoefs[0], order);
    return BlockPredictor(block, order, std::move(pcoefs));
  }

  explicit BlockPredictor(const T* block, int order,
                          const std::vector<float>& pcoefs)
      : position_(0), block_(block), order_(order), pcoefs_(pcoefs) {}

  float Predict() override {
    if (position_ == 0) {
      return 0;
    }
    if (position_ == 1) {
      return block_[position_ - 1];
    }
    if (position_ < order_ && order_ > 2) {
      return 2 * block_[position_ - 1] - block_[position_ - 2];
    }
    float prediction = 0.0f;
    for (int p = 0; p < order_; ++p) {
      prediction += pcoefs_[p] * block_[position_ - 1 - p];
    }
    return prediction;
  };

  void AddNewSample(float sample) override {
    CHECK_EQ(block_[position_], sample);
    position_++;
  };
  void Reset() override {};

 private:
  uint32_t position_;
  const T* block_;
  const int order_;
  const std::vector<float> pcoefs_;
};

}  // namespace ringli

#endif  // COMMON_BLOCK_PREDICTOR_H_
