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

#include "common/block_predictor.h"

#include <stdint.h>
#include <stdlib.h>

#include <cmath>

#include "Eigen/Dense"
#include "common/data_defs/constants.h"
#include "common/ringli_header.h"
#include "gtest/gtest.h"

namespace ringli {

TEST(BlockPredictorTest, DefaultLineSpectralFrequencies) {
  constexpr int kOrder = 8;
  uint16_t quant_lsf[kOrder];
  DefaultLineSpectralFrequencies(quant_lsf, kOrder);
  float pcoefs[kOrder];
  ComputeLinearPredictorCoeffs(quant_lsf, pcoefs, kOrder);
  EXPECT_NEAR(pcoefs[0], 1.0, 1e-6);
  for (int i = 1; i < kOrder; ++i) {
    EXPECT_NEAR(pcoefs[i], 0.0, 1e-6);
  }
}

TEST(BlockPredictorTest, Roundtrip) {
  // We start from line spectral frequencies because those always define a
  // stable predictor.
  constexpr int kOrder = 8;
  float lsf[kOrder];
  lsf[0] = 0.1;
  for (int i = 1; i < kOrder; ++i) {
    lsf[i] = lsf[i - 1] + 0.1 * (1.0 - lsf[i - 1]);
  }

  RingliPredictiveHeader header;

  uint16_t quant_lsf[kOrder];
  for (int i = 0; i < kOrder; ++i) {
    quant_lsf[i] = std::round(lsf[i] * kLSFQuant[i]);
  }
  float pcoefs[kOrder];
  ComputeLinearPredictorCoeffs(quant_lsf, pcoefs, kOrder);
  ComputePredictorParams(&header, pcoefs, kOrder);
  for (int i = 0; i < kOrder; ++i) {
    EXPECT_EQ(quant_lsf[i], header.quant_lsf[i]);
  }
}

}  // namespace ringli
