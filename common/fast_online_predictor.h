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

#ifndef COMMON_FAST_ONLINE_PREDICTOR_H_
#define COMMON_FAST_ONLINE_PREDICTOR_H_

#include <stddef.h>
#include <stdint.h>

#include <cstring>

#include "common/data_defs/constants.h"
#include "common/predictor.h"

namespace ringli {

class FastOnlinePredictor : public Predictor {
 public:
  FastOnlinePredictor() { FastOnlinePredictor::Reset(); }

  void Reset() override {
    position_ = 0;
    memset(data_buffer_, 0, sizeof(data_buffer_));
    memset(k_, 0, sizeof(k_));
    memset(g_, 0, sizeof(g_));
    memset(d_, 0, sizeof(d_));
  }

  float Predict() override;

  // Implementation of the adaptive lattice method for updating the reflection
  // coefficients based on the following paper:
  // J. Makhoul and R. Viswanathan, "Adaptive lattice methods for linear
  // prediction", Proc. ICASSP-78, vol. 3, pp. 83-86.
  void AddNewSample(float sample) override;

 private:
  uint32_t position_ = 0;
  float data_buffer_[kOnlinePredictorBufferSize] = {0};
  // The variable names for reflection coefficients, forward and backward
  // prediction errors, etc. follow the notations in the paper.
  float k_[kOnlinePredictorOrder + 1] = {0};
  float g_[kOnlinePredictorOrder + 1] = {0};
  float d_[kOnlinePredictorOrder + 1] = {0};
  const float beta_ = 0.999f;
  const float regul_ = 1.0f;
};

}  // namespace ringli

#endif  // COMMON_FAST_ONLINE_PREDICTOR_H_
