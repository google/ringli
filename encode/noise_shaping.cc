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

#include "encode/noise_shaping.h"

#include <stddef.h>

#include <vector>

#include "common/data_defs/constants.h"

namespace ringli {

float NoiseShaper::GetFilteredNoise() const {
  float output = 0;
  if (position_ == 0) {
    // do nothing
  } else if (position_ < kShapingFilterOrder) {
    // use order 1 filter
    output += kShapingO1Coeffs[0] *
              noise_buffer_[(position_ - 1) % kShapingFilterOrder];
  } else {  // full order filter
    for (int i = 1; i <= kShapingFilterOrder; ++i) {
      output += noise_buffer_[(position_ - i) % kShapingFilterOrder] *
                kShapingCoeffs[i - 1];
    }
  }
  return output;
}

void NoiseShaper::AddNewSample(float noise_sample) {
  noise_buffer_[position_ % kShapingFilterOrder] = noise_sample;
  position_++;
}

}  // namespace ringli
