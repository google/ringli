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

#include "decode/noise_filtering.h"

#include <stddef.h>

#include <vector>

#include "absl/log/check.h"
#include "common/data_defs/constants.h"

namespace ringli {

float SymNoiseFilter::GetFilteredSample() const {
  // Folded FIR filter implementation
  CHECK_EQ(SymNoiseFilter::kDelay, (kNoiseFilterOrder - 1) / 2);
  float output = 0;

  for (int i = 1; i <= SymNoiseFilter::kDelay; ++i) {
    output += kNoiseFilterCoeffs[i - 1] *
              (sample_buffer_[(position_ - i) % kNoiseFilterOrder] +
               sample_buffer_[(position_ - kNoiseFilterOrder - 1 + i) %
                              kNoiseFilterOrder]);
  }
  output += kNoiseFilterCoeffs[SymNoiseFilter::kDelay] *
            sample_buffer_[(position_ - SymNoiseFilter::kDelay - 1) %
                           kNoiseFilterOrder];

  return output;
}

void SymNoiseFilter::AddNewSample(float sample) {
  sample_buffer_[position_ % kNoiseFilterOrder] = sample;
  position_++;
}

}  // namespace ringli
