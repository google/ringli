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

#ifndef ENCODE_NOISE_SHAPING_H_
#define ENCODE_NOISE_SHAPING_H_

#include <stddef.h>
#include <stdint.h>

#include <algorithm>
#include <vector>

#include "absl/log/check.h"
#include "common/data_defs/constants.h"

namespace ringli {

class NoiseShaper {
 public:
  NoiseShaper() : position_(0), noise_buffer_(kShapingFilterOrder, 0) {}

  void Reset() {
    position_ = 0;
    std::fill(noise_buffer_.begin(), noise_buffer_.end(), 0);
  }

  float GetFilteredNoise() const;

  void AddNewSample(float sample);

 private:
  uint32_t position_;
  std::vector<float> noise_buffer_;
};

}  // namespace ringli

#endif  // ENCODE_NOISE_SHAPING_H_
