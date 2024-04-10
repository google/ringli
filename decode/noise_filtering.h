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

#ifndef DECODE_NOISE_FILTERING_H_
#define DECODE_NOISE_FILTERING_H_

#include <stddef.h>
#include <stdint.h>

#include <algorithm>
#include <vector>

#include "absl/log/check.h"
#include "common/data_defs/constants.h"

namespace ringli {

class SymNoiseFilter {
 public:
  SymNoiseFilter() : position_(0), sample_buffer_(kNoiseFilterOrder, 0) {
    CHECK_EQ(kNoiseFilterOrder % 2, 1);
  }
  void Reset() {
    position_ = 0;
    std::fill(sample_buffer_.begin(), sample_buffer_.end(), 0);
  }

  float GetFilteredSample() const;

  void AddNewSample(float sample);

  int get_delay() const { return SymNoiseFilter::kDelay; }

 private:
  static constexpr int kDelay = (kNoiseFilterOrder - 1) / 2;
  uint32_t position_;
  std::vector<float> sample_buffer_;
};

}  // namespace ringli

#endif  // DECODE_NOISE_FILTERING_H_
