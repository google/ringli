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

#ifndef COMMON_ERROR_NORM_H_
#define COMMON_ERROR_NORM_H_

#include <limits>
#include <string>
#include <vector>

namespace ringli {

struct ErrorNorm {
  double psnr = 0;
  double p4_diff = 0;
  double block_ends_diff = 0;
};

struct ChannelComparator {
  double p4_diff_sum = 0;
  double p2_diff_sum = 0;
  double p2_be_diff_sum = 0;
  double data_min = std::numeric_limits<double>::max();
  double data_max = std::numeric_limits<double>::min();
  double max_diff = 0;
  size_t max_diff_pos = 0;
  size_t num_samples = 0;
  double prev_a;
  double prev_b;

  void AddSamples(double a, double b);

  double PSNR() const;
  double P4Diff() const;
  double BlockEndsDiff() const;
};

class StreamComparator {
 public:
  explicit StreamComparator(size_t num_channels) : channel_cmp_(num_channels) {}

  void AddSamples(size_t c, double a, double b) {
    channel_cmp_[c].AddSamples(a, b);
  }

  double PSNR() const;
  double P4Diff() const;
  double BlockEndsDiff() const;
  double MaxDiff(size_t* channel, size_t* pos) const;

 private:
  std::vector<ChannelComparator> channel_cmp_;
};

ErrorNorm CompareFiles(const std::string& wav1, const std::string& wav2);

}  // namespace ringli

#endif  // COMMON_ERROR_NORM_H_
