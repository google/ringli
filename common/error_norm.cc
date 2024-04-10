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

#include "common/error_norm.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>

#include "absl/flags/flag.h"
#include "absl/log/check.h"
#include "common/data_defs/constants.h"
#include "common/logging.h"
#include "common/wav_header.h"
#include "common/wav_reader.h"

namespace ringli {

void ChannelComparator::AddSamples(double a, double b) {
  data_min = std::min(data_min, std::min(a, b));
  data_max = std::max(data_max, std::max(a, b));
  const double diff = std::abs(a - b);
  if (diff > max_diff) {
    max_diff = diff;
    max_diff_pos = num_samples;
  }
  const double diff2 = diff * diff;
  const double diff4 = diff2 * diff2;
  p2_diff_sum += diff2;
  p4_diff_sum += diff4;
  if (num_samples > 0 && num_samples % kRingliBlockSize == 0) {
    const double be_diff = (a - prev_a) - (b - prev_b);
    p2_be_diff_sum += be_diff * be_diff;
  }
  prev_a = a;
  prev_b = b;
  ++num_samples;
}

double ChannelComparator::PSNR() const {
  const double mse = p2_diff_sum / num_samples;
  const double range = data_max - data_min;
  return 10 * log(range * range / mse) / log(10);
}

double ChannelComparator::P4Diff() const {
  return pow(p4_diff_sum / num_samples, 1.0 / 4.0);
}

double ChannelComparator::BlockEndsDiff() const {
  const size_t num_blocks = (num_samples - 1) / kRingliBlockSize;
  return pow(p2_be_diff_sum / num_blocks, 0.5);
}

double StreamComparator::PSNR() const {
  double min_psnr = std::numeric_limits<double>::infinity();
  for (const auto& cmp : channel_cmp_) {
    min_psnr = std::min(min_psnr, cmp.PSNR());
  }
  return min_psnr;
}

double StreamComparator::P4Diff() const {
  double max_p4diff = 0.0;
  for (const auto& cmp : channel_cmp_) {
    max_p4diff = std::max(max_p4diff, cmp.P4Diff());
  }
  return max_p4diff;
}

double StreamComparator::BlockEndsDiff() const {
  double max_be_diff = 0.0;
  for (const auto& cmp : channel_cmp_) {
    max_be_diff = std::max(max_be_diff, cmp.BlockEndsDiff());
  }
  return max_be_diff;
}

double StreamComparator::MaxDiff(size_t* channel, size_t* pos) const {
  double max_diff = 0.0;
  for (size_t c = 0; c < channel_cmp_.size(); ++c) {
    if (c == 0 || channel_cmp_[c].max_diff > max_diff) {
      max_diff = channel_cmp_[c].max_diff;
      *channel = c;
      *pos = channel_cmp_[c].max_diff_pos;
    }
  }
  return max_diff;
}

ErrorNorm CompareFiles(const std::string& wav1, const std::string& wav2) {
  ErrorNorm error_norm;
  const uint8_t* data1 = reinterpret_cast<const uint8_t*>(wav1.data());
  const uint8_t* data2 = reinterpret_cast<const uint8_t*>(wav2.data());
  size_t pos1 = 0;
  size_t pos2 = 0;
  WavHeader header1;
  WavHeader header2;
  CHECK(ParseWavHeader(data1, wav1.size(), &pos1, &header1));
  CHECK(ParseWavHeader(data2, wav2.size(), &pos2, &header2));
  const int num_channels = header1.format_chunk.number_of_channels;
  CHECK_EQ(num_channels, header2.format_chunk.number_of_channels);
  uint32_t data_length1 = wav1.size() - pos1;
  uint32_t data_length2 = wav2.size() - pos2;
  std::pair<size_t, size_t> num_frames =
      std::minmax(data_length1 / num_channels / sizeof(int16_t),
                  data_length2 / num_channels / sizeof(int16_t));
  StreamComparator cmp(num_channels);
  for (size_t i = 0; i < num_frames.first; ++i) {
    for (int c = 0; c < num_channels; ++c) {
      const int16_t value1 = *reinterpret_cast<const int16_t*>(data1 + pos1);
      const int16_t value2 = *reinterpret_cast<const int16_t*>(data2 + pos2);
      pos1 += sizeof(value1);
      pos2 += sizeof(value2);
      cmp.AddSamples(c, value1, value2);
    }
  }
  const uint8_t* longer_data = data_length1 > data_length2 ? data1 : data2;
  size_t longer_pos = data_length1 > data_length2 ? pos1 : pos2;
  for (size_t i = num_frames.first; i < num_frames.second; ++i) {
    for (int c = 0; c < num_channels; ++c) {
      const int16_t value =
          *reinterpret_cast<const int16_t*>(longer_data + longer_pos);
      longer_pos += sizeof(value);
      cmp.AddSamples(c, value, 0);
    }
  }
  if (absl::GetFlag(FLAGS_log_level) >= 1) {
    size_t channel, pos;
    double max_diff = cmp.MaxDiff(&channel, &pos);
    printf("max diff: %f [channel: %zu block %zu idx %zu]\n", max_diff, channel,
           pos / kRingliBlockSize, pos % kRingliBlockSize);
  }
  error_norm.psnr = cmp.PSNR();
  error_norm.p4_diff = cmp.P4Diff();
  error_norm.block_ends_diff = cmp.BlockEndsDiff();
  return error_norm;
}

}  // namespace ringli
