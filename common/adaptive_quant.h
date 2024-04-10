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

#ifndef COMMON_ADAPTIVE_QUANT_H_
#define COMMON_ADAPTIVE_QUANT_H_

#include <stddef.h>

#include <algorithm>
#include <vector>

#include "absl/types/span.h"
#include "common/data_defs/constants.h"

namespace ringli {

class MultiEMA {
 public:
  MultiEMA(const int num_channels, const float momentum)
      : num_channels_(num_channels), momentum_(momentum) {
    buffer_.resize(num_channels_);
  }

  void Reset() { std::fill(buffer_.begin(), buffer_.end(), 0); }

  void ProcessSample(absl::Span<const float> samples);

  const std::vector<float>& Buffer() const { return buffer_; }

 private:
  int num_channels_;
  float momentum_;
  std::vector<float> buffer_;
};

class IirMultiFilter {
 public:
  IirMultiFilter(const FilterCoeffArray bwd_coeffs,
                 const FilterCoeffArray fwd_coeffs)
      : idx_(0) {
    Reset();
    for (int i = 0; i < kNumIIRFilters; ++i) {
      rescale_coeffs_[i] = 1 / fwd_coeffs[i][0];
      for (int j = 0; j < kIIROrder; ++j) {
        bwd_coeffs_[i][j] = bwd_coeffs[i][j];
        fwd_coeffs_[i][j] = fwd_coeffs[i][j];
      }
    }
  }

  void ProcessSample(float sample, float out[kNumIIRFilters]);

  void Reset() {
    memset(backward_buffer_, 0, sizeof(backward_buffer_));
    memset(forward_buffers_, 0, sizeof(forward_buffers_));

    idx_ = 0;
  }

 private:
  constexpr static size_t kBackwardBufferSize = kIIROrder;
  constexpr static size_t kForwardBufferSize = kIIROrder - 1;

  FilterCoeffArray bwd_coeffs_;
  FilterCoeffArray fwd_coeffs_;
  int idx_;
  float backward_buffer_[kBackwardBufferSize];
  float forward_buffers_[kNumIIRFilters][kForwardBufferSize];
  float rescale_coeffs_[kNumIIRFilters];
};

float NoisePowerToQuantStep(float noise_power);

class AdaptiveQuantizer {
 public:
  AdaptiveQuantizer()
      : iir_filters_(kBwdCoeffs, kFwdCoeffs),
        ema_filters_(kNumIIRFilters, kEmaMomentum),
        idx_(0) {}

  void ProcessSample(float sample);

  float QuantStep() const;

  void Reset() {
    idx_ = 0;
    iir_filters_.Reset();
    ema_filters_.Reset();
    max_noise_power_ = 0;
  }

  float FrameEdgeScaling(float sample) const;

 private:
  float filter_outputs_[kNumIIRFilters] = {0};
  IirMultiFilter iir_filters_;
  MultiEMA ema_filters_;
  int idx_;
  float max_noise_power_ = 0;
};

}  // namespace ringli

#endif  // COMMON_ADAPTIVE_QUANT_H_
