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

#include "common/adaptive_quant.h"

#include <stddef.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "absl/types/span.h"
#include "common/data_defs/constants.h"

namespace ringli {

void MultiEMA::ProcessSample(absl::Span<const float> samples) {
  for (int i = 0; i < num_channels_; ++i) {
    buffer_[i] = momentum_ * buffer_[i] + (1 - momentum_) * samples[i];
  }
}

void IirMultiFilter::ProcessSample(float sample, float out[kNumIIRFilters]) {
  backward_buffer_[idx_ % kBackwardBufferSize] = sample;

  for (int j = 0; j < kNumIIRFilters; ++j) {  // Loop over filters
    float val = 0.0f;
    int buffer_index = idx_ % kBackwardBufferSize;
    for (int k = 0; k < kBackwardBufferSize; ++k) {
      val += backward_buffer_[buffer_index] * bwd_coeffs_[j][k];
      buffer_index--;
      if (buffer_index < 0) buffer_index += kBackwardBufferSize;
    }

    buffer_index = idx_ % kForwardBufferSize;
    for (int m = 0; m < kForwardBufferSize; ++m) {
      val -= forward_buffers_[j][buffer_index] * fwd_coeffs_[j][m + 1];
      buffer_index--;
      if (buffer_index < 0) buffer_index += kForwardBufferSize;
    }
    val *= rescale_coeffs_[j];

    out[j] = val;
    forward_buffers_[j][idx_ % kForwardBufferSize] = val;
  }
  idx_++;
}

float NoisePowerToQuantStep(float noise_power) {
  return sqrtf(noise_power * 12);
}

float AdaptiveQuantizer::FrameEdgeScaling(float sample) const {
  // downscales first and last samples in ringli block to smooth out artefacts
  // introduced by IIR filter warmup
  constexpr static int kRampDownThreshold =
      kRingliBlockSize - kAdaptiveQuantRampDownSteps - 1;
  constexpr static float kAdaptiveQuantRampDownStepsPlusOne =
      kAdaptiveQuantRampDownSteps + 1.0f;
  constexpr static float kAdaptiveQuantRampUpStepsPlusOne =
      kAdaptiveQuantRampUpSteps + 1.0f;

  if (idx_ < kAdaptiveQuantRampUpSteps) {
    const float scaling = (idx_ + 1.0f) / kAdaptiveQuantRampUpStepsPlusOne;
    return sample * scaling;
  } else if (idx_ >= kRampDownThreshold) {
    const float scaling =
        (kRingliBlockSize - idx_) / kAdaptiveQuantRampDownStepsPlusOne;
    return sample * scaling;
  } else {
    return sample;
  }
}

void AdaptiveQuantizer::ProcessSample(const float sample) {
  // IIR filters separate signal into bands of interest
  iir_filters_.ProcessSample(FrameEdgeScaling(sample), filter_outputs_);

  // Square filter outputs to get instant band power estimates
  for (int i = 0; i < kNumIIRFilters; ++i) {
    filter_outputs_[i] *= filter_outputs_[i];
  }

  // EMA estimates power averaged over time
  ema_filters_.ProcessSample(filter_outputs_);
  const std::vector<float>& band_power = ema_filters_.Buffer();

  // scale power into max noise estimate
  max_noise_power_ = std::numeric_limits<float>::max();
  for (int i = 0; i < kNumIIRFilters; ++i) {
    max_noise_power_ =
        std::min(max_noise_power_, band_power[i] * kMaskingFilterTargetNSR[i]);
  }
  idx_++;
}

float AdaptiveQuantizer::QuantStep() const {
  if (idx_ < kHQSteps) {
    return kQuantStart;
  } else {
    return std::clamp(NoisePowerToQuantStep(max_noise_power_), kQuantMin,
                      kQuantMax);
  }
}

}  // namespace ringli
