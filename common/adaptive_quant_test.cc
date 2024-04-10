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

#include <cmath>
#include <cstdint>
#include <vector>

#include "absl/log/check.h"
#include "common/data_defs/constants.h"
#include "gtest/gtest.h"

namespace ringli {

TEST(MultiEMA, EMAInputOnes) {
  const int num_channels = 4;
  const std::vector<float> sig(num_channels, 1);

  const float momentum = 0.5;
  const int num_steps = 20;

  MultiEMA ema_filter(num_channels, momentum);

  for (int i = 1; i <= num_steps; i++) {
    ema_filter.ProcessSample(sig);
    const std::vector<float>& out = ema_filter.Buffer();
    for (int j = 0; j < num_channels; j++) {
      EXPECT_EQ(out[j], 1 - std::pow(momentum, i));
    }
  }
}

TEST(MultiEMA, ResetTest) {
  const int num_channels = 2;
  const std::vector<float> sig(num_channels, 1);

  const float momentum = 0.5;

  MultiEMA ema_filter(num_channels, momentum);

  ema_filter.ProcessSample(sig);
  const std::vector<float> out0 = ema_filter.Buffer();
  ema_filter.ProcessSample(sig);
  const std::vector<float> out1 = ema_filter.Buffer();
  ema_filter.Reset();
  ema_filter.ProcessSample(sig);
  const std::vector<float> out2 = ema_filter.Buffer();
  for (int i = 0; i < num_channels; i++) {
    EXPECT_NE(out0[i], out1[i]);
    EXPECT_EQ(out0[i], out2[i]);
  }
}

TEST(IirMultiFilter, HighpassLowpassFIR) {
  const std::vector<uint16_t> sig(30, 1);

  FilterCoeffArray bwd_coeffs = {{-1, 1, 0}, {0.5, 0.5, 0}};

  FilterCoeffArray fwd_coeffs = {{1, 0, 0}, {1, 0, 0}};

  IirMultiFilter filter = IirMultiFilter(bwd_coeffs, fwd_coeffs);
  std::vector<float> out(kNumIIRFilters, 0);

  // we skip first entry since filter has 1 sample delay
  filter.ProcessSample(sig[0], &out[0]);

  for (int i = 1; i < sig.size() - 1; i++) {
    filter.ProcessSample(sig[i], &out[0]);
    EXPECT_EQ(0, out[0]);  // highpass returns 0
    EXPECT_EQ(1, out[1]);  // lowpass returns 1
  }
}

TEST(IirMultiFilter, HighpassLowpassIIR) {
  const std::vector<uint16_t> sig(30, 1);

  FilterCoeffArray bwd_coeffs = {{-1, 1, 0}, {0.5, 0.5, 0}};

  FilterCoeffArray fwd_coeffs = {{1, -0.5, 0}, {1, -0.5, 0}};

  IirMultiFilter filter = IirMultiFilter(bwd_coeffs, fwd_coeffs);
  std::vector<float> out(kNumIIRFilters, 0);

  for (int i = 0; i < sig.size() - 1; i++) {
    filter.ProcessSample(sig[i], &out[0]);
  }
  // iir highpass converges to 0
  EXPECT_NEAR(0, out[0], 0.0001);
  // iir lowpass converges to 2
  EXPECT_NEAR(2, out[1], 0.0001);
}

TEST(IirMultiFilter, ResetTest) {
  IirMultiFilter filter = IirMultiFilter(kBwdCoeffs, kFwdCoeffs);

  const float sig = 3.5f;
  std::vector<float> out0(4);
  filter.ProcessSample(sig, &out0[0]);
  std::vector<float> out1(4);
  filter.ProcessSample(sig, &out1[0]);
  filter.Reset();
  std::vector<float> out2(4);
  filter.ProcessSample(sig, &out2[0]);

  for (int i = 0; i < kNumIIRFilters; i++) {
    EXPECT_EQ(out0[i], out2[i]);
    EXPECT_NE(out0[i], out1[i]);
  }
}

TEST(AdaptiveQuantizer, FrameEdgeScalingTest) {
  AdaptiveQuantizer adaptive_quantizer;

  const float sig_value = 1.0f;

  const std::vector<float> sig(kRingliBlockSize, sig_value);

  for (int i = 0; i < kRingliBlockSize; i++) {
    float scale = adaptive_quantizer.FrameEdgeScaling(sig[i]);

    adaptive_quantizer.ProcessSample(sig[i]);  // increase idx_

    if (i < kAdaptiveQuantRampUpSteps) {
      CHECK_LT(scale, sig_value);
    } else if (i >= (kRingliBlockSize - kAdaptiveQuantRampDownSteps)) {
      CHECK_LT(scale, sig_value);
    } else {
      CHECK_EQ(scale, sig_value);
    }
  }
}

TEST(AdaptiveQuantizer, MinMaxTest) {
  AdaptiveQuantizer adaptive_quantizer;

  const std::vector<float> zero_sig(10, 0);
  const std::vector<float> large_sig = {-10, 10, 20, 30, 40, 50,
                                        60,  70, 80, 90, 100};
  EXPECT_LE(kQuantMin, kQuantStart);
  for (int i = 0; i < 10; i++) {
    float quant_step = adaptive_quantizer.QuantStep();
    adaptive_quantizer.ProcessSample(zero_sig[i]);
    EXPECT_GE(quant_step, kQuantMin);
  }

  adaptive_quantizer.Reset();

  for (int i = 0; i < 10; i++) {
    float quant_step = adaptive_quantizer.QuantStep();
    adaptive_quantizer.ProcessSample(large_sig[i] * 2000);
    EXPECT_LE(quant_step, kQuantMax);
  }
}

TEST(AdaptiveQuantizer, ResetTest) {
  AdaptiveQuantizer adaptive_quantizer;

  const std::vector<float> sig = {-10, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};

  for (int i = 0; i < 10; i++) {
    float quant_step = adaptive_quantizer.QuantStep();
    adaptive_quantizer.ProcessSample(sig[i] * 1000);
    if (i < kHQSteps) {  // at the beginning we use max resolution
      EXPECT_EQ(quant_step, kQuantStart);
    } else {  // we use adaptive mode after some samples
      EXPECT_NE(quant_step, kQuantStart);
    }
  }

  adaptive_quantizer.Reset();

  for (int i = 0; i < 10; i++) {
    float quant_step = adaptive_quantizer.QuantStep();
    adaptive_quantizer.ProcessSample(sig[i] * 1000);
    if (i < kHQSteps) {
      EXPECT_EQ(quant_step, kQuantStart);
    } else {
      EXPECT_NE(quant_step, kQuantStart);
    }
  }
}

}  // namespace ringli
