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

#include "common/fast_online_predictor.h"

#include "absl/log/check.h"
#include "common/covariance_lattice.h"
#include "common/data_defs/constants.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "common/fast_online_predictor.cc"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace ringli {
namespace HWY_NAMESPACE {

using DF = HWY_CAPPED(float, kOnlinePredictorOrder);
constexpr DF df;

float PredictFromBackwardErrors(const float* __restrict k,
                                const float* __restrict g, int order) {
  auto vkg = Zero(df);
  for (int i = 0; i < order; i += Lanes(df)) {
    const auto vk = LoadU(df, k + i);
    const auto vg = LoadU(df, g + i);
    vkg = MulAdd(vk, vg, vkg);
  }
  return -1.0f * GetLane(SumOfLanes(df, vkg));
}

void UpdatePredictorState(float sample, float* __restrict k,
                          float* __restrict f, float* __restrict g,
                          float* __restrict d, int order, float beta,
                          float regul) {
  DCHECK_EQ(order % Lanes(df), 0);
  f[0] = sample;
  for (int m = 0; m < order; ++m) {
    f[m + 1] = f[m] + k[m] * g[m];
  }
  const auto vbeta = Set(df, beta);
  const auto vregul = Set(df, regul);
  for (int i = order - Lanes(df); i >= 0; i -= Lanes(df)) {
    const auto vf = LoadU(df, f + i);
    const auto vf1 = LoadU(df, f + i + 1);
    const auto vg = LoadU(df, g + i);
    const auto vk = LoadU(df, k + i);
    const auto vff = Mul(vf, vf);
    const auto vgg = Mul(vg, vg);
    const auto vd = Add(MulAdd(LoadU(df, d + i), vbeta, vregul), Add(vff, vgg));
    const auto vg1 = MulAdd(vk, vf, vg);
    StoreU(Sub(vk, Div(Add(Mul(vf, vg1), Mul(vf1, vg)), vd)), df, k + i);
    StoreU(vg1, df, g + i + 1);
    StoreU(vd, df, d + i);
  }
  g[0] = sample;
}

// NOLINTNEXTLINE(google-readability-namespace-comments)
}  // namespace HWY_NAMESPACE
}  // namespace ringli
HWY_AFTER_NAMESPACE();

#if HWY_ONCE

namespace ringli {

HWY_EXPORT(PredictFromBackwardErrors);
HWY_EXPORT(UpdatePredictorState);

float FastOnlinePredictor::Predict() {
  float prediction = 0;
  if (position_ == 0) {
    // do nothing
  } else if (position_ == 1) {
    // predict first sample
    prediction = data_buffer_[0];
  } else if (position_ < 2 * kOnlinePredictorOrder) {
    // use precomputed optimal order 2 predictor
    prediction =
        kOnlineO2PredictorDefaultCoeffs[0] * data_buffer_[position_ - 1] -
        kOnlineO2PredictorDefaultCoeffs[1] * data_buffer_[position_ - 2];
  } else {
    prediction = HWY_DYNAMIC_DISPATCH(PredictFromBackwardErrors)(
        &k_[1], &g_[0], kOnlinePredictorOrder);
  }
  return prediction;
}

void FastOnlinePredictor::AddNewSample(float sample) {
  data_buffer_[position_ % kOnlinePredictorBufferSize] = sample;
  position_++;

  if (position_ < 2 * kOnlinePredictorOrder) {
    return;
  }
  float f[kOnlinePredictorOrder + 1];
  if (position_ == 2 * kOnlinePredictorOrder) {
    // initialize
    CovarianceLattice<float> covlattice(&data_buffer_[0], position_,
                                        kOnlinePredictorOrder, regul_);
    covlattice.FitReflectionCoeffs(&k_[0], kOnlinePredictorOrder);
    for (int i = 0; i < position_; ++i) {
      f[0] = data_buffer_[i];
      for (int m = 1; m <= kOnlinePredictorOrder; ++m) {
        f[m] = f[m - 1] + k_[m] * g_[m - 1];
      }
      for (int m = kOnlinePredictorOrder; m >= 1; --m) {
        g_[m] = k_[m] * f[m - 1] + g_[m - 1];
      }
      g_[0] = data_buffer_[i];
    }
  } else {
    HWY_DYNAMIC_DISPATCH(UpdatePredictorState)
    (sample, &k_[1], &f[0], &g_[0], &d_[1], kOnlinePredictorOrder, beta_,
     regul_);
  }
}
}  // namespace ringli
#endif  // HWY_ONCE
