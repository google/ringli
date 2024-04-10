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

#include "analysis/visqol.h"

#include "absl/log/check.h"
#include "analysis/resample.h"
#include "hwy/aligned_allocator.h"
#include "libsvm_nu_svr_model.h"
#include "visqol_api.h"

constexpr size_t SAMPLE_RATE = 48000;

namespace ringli {

ViSQOL::ViSQOL() {
  model_path_ = std::tmpnam(nullptr);
  std::ofstream output_stream(model_path_);
  CHECK(output_stream.good());
  output_stream.write(reinterpret_cast<char*>(visqol_model_bytes),
                      visqol_model_bytes_len);
  CHECK(output_stream.good());
  output_stream.close();
  CHECK(output_stream.good());
}

ViSQOL::~ViSQOL() { std::filesystem::remove(model_path_); }

float ViSQOL::MOS(const hwy::AlignedNDArray<float, 2>& reference,
                  const hwy::AlignedNDArray<float, 2>& degraded,
                  float sample_rate) const {
  CHECK_EQ(reference.shape()[0], degraded.shape()[0]);
  CHECK_EQ(reference.shape()[1], degraded.shape()[1]);
  hwy::AlignedNDArray<double, 2> resampled_reference =
      Resample<double>(reference, sample_rate, SAMPLE_RATE);
  hwy::AlignedNDArray<double, 2> resampled_degraded =
      Resample<double>(degraded, sample_rate, SAMPLE_RATE);

  Visqol::VisqolConfig config;
  config.mutable_options()->set_svr_model_path(model_path_);
  config.mutable_audio()->set_sample_rate(SAMPLE_RATE);

  // When running in audio mode, sample rates of 48k is recommended for
  // the input signals. Using non-48k input will very likely negatively
  // affect the comparison result. If, however, API users wish to run with
  // non-48k input, set this to true.
  config.mutable_options()->set_allow_unsupported_sample_rates(false);

  // ViSQOL will run in audio mode comparison by default.
  // If speech mode comparison is desired, set to true.
  config.mutable_options()->set_use_speech_scoring(false);

  // Speech mode will scale the MOS mapping by default. This means that a
  // perfect NSIM score of 1.0 will be mapped to a perfect MOS-LQO of 5.0.
  // Set to true to use unscaled speech mode. This means that a perfect
  // NSIM score will instead be mapped to a MOS-LQO of ~4.x.
  config.mutable_options()->set_use_unscaled_speech_mos_mapping(false);

  Visqol::VisqolApi visqol;
  CHECK_OK(visqol.Create(config));

  double sum_of_squares = 0.0;
  for (size_t channel_index = 0; channel_index < reference.shape()[0];
       ++channel_index) {
    const absl::Span<double> reference_span(
        resampled_reference[{channel_index}].data(),
        resampled_reference.shape()[1]);
    const absl::Span<double> degraded_span(
        resampled_degraded[{channel_index}].data(),
        resampled_degraded.shape()[1]);
    absl::StatusOr<Visqol::SimilarityResultMsg> comparison_status_or =
        visqol.Measure(reference_span, degraded_span);
    CHECK_OK(comparison_status_or);

    Visqol::SimilarityResultMsg similarity_result =
        comparison_status_or.value();

    const double moslqo = similarity_result.moslqo();
    sum_of_squares += moslqo * moslqo;
  }
  return std::sqrt(sum_of_squares / reference.shape()[0]);
}

}  // namespace ringli