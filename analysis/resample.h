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

#ifndef ANALYSIS_RESAMPLE_H_
#define ANALYSIS_RESAMPLE_H_

#include "hwy/aligned_allocator.h"
#include "samplerate.h"

namespace ringli {

template <typename O, typename I>
hwy::AlignedNDArray<O, 2> Convert(const hwy::AlignedNDArray<I, 2>& input) {
  if constexpr (std::is_same<O, I>::value) {
    hwy::AlignedNDArray<O, 2> result(input.shape());
    hwy::CopyBytes(input.data(), result.data(),
                   input.memory_size() * sizeof(I));
    return result;
  }
  hwy::AlignedNDArray<O, 2> output(input.shape());
  for (size_t channel_index = 0; channel_index < input.shape()[0];
       ++channel_index) {
    O* output_ptr = output[{channel_index}].data();
    const I* input_ptr = input[{channel_index}].data();
    for (size_t sample_index = 0; sample_index < input.shape()[1];
         ++sample_index) {
      *(output_ptr + sample_index) = *(input_ptr + sample_index);
    }
  }
  return output;
}

template <typename O, typename I>
hwy::AlignedNDArray<O, 2> Resample(const hwy::AlignedNDArray<I, 2>& samples,
                                   float in_sample_rate,
                                   float out_sample_rate) {
  if (in_sample_rate == out_sample_rate) {
    return Convert<O>(samples);
  }

  const hwy::AlignedNDArray<float, 2> samples_as_floats =
      Convert<float>(samples);
  hwy::AlignedNDArray<float, 2> result_as_floats(
      {samples.shape()[0],
       static_cast<size_t>(samples.size() * out_sample_rate / in_sample_rate)});
  for (size_t channel_index = 0; channel_index < samples.shape()[0];
       ++channel_index) {
    SRC_DATA resample_data = {
        .data_in = samples_as_floats[{channel_index}].data(),
        .data_out = result_as_floats[{channel_index}].data(),
        .input_frames = static_cast<long>(samples_as_floats.shape()[1]),
        .output_frames = static_cast<long>(result_as_floats.shape()[1]),
        .src_ratio = out_sample_rate / in_sample_rate,
    };
    CHECK_EQ(src_simple(&resample_data, SRC_SINC_BEST_QUALITY, 1), 0);
  }
  return Convert<O>(result_as_floats);
}

}  // namespace ringli

#endif  // ANALYSIS_RESAMPLE_H_
