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

#ifndef COMMON_WAV_WRITER_H_
#define COMMON_WAV_WRITER_H_

#include <ostream>
#include <string>

#include "common/data_defs/constants.h"
#include "common/wav_header.h"

namespace ringli {

void WriteWavHeader(const WavHeader& wav_header, std::string* output);

void WriteWavBlock(const AudioBlock& block, size_t max_samples,
                   std::string* output);

void WriteSamples(const int16_t* samples, size_t num_samples,
                  std::string* output);

void WriteSamples(const int32_t* samples, size_t num_samples,
                  std::string* output);

}  // namespace ringli

#endif  // COMMON_WAV_WRITER_H_
