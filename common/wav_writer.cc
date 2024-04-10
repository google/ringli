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

#include "common/wav_writer.h"

#include <stddef.h>
#include <stdint.h>

#include <algorithm>
#include <limits>
#include <string>

#include "common/data_defs/constants.h"
#include "common/wav_header.h"

namespace ringli {

namespace {

template <typename T>
void Write(const T& t, std::string* output) {
  const size_t pos = output->size();
  output->resize(pos + sizeof(T));
  memcpy(&(*output)[pos], &t, sizeof(T));
}

void Write(const char id[4], std::string* output) {
  output->append(std::string(id, id + 4));
}

}  // namespace

void WriteWavHeader(const WavHeader& wav_header, std::string* output) {
  Write(wav_header.riff_chunk, output);
  const FormatChunk& format = wav_header.format_chunk;
  Write(format.format_chunk_id, output);
  Write(format.format_chunk_size, output);
  Write(format.audio_format, output);
  Write(format.number_of_channels, output);
  Write(format.sampling_frequency, output);
  Write(format.byte_rate, output);
  Write(format.block_align, output);
  Write(format.bits_per_sample, output);
  Write(wav_header.channel_id, output);
  Write(wav_header.channel_data_length, output);
}

void WriteWavBlock(const AudioBlock& block, size_t max_samples,
                   std::string* output) {
  size_t max_time_slots =
      std::min(kRingliBlockSize, max_samples / block.GetChannels().size());
  // TODO(szabadka) Make this work with other wav formats as well.
  output->reserve(output->size() + max_samples * sizeof(int16_t));
  const int32_t minval = std::numeric_limits<int16_t>::min();
  const int32_t maxval = std::numeric_limits<int16_t>::max();
  for (int i = 0; i < max_time_slots; i++) {
    for (const RingliVector& raw_block : block.GetChannels()) {
      const int32_t raw_value = raw_block[i];
      const int16_t clamped_value = std::clamp(raw_value, minval, maxval);
      Write(clamped_value, output);
    }
  }
}

void WriteSamples(const int16_t* samples, size_t num_samples,
                  std::string* output) {
  for (size_t i = 0; i < num_samples; ++i) {
    Write(samples[i], output);
  }
}

void WriteSamples(const int32_t* samples, size_t num_samples,
                  std::string* output) {
  const int32_t minval = std::numeric_limits<int16_t>::min();
  const int32_t maxval = std::numeric_limits<int16_t>::max();
  for (size_t i = 0; i < num_samples; ++i) {
    const int16_t value = std::clamp(samples[i], minval, maxval);
    Write(value, output);
  }
}

}  // namespace ringli
