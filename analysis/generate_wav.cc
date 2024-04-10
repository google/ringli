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

#include "analysis/generate_wav.h"

#include <cmath>

#include "common/wav_header.h"
#include "common/wav_writer.h"

namespace ringli {

std::string GenerateWav(const std::vector<Waveform>& waveforms,
                        float sample_rate, float seconds, float noise) {
  const size_t num_samples = static_cast<size_t>(sample_rate * seconds);
  const float period = 1.0f / sample_rate;
  std::vector<int16_t> samples(num_samples);
  std::srand(0);
  const float rand_max_reciprocal =
      2 * std::numeric_limits<int16_t>::max() / static_cast<float>(RAND_MAX);
  for (size_t i = 0; i < num_samples; ++i) {
    const float t = i * period;
    samples[i] = noise * (std::rand() - RAND_MAX * 0.5) * rand_max_reciprocal;
    for (const auto& waveform : waveforms) {
      samples[i] +=
          std::sin(t * 2 * M_PI * waveform.frequency + waveform.phase) *
          (std::numeric_limits<int16_t>::max() * waveform.amplitude);
    }
  }
  const size_t data_size = samples.size() * sizeof(int16_t);
  const size_t num_channels = 1;
  std::string wav_data;
  WavHeader wav_header;
  // RIFF chunk
  memcpy(wav_header.riff_chunk.riff_chunk_id, "RIFF", 4);
  wav_header.riff_chunk.riff_chunk_size = data_size + kWavHeaderSize - 8;
  memcpy(wav_header.riff_chunk.wave_format, "WAVE", 4);
  // Format chunk
  memcpy(wav_header.format_chunk.format_chunk_id, "fmt ", 4);
  wav_header.format_chunk.format_chunk_size = 16;
  wav_header.format_chunk.audio_format = 1;
  wav_header.format_chunk.number_of_channels = 1;
  wav_header.format_chunk.sampling_frequency = sample_rate;
  wav_header.format_chunk.byte_rate =
      num_channels * sample_rate * sizeof(int16_t);
  wav_header.format_chunk.block_align = num_channels * 2;
  wav_header.format_chunk.bits_per_sample = sizeof(int16_t) * 8;
  // Data chunk
  memcpy(wav_header.channel_id, "data", 4);
  wav_header.channel_data_length = data_size;
  WriteWavHeader(wav_header, &wav_data);
  wav_data.reserve(wav_data.size() + wav_header.channel_data_length);
  WriteSamples(samples.data(), samples.size(), &wav_data);
  return wav_data;
}

}  // namespace ringli