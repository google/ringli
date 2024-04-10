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

#ifndef COMMON_WAV_HEADER_H_
#define COMMON_WAV_HEADER_H_

#include <stddef.h>
#include <stdint.h>

namespace ringli {

// http://soundfile.sapp.org/doc/WaveFormat/
struct RiffChunk {
  char riff_chunk_id[4] = "   ";
  uint32_t riff_chunk_size;
  char wave_format[4] = "   ";
} __attribute__((packed));

struct FormatChunk {
  char format_chunk_id[4] = "   ";
  uint32_t format_chunk_size;
  uint16_t audio_format;  // Audio format 1=PCM,6=mulaw,7=alaw,     257=IBM
                          // Mu-Law, 258=IBM A-Law, 259=ADPCM
  uint16_t number_of_channels;  // Number of channels 1=Mono 2=Sterio
  uint32_t sampling_frequency;  // Sampling Frequency in Hz
  uint32_t byte_rate;           // bytes per second
  uint16_t block_align;         // 2=16-bit mono, 4=16-bit stereo
  uint16_t bits_per_sample;
} __attribute__((packed));

struct WavHeader {
  uint16_t NumChannels() const { return format_chunk.number_of_channels; }
  size_t NumSamples() const {
    return channel_data_length / (format_chunk.number_of_channels *
                                  (format_chunk.bits_per_sample / 8));
  }
  uint32_t SampleRate() const { return format_chunk.sampling_frequency; }

  RiffChunk riff_chunk;
  FormatChunk format_chunk;
  size_t header_size;

  char channel_id[4] = "   ";
  uint32_t channel_data_length;

  int LengthSeconds() const {
    return channel_data_length * 1.0 / format_chunk.byte_rate;
  }
} __attribute__((packed));

// Minimum size of WAV header, this is what the ringli decoder produces.
constexpr size_t kWavHeaderSize = 44;

}  // namespace ringli

#endif  // COMMON_WAV_HEADER_H_
