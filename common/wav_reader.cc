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

#include "common/wav_reader.h"

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include <algorithm>
#include <iostream>
#include <istream>
#include <ostream>
#include <string>

#include "absl/log/check.h"
#include "common/wav_header.h"

namespace ringli {

namespace {

std::string to_string_4(const char word[4]) {
  return std::string(word, word + 4);
}

template <typename T>
bool Read(const uint8_t* data, size_t len, size_t* pos, T* t) {
  if (*pos + sizeof(T) > len) return false;
  memcpy(t, data + (*pos), sizeof(T));
  *pos += sizeof(T);
  return true;
}

bool Read(const uint8_t* data, size_t len, size_t* pos, char* id) {
  if (*pos + 4 > len) return false;
  memcpy(id, data + (*pos), 4);
  *pos += 4;
  return true;
}

// Verifies that the RIFF container has a wav stream inside.
bool VerifyRIFFChunk(void* opaque, const uint8_t* data, size_t len,
                     size_t chunk_pos, size_t chunk_size) {
  return len == 4 && chunk_size == 4 && memcmp("WAVE", data, len) == 0;
}

}  // namespace

void DescribeWavHeader(const WavHeader& header) {
  std::cout << "================ R I F F =================" << std::endl;
  std::cout << "Riff chunk id:\t\t"
            << to_string_4(header.riff_chunk.riff_chunk_id) << std::endl;
  std::cout << "Riff chunk size:\t" << header.riff_chunk.riff_chunk_size
            << std::endl;
  std::cout << "Riff format:\t\t" << to_string_4(header.riff_chunk.wave_format)
            << std::endl;
  std::cout << "============= F O R M A T ================" << std::endl;
  std::cout << "Format chunk id:\t"
            << to_string_4(header.format_chunk.format_chunk_id) << std::endl;
  std::cout << "Format chunk size:\t" << header.format_chunk.format_chunk_size
            << std::endl;
  std::cout << "Number of channels:\t" << header.format_chunk.number_of_channels
            << std::endl;
  std::cout << "Audio format:\t\t" << header.format_chunk.audio_format
            << std::endl;
  std::cout << "Sample rate:\t\t" << header.format_chunk.sampling_frequency
            << std::endl;
  std::cout << "Byte rate:\t\t" << header.format_chunk.byte_rate << std::endl;
  std::cout << "Block align:\t\t" << header.format_chunk.block_align
            << std::endl;
  std::cout << "Bits per sample:\t" << header.format_chunk.bits_per_sample
            << std::endl;
  std::cout << "=============== D A T A ==================" << std::endl;
  std::cout << "Channel id:\t" << to_string_4(header.channel_id) << std::endl;
  std::cout << "Data length:\t" << header.channel_data_length << std::endl;
  std::cout << "Track length:\t" << header.LengthSeconds() << std::endl;
}

bool ParseFormatChunk(const uint8_t* data, size_t len, size_t* pos,
                      FormatChunk* format) {
  uint16_t audio_format;
  uint16_t number_of_channels;
  uint32_t sampling_frequency;
  uint32_t byte_rate;
  uint16_t block_align;
  uint16_t bits_per_sample;
  bool ok =
      (Read(data, len, pos, &audio_format) &&
       Read(data, len, pos, &number_of_channels) &&
       Read(data, len, pos, &sampling_frequency) &&
       Read(data, len, pos, &byte_rate) && Read(data, len, pos, &block_align) &&
       Read(data, len, pos, &bits_per_sample));
  if (ok) {
    format->audio_format = audio_format;
    format->number_of_channels = number_of_channels;
    format->sampling_frequency = sampling_frequency;
    format->byte_rate = byte_rate;
    format->block_align = block_align;
    format->bits_per_sample = bits_per_sample;
  }
  return ok;
}

bool ParseWavHeader(const uint8_t* data, size_t len, size_t* pos,
                    WavHeader* header) {
  RiffChunk riff_chunk;
  FormatChunk format_chunk;
  char format_chunk_id[4];
  uint32_t format_chunk_size;
  if (!Read(data, len, pos, &riff_chunk) ||
      !Read(data, len, pos, format_chunk_id) ||
      !Read(data, len, pos, &format_chunk_size) ||
      !ParseFormatChunk(data, len, pos, &header->format_chunk)) {
    return false;
  }
  header->riff_chunk = riff_chunk;
  memcpy(header->format_chunk.format_chunk_id, format_chunk_id, 4);
  header->format_chunk.format_chunk_size = format_chunk_size;

  bool found_data_chunk = false;
  while (!found_data_chunk) {
    char channel_id[4] = "   ";
    uint32_t channel_data_length;
    if (!Read(data, len, pos, channel_id) ||
        !Read(data, len, pos, &channel_data_length)) {
      return false;
    }
    memcpy(header->channel_id, channel_id, 4);
    header->channel_data_length = channel_data_length;
    if (to_string_4(header->channel_id) == "data") {
      found_data_chunk = true;
    } else {
      header->riff_chunk.riff_chunk_size -= header->channel_data_length + 8;
      *pos += header->channel_data_length;
    }
  }
  header->header_size = *pos;
  return true;
}

void VerifyWavHeader(const WavHeader& header, const uint8_t* data, size_t len) {
  size_t pos = 0;
  RiffChunk riff;
  CHECK(Read(data, len, &pos, &riff));
  FormatChunk format;
  CHECK(Read(data, len, &pos, format.format_chunk_id));
  uint32_t format_chunk_size;
  CHECK(Read(data, len, &pos, &format_chunk_size));
  format.format_chunk_size = format_chunk_size;
  CHECK(ParseFormatChunk(data, len, &pos, &format));
  char channel_id[4];
  uint32_t channel_data_len;
  CHECK(Read(data, len, &pos, channel_id));
  CHECK(Read(data, len, &pos, &channel_data_len));
  // Don't compare riff_chunk_size because it can be different if input had
  // some metadata chunks.
  CHECK_EQ(memcmp(&header.riff_chunk, &riff, sizeof(riff.riff_chunk_id)), 0);
  CHECK_EQ(memcmp(&header.format_chunk, &format, sizeof(format)), 0);
  CHECK_EQ(memcmp(header.channel_id, channel_id, 4), 0);
  CHECK_EQ(header.channel_data_length, channel_data_len);
}

StreamingWavReader::StreamingWavReader() {
  RegisterCallback("RIFF", nullptr, VerifyRIFFChunk, 4);
}

void StreamingWavReader::Reset() {
  stream_pos_ = 0;
  stream_size_ = 0;
  chunk_start_ = 0;
}

bool StreamingWavReader::ProcessInput(const uint8_t* data, size_t len) {
  size_t pos = 0;
  while (pos < len) {
    const size_t chunk_data_start = chunk_start_ + 8;
    if (stream_pos_ < chunk_data_start) {
      // Parse the 8-byte chunk header (chunk id and chunk size).
      const size_t n =
          std::min<size_t>(len - pos, chunk_data_start - stream_pos_);
      memcpy(&chunk_id_[stream_pos_ - chunk_start_], &data[pos], n);
      stream_pos_ += n;
      pos += n;
      // RIFF chunk size contains the size of the whole stream, not the
      // size of the 4-byte format string, so we save the stream size and
      // set the chunk size to 4, so that the parser can skip to the next
      // chunk.
      if (stream_pos_ == chunk_data_start &&
          memcmp(chunk_id_, "RIFF", 4) == 0) {
        stream_size_ = chunk_size_ + 8;
        chunk_size_ = 4;
      }
      continue;
    }
    const size_t chunk_pos = stream_pos_ - chunk_data_start;
    if (chunk_pos == 0) {
      // We are the start of the chunk data, so we look up the appropriate
      // chunk processing callback and initialize the chunk data buffer.
      auto it = callbacks_.find(to_string_4(chunk_id_));
      if (it != callbacks_.end()) {
        current_cb_ = it->second;
      } else {
        current_cb_.callback = nullptr;
      }
      buffer_pos_ = 0;
    }
    if (chunk_pos < chunk_size_) {
      size_t n = std::min<size_t>(len - pos, chunk_size_ - chunk_pos);
      if (current_cb_.callback) {
        // If we have a callback registered for this chunk, then we
        // either buffer as many bytes as the callback registerer
        // requested, or call the callback on the input buffer directly.
        const bool buffered = current_cb_.buffer_size > 0;
        if (buffered) {
          CHECK_LT(buffer_pos_, current_cb_.buffer_size);
          n = std::min<size_t>(n, current_cb_.buffer_size - buffer_pos_);
          memcpy(&chunk_buffer_[buffer_pos_], &data[pos], n);
          buffer_pos_ += n;
        }
        if (buffer_pos_ == current_cb_.buffer_size ||
            chunk_pos + n == chunk_size_) {
          const uint8_t* buffer = buffered ? &chunk_buffer_[0] : &data[pos];
          const size_t buflen = buffered ? buffer_pos_ : n;
          if (!current_cb_.callback(current_cb_.opaque, buffer, buflen,
                                    chunk_pos + n - buflen, chunk_size_)) {
            return false;
          }
          buffer_pos_ = 0;
        }
      }
      stream_pos_ += n;
      pos += n;
    } else {
      // Start next chunk.
      chunk_start_ = stream_pos_;
    }
  }
  return true;
}

void StreamingWavReader::RegisterCallback(const char* id, void* opaque,
                                          ProcessDataCallback callback,
                                          size_t buffer_size) {
  CallbackConfig config = {opaque, callback, buffer_size};
  callbacks_[std::string(id)] = config;
  if (buffer_size > max_buffer_size_) {
    chunk_buffer_.resize(buffer_size);
    max_buffer_size_ = buffer_size;
  }
}

}  // namespace ringli
