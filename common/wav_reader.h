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

#ifndef COMMON_WAV_READER_H_
#define COMMON_WAV_READER_H_

#include <istream>
#include <map>
#include <string>
#include <vector>

#include "common/data_defs/constants.h"
#include "common/wav_header.h"

namespace ringli {

// A streaming wav reader which calls user registered callbacks for chunk data
// processing.
class StreamingWavReader {
 public:
  StreamingWavReader();

  typedef bool (*ProcessDataCallback)(void* opaque, const uint8_t* data,
                                      size_t len, size_t chunk_pos,
                                      size_t chunk_size);
  // Register a processing callback for a specific 4-byte chunk id.
  // opaque will be provided as first argument to callback.
  // If buffer_size is nonzero, this many bytes of chunk data will be buffered
  // before calling callback, except at the end of the chunk where the
  // callback will be called with the rest of the buffered chunk data, if any.
  void RegisterCallback(const char* id, void* opaque,
                        ProcessDataCallback callback, size_t buffer_size);

  void Reset();

  bool ProcessInput(const uint8_t* data, size_t len);

 private:
  size_t stream_pos_ = 0;
  size_t stream_size_ = 0;
  uint32_t chunk_start_ = 0;
  char chunk_id_[4];
  uint32_t chunk_size_;
  size_t max_buffer_size_ = 0;
  std::vector<uint8_t> chunk_buffer_;
  size_t buffer_pos_;

  struct CallbackConfig {
    void* opaque;
    ProcessDataCallback callback = nullptr;
    size_t buffer_size = 0;
  };
  std::map<std::string, CallbackConfig> callbacks_;
  CallbackConfig current_cb_;
};

void DescribeWavHeader(const WavHeader& header);

bool ParseFormatChunk(const uint8_t* data, size_t len, size_t* pos,
                      FormatChunk* format);

bool ParseWavHeader(const uint8_t* data, size_t len, size_t* pos,
                    WavHeader* header);

void VerifyWavHeader(const WavHeader& header, const uint8_t* data, size_t len);

}  // namespace ringli

#endif  // COMMON_WAV_READER_H_
