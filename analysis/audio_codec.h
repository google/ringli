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

#ifndef ANALYSIS_AUDIO_CODEC_H_
#define ANALYSIS_AUDIO_CODEC_H_

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "absl/log/log.h"
#include "absl/strings/str_cat.h"
#include "absl/types/span.h"
#include "common/streaming.h"

namespace ringli {

class AudioCodecBase {
 public:
  virtual ~AudioCodecBase() = default;

  bool ParseParams(absl::Span<const std::string> params) {
    if (params[0] != Name()) {
      return false;
    }
    // We skip params[0] because that is the type of the codec.
    for (size_t i = 1; i < params.size(); ++i) {
      if (!ParseParam(params[i])) {
        LOG(ERROR) << "Unknown param " << params[i];
        return false;
      }
    }
    return true;
  }

  virtual std::string Name() const = 0;
  virtual std::string ParamsToString() const = 0;

  std::string ToString() const {
    return absl::StrCat(Name(), ":", ParamsToString());
  }

 protected:
  virtual bool ParseParam(const std::string& param) = 0;
};

class AudioCodec : public AudioCodecBase {
 public:
  ~AudioCodec() override = default;

  virtual bool Compress(const std::string& wav_data,
                        std::string* compressed_data) = 0;
  virtual bool Decompress(const std::string& compressed_data,
                          std::string* wav_data) = 0;
};

class StreamingAudioCodec : public AudioCodec {
 public:
  ~StreamingAudioCodec() override = default;

  virtual StreamingInterface* encoder() = 0;
  virtual StreamingInterface* decoder() = 0;

  bool Compress(const std::string& wav_data,
                std::string* compressed_data) override {
    return Process(encoder(), wav_data, compressed_data);
  }
  bool Decompress(const std::string& compressed_data,
                  std::string* wav_data) override {
    return Process(decoder(), compressed_data, wav_data);
  }

 private:
  static bool Process(StreamingInterface* processor, const std::string& input,
                      std::string* output) {
    processor->Reset();
    if (!processor->ProcessInput(reinterpret_cast<const uint8_t*>(input.data()),
                                 input.size()) ||
        !processor->Flush()) {
      return false;
    }
    output->resize(processor->OutputSize());
    processor->CopyOutput(reinterpret_cast<uint8_t*>(output->data()),
                          output->size());
    return true;
  }
};

}  // namespace ringli

#endif  // ANALYSIS_AUDIO_CODEC_H_
