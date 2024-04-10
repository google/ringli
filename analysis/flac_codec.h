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

#ifndef ANALYSIS_FLAC_CODEC_H_
#define ANALYSIS_FLAC_CODEC_H_

#include <memory>
#include <string>

#include "absl/strings/substitute.h"
#include "analysis/audio_codec.h"
#include "common/streaming.h"

namespace ringli {

class StreamingFlacEncoder : public StreamingInterface {
 public:
  explicit StreamingFlacEncoder(int compression_level)
      : compression_level_(compression_level) {
    Reset();
  }

  void Reset() override;
  bool ProcessInput(const uint8_t* data, size_t len) override;
  bool Flush() override;
  size_t OutputSize() const override;
  size_t CopyOutput(uint8_t* buffer, size_t len) override;

 private:
  int compression_level_;
  std::string wav_data_;
  std::string flac_data_;
  size_t output_pos_;
};

class StreamingFlacDecoder : public StreamingInterface {
 public:
  StreamingFlacDecoder() { Reset(); }

  void Reset() override;
  bool ProcessInput(const uint8_t* data, size_t len) override;
  bool Flush() override;
  size_t OutputSize() const override;
  size_t CopyOutput(uint8_t* buffer, size_t len) override;

 private:
  std::string flac_data_;
  std::string wav_data_;
  size_t output_pos_;
};

class StreamingFlacCodec : public StreamingAudioCodec {
 public:
  StreamingInterface* encoder() override {
    if (!encoder_) {
      encoder_ = std::make_unique<StreamingFlacEncoder>(compression_level_);
    }
    return encoder_.get();
  }
  StreamingInterface* decoder() override {
    if (!decoder_) {
      decoder_ = std::make_unique<StreamingFlacDecoder>();
    }
    return decoder_.get();
  }
  std::string Name() const override { return "flac"; }

 protected:
  bool ParseParam(const std::string& param) override {
    if (param[0] == 'c') {
      compression_level_ = std::stoi(param.substr(1));
      return true;
    }
    return false;
  }

  std::string ParamsToString() const override {
    return absl::Substitute("c$0", compression_level_);
  }

 private:
  int compression_level_ = 8;
  std::unique_ptr<StreamingFlacEncoder> encoder_;
  std::unique_ptr<StreamingFlacDecoder> decoder_;
};

}  // namespace ringli

#endif  // ANALYSIS_FLAC_CODEC_H_
