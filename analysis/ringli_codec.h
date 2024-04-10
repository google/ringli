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

#ifndef ANALYSIS_RINGLI_CODEC_H_
#define ANALYSIS_RINGLI_CODEC_H_

#include <memory>
#include <string>

#include "analysis/audio_codec.h"
#include "common/streaming.h"
#include "decode/ringli_decoder.h"
#include "encode/ringli_encoder.h"

namespace ringli {

class StreamingRingliCodec : public StreamingAudioCodec {
 public:
  StreamingInterface* encoder() override {
    if (!encoder_) {
      encoder_ = std::make_unique<StreamingRingliEncoder>(config_);
    }
    return encoder_.get();
  }
  StreamingInterface* decoder() override {
    if (!decoder_) {
      decoder_ = std::make_unique<StreamingRingliDecoder>();
    }
    return decoder_.get();
  }

  std::string Name() const override { return "ringli"; }

  std::string ParamsToString() const override;

 protected:
  bool ParseParam(const std::string& param) override;

 private:
  RingliEncoderConfig config_;
  std::unique_ptr<StreamingRingliEncoder> encoder_;
  std::unique_ptr<StreamingRingliDecoder> decoder_;
};

}  // namespace ringli

#endif  // ANALYSIS_RINGLI_CODEC_H_
