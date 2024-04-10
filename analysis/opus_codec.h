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

#ifndef ANALYSIS_OPUS_CODEC_H_
#define ANALYSIS_OPUS_CODEC_H_

#include <memory>
#include <string>
#include <vector>

#include "absl/strings/substitute.h"
#include "analysis/audio_codec.h"
#include "common/streaming.h"
#include "common/wav_header.h"
#include "common/wav_reader.h"
#include "opus_multistream.h"

namespace ringli {

class StreamingOpusEncoder : public StreamingInterface {
 public:
  explicit StreamingOpusEncoder(int bitrate);

  void Reset() override;
  bool ProcessInput(const uint8_t* data, size_t len) override;
  bool Flush() override;
  size_t OutputSize() const override;
  size_t CopyOutput(uint8_t* buffer, size_t len) override;

 private:
  bool InitForFormat();
  bool ProcessData(const uint8_t* data, size_t len, size_t chunk_pos,
                   size_t chunk_size);
  void WriteHeader(size_t chunk_size);
  void CopyBlock(const uint8_t* data, size_t len);
  bool ProcessBlock();

  // wav reader callbacks
  static bool ParseFormatCb(void* opaque, const uint8_t* data, size_t len,
                            size_t chunk_pos, size_t chunk_size);
  static bool ProcessDataCb(void* opaque, const uint8_t* data, size_t len,
                            size_t chunk_pos, size_t chunk_size);

  int bitrate_;
  std::string opus_data_;
  StreamingWavReader wav_reader_;
  FormatChunk format_;
  size_t output_pos_;

  std::unique_ptr<OpusMSEncoder, void (*)(OpusMSEncoder*)> enc_;
  std::vector<int16_t> block_in_;
  int pre_skip_;
  size_t framesize_;
  size_t input_samples_;
  const uint32_t serial_number_ = 0;
  uint32_t sequence_number_;
  std::vector<uint8_t> page_data_;
  uint64_t total_encoded_;
  size_t page_pos_;
  int num_segments_;
  std::vector<uint16_t> packet_lengths_;
};

class StreamingOpusDecoder : public StreamingInterface {
 public:
  StreamingOpusDecoder() { Reset(); }

  void Reset() override;
  bool ProcessInput(const uint8_t* data, size_t len) override;
  bool Flush() override;
  size_t OutputSize() const override;
  size_t CopyOutput(uint8_t* buffer, size_t len) override;

 private:
  std::string opus_data_;
  std::string wav_data_;
  size_t output_pos_;
};

class StreamingOpusCodec : public StreamingAudioCodec {
 public:
  StreamingInterface* encoder() override {
    if (!encoder_) {
      encoder_ = std::make_unique<StreamingOpusEncoder>(bitrate_);
    }
    return encoder_.get();
  }
  StreamingInterface* decoder() override {
    if (!decoder_) {
      decoder_ = std::make_unique<StreamingOpusDecoder>();
    }
    return decoder_.get();
  }
  std::string Name() const override { return "opus"; }

 protected:
  bool ParseParam(const std::string& param) override {
    if (param.compare(0, 2, "br") == 0) {
      bitrate_ = std::stod(param.substr(2));
    } else {
      return false;
    }
    return true;
  }

  std::string ParamsToString() const override {
    return absl::Substitute("br$0", bitrate_);
  }

 private:
  int bitrate_ = 128;
  std::unique_ptr<StreamingOpusEncoder> encoder_;
  std::unique_ptr<StreamingOpusDecoder> decoder_;
};

}  // namespace ringli

#endif  // ANALYSIS_OPUS_CODEC_H_
