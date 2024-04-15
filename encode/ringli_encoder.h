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

#ifndef ENCODE_RINGLI_ENCODER_H_
#define ENCODE_RINGLI_ENCODER_H_

#include <stdbool.h>

#include <memory>
#include <string>
#include <vector>

#include "common/adaptive_quant.h"
#include "common/data_defs/constants.h"
#include "common/data_defs/data_vector.h"
#include "common/dct.h"
#include "common/predictor.h"
#include "common/ringli_header.h"
#include "common/segment_curve.h"
#include "common/streaming.h"
#include "common/wav_header.h"
#include "common/wav_reader.h"
#include "encode/entropy_encode.h"
#include "encode/noise_shaping.h"

namespace ringli {

enum QuantizationType {
  CONST = 0,
  ADAPTIVE_PER_PACKET = 1,
  ADAPTIVE_PER_BAND = 2,
};

struct RingliEncoderConfig {
  SegmentCurve quantization_curve =
      SegmentCurve(1.0);  // No quantization by default.
  QuantizationType quantization_type = QuantizationType::CONST;
  uint8_t pred_order = 8;
  uint8_t pred_order_min = 2;
  uint8_t pred_order_max = kMaxPredictorOrder;
  bool use_noise_shaping = false;
  RingliDecoderConfig dconfig;
};

class StreamingRingliEncoder : public StreamingInterface {
 public:
  explicit StreamingRingliEncoder(const RingliEncoderConfig& config);

  void Reset() override;
  bool ProcessInput(const uint8_t* data, size_t len) override;
  bool Flush() override;
  size_t OutputSize() const override;
  size_t CopyOutput(uint8_t* buffer, size_t len) override;

 private:
  void InitForFormat();
  bool ProcessData(const uint8_t* data, size_t len, size_t chunk_pos,
                   size_t chunk_size);
  void WriteHeader(size_t chunk_size);
  void CopyBlock(const uint8_t* data, size_t len, AudioBlock* block);
  bool ProcessBlock(const RingliBlock& ringli_block);

  // wav reader callbacks
  static bool ParseFormatCb(void* opaque, const uint8_t* data, size_t len,
                            size_t chunk_pos, size_t chunk_size);
  static bool ProcessDataCb(void* opaque, const uint8_t* data, size_t len,
                            size_t chunk_pos, size_t chunk_size);

  RingliEncoderConfig config_;
  std::string ringli_data_;
  StreamingWavReader wav_reader_;
  FormatChunk format_;
  size_t data_len_;
  size_t output_pos_ = 0;
  size_t idx_ = 0;

  std::unique_ptr<AudioBlock> prev_;
  std::unique_ptr<AudioBlock> current_;
  std::unique_ptr<AudioBlock> next_;
  std::unique_ptr<EntropyCoder> entropy_coder_;
  std::vector<RingliBlock> ringli_blocks_;
  std::vector<RingliBlockHeader> ringli_headers_;
  std::vector<std::unique_ptr<Predictor>> predictors_;
  std::vector<NoiseShaper> noise_shapers_;
  std::vector<AdaptiveQuantizer> adaptive_quantizers_;
};

void RingliCompress(const std::string& wav_data,
                    const RingliEncoderConfig& config,
                    std::string* ringli_data);

}  // namespace ringli
#endif  // ENCODE_RINGLI_ENCODER_H_
