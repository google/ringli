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

#ifndef DECODE_RINGLI_DECODER_H_
#define DECODE_RINGLI_DECODER_H_

#include <memory>
#include <string>
#include <vector>

#include "common/adaptive_quant.h"
#include "common/data_defs/constants.h"
#include "common/dct.h"
#include "common/predictor.h"
#include "common/ringli_header.h"
#include "common/streaming.h"
#include "decode/entropy_decode.h"
#include "decode/noise_filtering.h"

namespace ringli {

class StreamingRingliDecoder : public StreamingInterface {
 public:
  StreamingRingliDecoder() { StreamingRingliDecoder::Reset(); }

  void Reset() override;
  bool ProcessInput(const uint8_t* data, size_t len) override;
  bool Flush() override;
  size_t OutputSize() const override;
  size_t CopyOutput(uint8_t* buffer, size_t len) override;

 private:
  bool ProcessHeader(const uint8_t* data, size_t len);
  bool ProcessBlock(const RingliBlock& block);
  bool ProcessSamples(const int* samples);
  void WriteBlock(const AudioBlock& block);

  static bool ProcessSamplesCb(void* opaque, const int* samples) {
    return reinterpret_cast<StreamingRingliDecoder*>(opaque)->ProcessSamples(
        samples);
  }

  std::string ringli_data_;
  std::string wav_data_;
  RingliHeader ringli_header_;
  std::unique_ptr<DCT<kDctLength>> dct_;
  std::unique_ptr<RingliBlock> prev_;
  std::unique_ptr<RingliBlock> current_;
  std::unique_ptr<RingliBlock> next_;
  std::unique_ptr<EntropyDecoder> entropy_decoder_;
  std::vector<std::unique_ptr<Predictor>> predictors_;
  std::vector<SymNoiseFilter> noise_filters_;
  std::vector<AdaptiveQuantizer> adaptive_quantizers_;
  size_t idx_;
  size_t num_blocks_;
  size_t input_pos_;
  size_t output_pos_;
  size_t samples_written_;
  size_t remaining_samples_;
};

bool RingliDecompress(const std::string& ringli_data, std::string* wav_data);

}  // namespace ringli

#endif  // DECODE_RINGLI_DECODER_H_
