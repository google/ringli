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

// Library to decode the ANS population counts from the bit-stream and build a
// decoding table from them.

#ifndef DECODE_ANS_DECODE_H_
#define DECODE_ANS_DECODE_H_

#include <stdint.h>

#include "common/ans_params.h"
#include "decode/bit_reader.h"
#include "decode/ringli_input.h"

namespace ringli {

typedef struct {
  uint16_t offset_;
  uint16_t freq_;
  uint8_t symbol_;
} ANSSymbolInfo;

struct ANSDecodingData {
  ANSDecodingData() {}

  bool ReadFromBitStream(RingliBitReader* br);

  ANSSymbolInfo map_[ANS_TAB_SIZE];
};

class ANSDecoder {
 public:
  ANSDecoder() : state_(0) {}

  void Init(RingliInput* in) {
    state_ = in->GetNextWord();
    state_ = (state_ << 16) | in->GetNextWord();
  }

  int ReadSymbol(const ANSDecodingData& code, RingliInput* in) {
    const uint32_t res = state_ & (ANS_TAB_SIZE - 1);
    const ANSSymbolInfo& s = code.map_[res];
    state_ = s.freq_ * (state_ >> ANS_LOG_TAB_SIZE) + s.offset_;
    if (state_ < (1u << 16)) {
      state_ = (state_ << 16) | in->GetNextWord();
    }
    return s.symbol_;
  }
  bool CheckCRC() const { return state_ == (ANS_SIGNATURE << 16); }

 private:
  uint32_t state_;
};

}  // namespace ringli

#endif  // DECODE_ANS_DECODE_H_
