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

#ifndef DECODE_RINGLI_INPUT_H_
#define DECODE_RINGLI_INPUT_H_

#include <stddef.h>
#include <stdint.h>

namespace ringli {

static const int kBitMask[] = {
    0, 1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047, 4095, 8191, 16383, 32767,
};

class RingliInput {
 public:
  RingliInput(const uint8_t* data, size_t len)
      : data_(data), len_(len), pos_(0), val_(0), bit_pos_(0), error_(0) {}

  void InitBitReader() { val_ = GetNextWord(); }

  uint16_t GetNextWord() {
    uint16_t val = 0;
    if (pos_ + 1 < len_) {
      val = data_[pos_] + (data_[pos_ + 1] << 8);
    } else {
      error_ = 1;
    }
    pos_ += 2;
    return val;
  }

  int ReadBits(int nbits) {
    if (bit_pos_ + nbits > 16) {
      uint32_t new_bits = GetNextWord();
      val_ |= new_bits << 16;
    }
    int retval = (val_ >> bit_pos_) & kBitMask[nbits];
    bit_pos_ += nbits;
    if (bit_pos_ > 16) {
      bit_pos_ -= 16;
      val_ >>= 16;
    }
    return retval;
  }

  bool ok() const { return !error_; }

 private:
  const uint8_t* data_;
  size_t len_;
  size_t pos_;
  uint32_t val_;
  int bit_pos_;
  int error_;
};

}  // namespace ringli

#endif  // DECODE_RINGLI_INPUT_H_
