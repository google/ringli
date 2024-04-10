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

#ifndef DECODE_ARITH_DECODE_H_
#define DECODE_ARITH_DECODE_H_

#include <stdint.h>
#include <stdio.h>

#include "decode/ringli_input.h"

namespace ringli {

// A class used for entropy decoding a sequence of binary values.
// skal@ wrote the original version, szabadka@ ported it for brunsli.
class BinaryArithmeticDecoder {
 public:
  BinaryArithmeticDecoder() : low_(0), high_(0), value_(0) {}

  void Init(RingliInput* in) {
    value_ = in->GetNextWord();
    value_ = (value_ << 16) | in->GetNextWord();
    low_ = 0;
    high_ = ~0;
  }

  bool HasBit() const { return ((low_ ^ high_) >> 16) != 0; }

  void Fill(uint16_t next_word) {
    value_ = (value_ << 16) | next_word;
    low_ <<= 16;
    high_ <<= 16;
    high_ |= 0xffff;
  }

  // Returns the next bit decoded from the bit stream, based on the given
  // 8-bit precision probability, i.e. P(bit = 0) = prob / 256. This
  // probability must be the same as the one used by the encoder. Can be
  // called only when HasBit() returns true.
  int ReadBitNoFill(int prob) {
    const uint32_t diff = high_ - low_;
    const uint32_t split = low_ + (((uint64_t)diff * prob) >> 8);
    int bit;
    if (value_ > split) {
      low_ = split + 1;
      bit = 1;
    } else {
      high_ = split;
      bit = 0;
    }
    return bit;
  }

  // Same as above but fills in the state if necessary with new input values.
  int ReadBit(int prob, RingliInput* in) {
    while (!HasBit()) {
      Fill(in->GetNextWord());
    }
    return ReadBitNoFill(prob);
  }

  int ReadBits(int nbits, RingliInput* in) {
    int val = 0;
    for (int b = 0; b < nbits; ++b) {
      val |= ReadBit(128, in) << b;
    }
    return val;
  }

 private:
  uint32_t low_;
  uint32_t high_;
  uint32_t value_;
};

}  // namespace ringli

#endif  // DECODE_ARITH_DECODE_H_
