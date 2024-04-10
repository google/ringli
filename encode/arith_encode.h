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

#ifndef ENCODE_ARITH_ENCODE_H_
#define ENCODE_ARITH_ENCODE_H_

#include <stdint.h>

namespace ringli {

class BinaryArithmeticEncoder {
 public:
  BinaryArithmeticEncoder() { Reset(); }

  void Reset() {
    low_ = 0;
    high_ = ~0;
  }

  typedef void (*Output)(void* opaque, uint16_t val);

  void AddBit(uint8_t prob, int bit, void* opaque, Output output) {
    while (((low_ ^ high_) >> 16) == 0) {
      output(opaque, high_ >> 16);
      low_ <<= 16;
      high_ <<= 16;
      high_ |= 0xffff;
    }
    const uint32_t diff = high_ - low_;
    const uint32_t split = low_ + (((uint64_t)diff * prob) >> 8);
    if (bit) {
      low_ = split + 1;
    } else {
      high_ = split;
    }
  }

  void Flush(void* opaque, Output output) {
    output(opaque, high_ >> 16);
    output(opaque, high_ & 0xffff);
    Reset();
  }

 private:
  uint32_t low_;
  uint32_t high_;
};

}  // namespace ringli

#endif  // ENCODE_ARITH_ENCODE_H_
