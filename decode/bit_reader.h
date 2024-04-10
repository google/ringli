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

/* Copyright 2013 Google Inc. All Rights Reserved.
   Author: szabadka@google.com (Zoltan Szabadka)

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

/* Bit reading helpers */

#ifndef DECODE_BIT_READER_H_
#define DECODE_BIT_READER_H_

#include <stddef.h>
#include <stdint.h>

#include <iomanip>
#include <ios>

#include "absl/log/log.h"

namespace ringli {

// TODO(szabadka) Add back the little endian 64-bit optimizations.

#define RINGLI_BIT_READER_SLACK (32 + 8)

/* Masking with this expression turns to a single "Unsigned Bit Field Extract"
   UBFX instruction on ARM. */
static inline uint32_t RingliBitMask(int n) { return ~((0xffffffff) << n); }

typedef struct {
  uint32_t val;       /* pre-fetched bits */
  const uint8_t* src; /* next input */
  int available;      /* number of bytes available to read */
  uint32_t bit_pos;   /* current bit-reading position in val_ */

  /* Tail is used as data source when input is almost depleted. At most
   * RINGLI_BIT_READER_SLACK are user input, the remaining are zeros to make
   * a guarantee that RINGLI_BIT_READER_SLACK bytes would be safely read,
   * plus extra space for bulk (uint64_t) reads. */
  uint8_t tail[RINGLI_BIT_READER_SLACK * 2];
} RingliBitReader;

/* Initializes the bitreader fields. Returns 0 in case of failure. */
void RingliBitReaderInit(RingliBitReader* br, const uint8_t* buffer,
                         size_t length);

/* Prepare for further input reading. */
static inline int RingliBitReaderReadMoreInput(RingliBitReader* br) {
  if (br->available <= -8) {
    return 0;
  }
  int tail_offset = RINGLI_BIT_READER_SLACK - br->available;
  if (tail_offset >= 0) {
    br->src = &br->tail[tail_offset];
  }
  return 1;
}

/* Guarantees that there are at least n_bits in the buffer.
   n_bits should be in the range [1..24] */
static inline void RingliBitReaderFillWindow(RingliBitReader* br, int n_bits) {
  while (br->bit_pos >= 8) {
    br->val >>= 8;
    br->val |= ((uint32_t)*br->src) << 24;
    ++br->src;
    --br->available;
    br->bit_pos -= 8;
  }
}

/* Reads the specified number of bits from Read Buffer. */
static inline uint32_t RingliBitReaderReadBits(RingliBitReader* br,
                                               int n_bits) {
  uint32_t val;
  RingliBitReaderFillWindow(br, n_bits);
  val = (uint32_t)(br->val >> br->bit_pos) & RingliBitMask(n_bits);
  VLOG(1) << "[RingliReadBits]  " << br->bit_pos << " " << std::setw(2)
          << n_bits << "  val: " << std::setw(6) << std::hex << val;
  br->bit_pos += (uint32_t)n_bits;
  return val;
}

/* Returns the number of bytes available skipping bits. */
static inline size_t RingliBitReaderJumpToByteBoundary(RingliBitReader* br) {
  int nbits = br->bit_pos & 7;
  if (nbits > 0) {
    RingliBitReaderReadBits(br, 8 - nbits);
  }
  return br->available + sizeof(br->val) - (br->bit_pos >> 3);
}

}  // namespace ringli

#endif  // DECODE_BIT_READER_H_
