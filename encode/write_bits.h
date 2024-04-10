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

/* Write bits into a byte array. */

#ifndef ENCODE_WRITE_BITS_H_
#define ENCODE_WRITE_BITS_H_

#include <stddef.h>
#include <stdint.h>

#include <iomanip>
#include <ios>

#include "absl/log/check.h"
#include "absl/log/log.h"

namespace ringli {

/* This function writes bits into bytes in increasing addresses, and within
   a byte least-significant-bit first.

   The function can write up to 56 bits in one go with WriteBits
   Example: let's assume that 3 bits (Rs below) have been written already:

   BYTE-0     BYTE+1       BYTE+2

   0000 0RRR    0000 0000    0000 0000

   Now, we could write 5 or less bits in MSB by just shifting by 3
   and OR'ing to BYTE-0.

   For n bits, we take the last 5 bits, OR that with high bits in BYTE-0,
   and locate the rest in BYTE+1, BYTE+2, etc. */
inline void WriteBits(size_t n_bits, uint64_t bits, size_t* pos,
                      uint8_t* array) {
  VLOG(1) << "WriteBits  " << std::setw(2) << n_bits << "  " << std::hex
          << std::setw(16) << bits << "  " << std::dec << std::setw(10) << *pos;
  DCHECK_EQ(bits >> n_bits, 0);
  DCHECK_LE(n_bits, 56);
  // TODO(szabadka) Add back the little-endian optimization
  /* implicit & 0xFF is assumed for uint8_t arithmetics */
  uint8_t* array_pos = &array[*pos >> 3];
  const size_t bits_reserved_in_first_byte = (*pos & 7);
  bits <<= bits_reserved_in_first_byte;
  *array_pos++ |= static_cast<uint8_t>(bits);
  for (size_t bits_left_to_write = n_bits + bits_reserved_in_first_byte;
       bits_left_to_write >= 9; bits_left_to_write -= 8) {
    bits >>= 8;
    *array_pos++ = static_cast<uint8_t>(bits);
  }
  *array_pos = 0;
  *pos += n_bits;
}

inline void WriteBitsPrepareStorage(size_t pos, uint8_t* array) {
  VLOG(1) << "WriteBitsPrepareStorage            " << std::setw(10) << pos;
  DCHECK_EQ(pos & 7, 0);
  array[pos >> 3] = 0;
}

}  // namespace ringli

#endif  // ENCODE_WRITE_BITS_H_
