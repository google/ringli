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

#include "decode/histogram_decode.h"

#include <stdint.h>

#include "absl/log/log.h"
#include "common/ans_params.h"
#include "common/entropy_coding.h"
#include "common/log2floor.h"
#include "decode/bit_reader.h"
#include "decode/huffman_table.h"

namespace ringli {

bool ReadHistogram(int precision_bits, int* counts, RingliBitReader* br) {
  if (!RingliBitReaderReadMoreInput(br)) {
    LOG(ERROR) << "[ReadHistogram] Unexpected end of input.";
    return false;
  }
  const int max_bits = 1 + Log2FloorNonZero(MAX_SYMBOLS - 1);
  memset(counts, 0, MAX_SYMBOLS * sizeof(counts[0]));
  int simple_code = RingliBitReaderReadBits(br, 1);
  if (simple_code == 1) {
    int i;
    int symbols[2] = {0};
    const int num_symbols = RingliBitReaderReadBits(br, 1) + 1;
    for (i = 0; i < num_symbols; ++i) {
      symbols[i] = RingliBitReaderReadBits(br, max_bits) % MAX_SYMBOLS;
    }
    if (num_symbols == 1) {
      counts[symbols[0]] = 1 << precision_bits;
    } else {
      if (symbols[0] == symbols[1]) {  // corrupt data
        return false;
      }
      counts[symbols[0]] = RingliBitReaderReadBits(br, precision_bits);
      counts[symbols[1]] = (1 << precision_bits) - counts[symbols[0]];
    }
  } else {
    // Non-small tree means num_symbols >= 3, therefore length >= 3,
    // and we decode length - 3 from the bitstream.
    int length = 3 + RingliBitReaderReadBits(br, max_bits);
    int total_count = 0;
    static const HuffmanCode huff[64] = {
        {2, 6}, {3, 7}, {3, 4}, {4, 1}, {2, 6}, {3, 8}, {3, 5}, {4, 3},
        {2, 6}, {3, 7}, {3, 4}, {4, 2}, {2, 6}, {3, 8}, {3, 5}, {5, 0},
        {2, 6}, {3, 7}, {3, 4}, {4, 1}, {2, 6}, {3, 8}, {3, 5}, {4, 3},
        {2, 6}, {3, 7}, {3, 4}, {4, 2}, {2, 6}, {3, 8}, {3, 5}, {6, 9},
        {2, 6}, {3, 7}, {3, 4}, {4, 1}, {2, 6}, {3, 8}, {3, 5}, {4, 3},
        {2, 6}, {3, 7}, {3, 4}, {4, 2}, {2, 6}, {3, 8}, {3, 5}, {5, 0},
        {2, 6}, {3, 7}, {3, 4}, {4, 1}, {2, 6}, {3, 8}, {3, 5}, {4, 3},
        {2, 6}, {3, 7}, {3, 4}, {4, 2}, {2, 6}, {3, 8}, {3, 5}, {6, 10},
    };
    int logcounts[MAX_SYMBOLS];
    int omit_log = -1;
    int omit_pos = -1;
    for (int i = 0; i < length; ++i) {
      const HuffmanCode* p = huff;
      RingliBitReaderFillWindow(br, 6);
      p += (br->val >> br->bit_pos) & 63;
      br->bit_pos += p->bits;
      logcounts[i] = p->value;
      if (logcounts[i] > omit_log) {
        omit_log = logcounts[i];
        omit_pos = i;
      }
    }
    for (int i = 0; i < length; ++i) {
      int code = logcounts[i];
      if (i == omit_pos) {
        continue;
      } else if (code == 0) {
        continue;
      } else if (code == 1) {
        counts[i] = 1;
      } else {
        int bitcount = GetPopulationCountPrecision(code - 1);
        counts[i] = (1 << (code - 1)) + (RingliBitReaderReadBits(br, bitcount)
                                         << (code - 1 - bitcount));
      }
      total_count += counts[i];
    }
    if (omit_pos < 0 || (1 << precision_bits) <= total_count) {
      // The histogram we've read sums to more than total_count (including
      // at least 1 for the omitted value).
      return false;
    }
    counts[omit_pos] = (1 << precision_bits) - total_count;
  }
  return true;
}

}  // namespace ringli
