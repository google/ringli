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

// Library to decode the Huffman code lengths from the bit-stream and build a
// decoding table from them.

#ifndef DECODE_HUFFMAN_DECODE_H_
#define DECODE_HUFFMAN_DECODE_H_

#include <vector>

#include "decode/bit_reader.h"
#include "decode/huffman_table.h"

namespace ringli {

static const int kHuffmanTableMask = 0xff;
static const int kHuffmanTableBits = 8;
static const int kMaxHuffmanTableSize = 2048;

struct HuffmanDecodingData {
  HuffmanDecodingData() : table_(kMaxHuffmanTableSize) {}

  // Decodes the Huffman code lengths from the bit-stream and fills in the
  // pre-allocated table with the corresponding 2-level Huffman decoding
  // table. Returns false if the Huffman code lengths can not de decoded.
  bool ReadFromBitStream(int alphabet_size, RingliBitReader* br);

  std::vector<HuffmanCode> table_;
};

struct HuffmanDecoder {
  // Decodes the next Huffman coded symbol from the bit-stream.
  int ReadSymbol(const HuffmanDecodingData& code, RingliBitReader* br) {
    int nbits;
    RingliBitReaderFillWindow(br, 16);
    const HuffmanCode* table = &code.table_[0];
    table += (br->val >> br->bit_pos) & kHuffmanTableMask;
    nbits = table->bits - kHuffmanTableBits;
    if (nbits > 0) {
      br->bit_pos += kHuffmanTableBits;
      table += table->value;
      table += (br->val >> br->bit_pos) & ((1 << nbits) - 1);
    }
    br->bit_pos += table->bits;
    return table->value;
  }
};

}  // namespace ringli

#endif  // DECODE_HUFFMAN_DECODE_H_
