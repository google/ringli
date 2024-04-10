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

#include "decode/context_map_decode.h"

#include <stddef.h>
#include <stdint.h>

#include <cstring> /* for memset */
#include <vector>

#include "absl/log/log.h"
#include "decode/bit_reader.h"
#include "decode/huffman_decode.h"
#include "decode/huffman_table.h"

namespace ringli {

namespace {

// Decodes a number in the range [0..255], by reading 1 - 11 bits.
inline int DecodeVarLenUint8(RingliBitReader* br) {
  if (RingliBitReaderReadBits(br, 1)) {
    int nbits = static_cast<int>(RingliBitReaderReadBits(br, 3));
    if (nbits == 0) {
      return 1;
    } else {
      return static_cast<int>(RingliBitReaderReadBits(br, nbits)) +
             (1 << nbits);
    }
  }
  return 0;
}

void MoveToFront(uint8_t* v, uint8_t index) {
  uint8_t value = v[index];
  uint8_t i = index;
  for (; i; --i) v[i] = v[i - 1];
  v[0] = value;
}

void InverseMoveToFrontTransform(uint8_t* v, int v_len) {
  uint8_t mtf[256];
  int i;
  for (i = 0; i < 256; ++i) {
    mtf[i] = (uint8_t)i;
  }
  for (i = 0; i < v_len; ++i) {
    uint8_t index = v[i];
    v[i] = mtf[index];
    if (index) MoveToFront(mtf, index);
  }
}

}  // namespace

bool DecodeContextMap(int context_map_size, uint8_t* context_map,
                      int* num_histograms, RingliBitReader* br) {
  *num_histograms = 1 + DecodeVarLenUint8(br);
  if (*num_histograms == 1) {
    memset(context_map, 0, (size_t)context_map_size);
    return true;
  }

  int max_run_length_prefix = 0;
  int use_rle_for_zeros = (int)RingliBitReaderReadBits(br, 1);
  if (use_rle_for_zeros) {
    max_run_length_prefix = (int)RingliBitReaderReadBits(br, 4) + 1;
  }
  std::vector<HuffmanCode> table(kMaxHuffmanTableSize);
  HuffmanDecodingData entropy;
  if (!entropy.ReadFromBitStream(*num_histograms + max_run_length_prefix, br)) {
    return false;
  }
  HuffmanDecoder decoder;
  int i;
  for (i = 0; i < context_map_size;) {
    int code;
    if (!RingliBitReaderReadMoreInput(br)) {
      LOG(ERROR) << "[DecodeContextMap] Unexpected end of input.";
      return false;
    }
    code = decoder.ReadSymbol(entropy, br);
    if (code == 0) {
      context_map[i] = 0;
      ++i;
    } else if (code <= max_run_length_prefix) {
      int reps = 1 + (1 << code) + (int)RingliBitReaderReadBits(br, code);
      while (--reps) {
        if (i >= context_map_size) {
          return false;
        }
        context_map[i] = 0;
        ++i;
      }
    } else {
      context_map[i] = (uint8_t)(code - max_run_length_prefix);
      ++i;
    }
  }
  if (RingliBitReaderReadBits(br, 1)) {
    InverseMoveToFrontTransform(context_map, context_map_size);
  }
  return true;
}

}  // namespace ringli
