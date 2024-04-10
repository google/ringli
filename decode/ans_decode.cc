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

#include "decode/ans_decode.h"

#include <vector>

#include "common/ans_params.h"
#include "common/entropy_coding.h"
#include "decode/bit_reader.h"
#include "decode/histogram_decode.h"

namespace ringli {

namespace {

bool ANSBuildMapTable(const int counts[MAX_SYMBOLS],
                      ANSSymbolInfo map[ANS_TAB_SIZE]) {
  int i;
  int pos = 0;
  for (i = 0; i < MAX_SYMBOLS; ++i) {
    int j;
    for (j = 0; j < counts[i]; ++j, ++pos) {
      map[pos].symbol_ = i;
      map[pos].freq_ = counts[i];
      map[pos].offset_ = j;
    }
  }
  return (pos == ANS_TAB_SIZE);
}

}  // namespace

bool ANSDecodingData::ReadFromBitStream(RingliBitReader* br) {
  int counts[MAX_SYMBOLS];
  return (ReadHistogram(ANS_LOG_TAB_SIZE, counts, br) &&
          ANSBuildMapTable(counts, map_));
}

}  // namespace ringli
