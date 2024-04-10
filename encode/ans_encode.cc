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

#include "encode/ans_encode.h"

#include <stddef.h>
#include <stdint.h>

#include <vector>

#include "common/ans_params.h"
#include "common/entropy_coding.h"
#include "encode/histogram_encode.h"

namespace ringli {

namespace {

void ANSBuildInfoTable(const int* counts, int alphabet_size,
                       ANSEncSymbolInfo info[MAX_SYMBOLS]) {
  int total = 0;
  for (int s = 0; s < alphabet_size; ++s) {
    const uint32_t freq = counts[s];
    info[s].freq_ = counts[s];
    info[s].start_ = total;
    total += freq;
#ifdef USE_MULT_BY_RECIPROCAL
    if (freq != 0) {
      info[s].ifreq_ =
          ((1ull << RECIPROCAL_PRECISION) + info[s].freq_ - 1) / info[s].freq_;
    } else {
      info[s].ifreq_ = 1;  // shouldn't matter (symbol shoudln't occur), but...
    }
#endif
  }
}

}  // namespace

void BuildAndStoreANSEncodingData(const int* histogram, ANSTable* table,
                                  size_t* storage_ix, uint8_t* storage) {
  int num_symbols;
  int symbols[kMaxNumSymbolsForSmallCode] = {0};
  std::vector<int> counts(histogram, histogram + MAX_SYMBOLS);
  int omit_pos = 0;
  NormalizeCounts(&counts[0], &omit_pos, MAX_SYMBOLS, ANS_LOG_TAB_SIZE,
                  &num_symbols, symbols);
  ANSBuildInfoTable(&counts[0], MAX_SYMBOLS, table->info_);
  EncodeCounts(&counts[0], omit_pos, num_symbols, symbols, storage_ix, storage);
}

}  // namespace ringli
