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

// Library to store a histogram to the bit-stream.

#ifndef ENCODE_HISTOGRAM_ENCODE_H_
#define ENCODE_HISTOGRAM_ENCODE_H_

#include <stddef.h>
#include <stdint.h>

#include "common/ans_params.h"

namespace ringli {

static const int kMaxNumSymbolsForSmallCode = 4;

// Normalizes the population counts in counts[0 .. length) so that the sum of
// all counts will be 1 << precision_bits.
// Sets *num_symbols to the number of symbols in the range [0 .. length) with
// non-zero population counts.
// Fills in symbols[0 .. kMaxNumSymbolsForSmallCode) with the first few symbols
// with non-zero population counts.
// Each count will all be rounded to multiples of
// 1 << GetPopulationCountPrecision(count), except possibly for one. The index
// of that count will be stored in *omit_pos.
void NormalizeCounts(int* counts, int* omit_pos, const int length,
                     const int precision_bits, int* num_symbols, int* symbols);

// Stores a histogram in counts[0 .. MAX_SYMBOLS) to the bit-stream where
// the sum of all population counts is ANS_TAB_SIZE and the number of symbols
// with non-zero counts is num_symbols.
// symbols[0 .. kMaxNumSymbolsForSmallCode) contains the first few symbols
// with non-zero population counts.
// Each count must be rounded to a multiple of
// 1 << GetPopulationCountPrecision(count), except possibly counts[omit_pos].
void EncodeCounts(const int* counts, const int omit_pos, const int num_symbols,
                  const int* symbols, size_t* storage_ix, uint8_t* storage);

// Returns an estimate of the number of bits required to encode the given
// histogram (header bits plus data bits).
double PopulationCost(const int* data, int total_count);

}  // namespace ringli

#endif  // ENCODE_HISTOGRAM_ENCODE_H_
