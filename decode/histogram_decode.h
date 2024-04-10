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

// Library to decode a histogram from the bit-stream.

#ifndef DECODE_HISTOGRAM_DECODE_H_
#define DECODE_HISTOGRAM_DECODE_H_

#include "decode/bit_reader.h"

namespace ringli {

// Decodes a histogram from the bit-stream where the sum of all population
// counts is 1 << precision_bits.
// Fills in counts[0 .. MAX_SYMBOLS) with the decoded population count values.
// Returns false on decoding error.
bool ReadHistogram(int precision_bits, int* counts, RingliBitReader* br);

}  // namespace ringli

#endif  // DECODE_HISTOGRAM_DECODE_H_
