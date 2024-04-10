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

#ifndef DECODE_CONTEXT_MAP_DECODE_H_
#define DECODE_CONTEXT_MAP_DECODE_H_

#include <stdint.h>

#include "decode/bit_reader.h"

namespace ringli {

// Reads the context map from the bit stream. The context map is an array of
// context_map_size histogram ids. The number of different histogram ids is
// returned in *num_histograms.
bool DecodeContextMap(int context_map_size, uint8_t* context_map,
                      int* num_histograms, RingliBitReader* br);

}  // namespace ringli

#endif  // DECODE_CONTEXT_MAP_DECODE_H_
