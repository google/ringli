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

#ifndef COMMON_ENTROPY_CODING_H_
#define COMMON_ENTROPY_CODING_H_

#include <stdint.h>

namespace ringli {

// TODO(szabadka) Send these parameters in the bit-stream so we can choose
// different ones based on coding method and stream statistics.
constexpr int NUM_DIRECT_CODES = 4;
constexpr int MAX_SYMBOLS = 256;

constexpr int kPredNumDirectAbsval = 4;

struct EntropyCodingParams {
  uint8_t arithmetic_only = 0;
} __attribute__((packed));

}  // namespace ringli

#endif  // COMMON_ENTROPY_CODING_H_
