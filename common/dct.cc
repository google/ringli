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

#include "common/dct.h"

namespace ringli {

constexpr float WcMultipliers<4>::kMultipliers[];
constexpr float WcMultipliers<8>::kMultipliers[];
constexpr float WcMultipliers<16>::kMultipliers[];
constexpr float WcMultipliers<32>::kMultipliers[];
constexpr float WcMultipliers<64>::kMultipliers[];

}  // namespace ringli
