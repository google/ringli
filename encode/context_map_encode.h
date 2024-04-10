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

#ifndef ENCODE_CONTEXT_MAP_ENCODE_H_
#define ENCODE_CONTEXT_MAP_ENCODE_H_

#include <stddef.h>
#include <stdint.h>

#include <vector>

namespace ringli {

// Encodes the given context map to the bit stream. The number of different
// histogram ids is given by num_clusters.
void EncodeContextMap(const std::vector<uint32_t>& context_map,
                      size_t num_clusters, size_t* storage_ix,
                      uint8_t* storage);

}  // namespace ringli

#endif  // ENCODE_CONTEXT_MAP_ENCODE_H_
