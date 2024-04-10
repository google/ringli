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

/* Copyright 2013 Google Inc. All Rights Reserved.
   Author: szabadka@google.com (Zoltan Szabadka)

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

/* Bit reading helpers */

#include "decode/bit_reader.h"

#include <cstring> /* for memset, memcpy */

namespace ringli {

void RingliBitReaderInit(RingliBitReader* br, const uint8_t* buffer,
                         size_t length) {
  br->src = buffer;
  br->available = length;

  memset(br->tail, 0, sizeof(br->tail));
  size_t tail_size =
      length > RINGLI_BIT_READER_SLACK ? RINGLI_BIT_READER_SLACK : length;
  memcpy(&br->tail[RINGLI_BIT_READER_SLACK - tail_size],
         &buffer[length - tail_size], tail_size);
  // Switch to tail right away, if input is too short.
  int tail_offset = RINGLI_BIT_READER_SLACK - br->available;
  if (tail_offset >= 0) {
    br->src = &br->tail[tail_offset];
  }

  br->val = 0;
  for (size_t i = 0; i < sizeof(br->val); ++i) {
    br->val |= ((uint32_t)(*br->src)) << (8 * i);
    ++br->src;
    --br->available;
  }
  br->bit_pos = 0;
}

}  // namespace ringli
