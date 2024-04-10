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

#ifndef COMMON_LOGGING_H_
#define COMMON_LOGGING_H_

#include <stddef.h>
#include <stdio.h>

#include "absl/flags/declare.h"
#include "absl/flags/flag.h"

ABSL_DECLARE_FLAG(int, log_level);

namespace ringli {

static inline void PrintSize(const char* name, size_t len, double duration) {
  if (absl::GetFlag(FLAGS_log_level) >= 1) {
    printf("%-25s %10zu bytes  %8.3f kbps\n", name, len,
           len * 8.0 / 1024.0 / duration);
  }
}

template <typename T>
void PrintHistogram(const char* name, T* data, size_t len) {
  if (absl::GetFlag(FLAGS_log_level) >= 1) {
    printf("%s histogram: ", name);
    for (int i = 0; i < len; ++i) {
      if (data[i]) printf(" %d:%d", i, data[i]);
    }
    printf("\n");
  }
}

}  // namespace ringli

#endif  // COMMON_LOGGING_H_
