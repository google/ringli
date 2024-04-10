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

#ifndef COMMON_STREAMING_H_
#define COMMON_STREAMING_H_

#include <stddef.h>
#include <stdint.h>

namespace ringli {

class StreamingInterface {
 public:
  virtual ~StreamingInterface() = default;

  virtual void Reset() = 0;

  virtual bool ProcessInput(const uint8_t* data, size_t len) = 0;

  virtual bool Flush() = 0;

  virtual size_t OutputSize() const = 0;

  virtual size_t CopyOutput(uint8_t* buffer, size_t len) = 0;
};

}  // namespace ringli

#endif  // COMMON_STREAMING_H_
