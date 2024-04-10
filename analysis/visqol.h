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

#ifndef ANALYSIS_VISQOL_H_
#define ANALYSIS_VISQOL_H_

#include <filesystem>

#include "hwy/aligned_allocator.h"

namespace ringli {

class ViSQOL {
 public:
  ViSQOL();
  ~ViSQOL();
  float MOS(const hwy::AlignedNDArray<float, 2>& reference,
            const hwy::AlignedNDArray<float, 2>& degraded,
            float sample_rate) const;

 private:
  std::filesystem::path model_path_;
};

}  // namespace ringli

#endif  // ANALYSIS_VISQOL_H_
