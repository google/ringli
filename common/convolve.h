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

#ifndef COMMON_CONVOLVE_H_
#define COMMON_CONVOLVE_H_

#include <algorithm>
#include <cmath>

#include "common/data_defs/data_vector.h"

namespace ringli {

template <int K>
DataVector<double, 2 * K + 1> GaussianKernel(double sigma) {
  DataVector<double, 2 * K + 1> res;
  double alpha = -0.5 / sigma / sigma;
  double sum = 0.0;
  for (int j = -K; j <= K; ++j) {
    res[K + j] = std::exp(alpha * j * j);
    sum += res[K + j];
  }
  for (int j = -K; j <= K; ++j) {
    res[K + j] /= sum;
  }
  return res;
}

template <int W, int SIZE>
DataVector<double, SIZE> Convolve(const DataVector<double, W>& kernel,
                                  const DataVector<double, SIZE>& data) {
  static_assert(SIZE >= W);
  static_assert(W % 2 == 1);
  constexpr int K = W / 2;
  DataVector<double, SIZE> res;
  for (int i = 0; i < SIZE; ++i) {
    double sum = 0;
    for (int j = -K; j <= K; ++j) {
      int idx = std::max(0, std::min(SIZE - 1, i + j));
      sum += data[idx] * kernel[K + j];
    }
    res[i] = sum;
  }
  return res;
}

}  // namespace ringli

#endif  // COMMON_CONVOLVE_H_
