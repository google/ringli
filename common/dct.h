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

#ifndef COMMON_DCT_H_
#define COMMON_DCT_H_

#include <math.h>

#include "common/data_defs/data_matrix.h"
#include "common/data_defs/data_vector.h"

namespace ringli {

template <int SIZE>
class DCT {
 public:
  DataVector<double, SIZE> ApplyDirectDCT(
      const DataVector<double, SIZE>& input) const {
    return direct_matrix_ * input;
  }

  DataVector<double, SIZE> ApplyInverseDCT(
      const DataVector<double, SIZE>& input) const {
    return inverse_matrix_ * input;
  }

 private:
  static DataMatrix<double, SIZE> CreateDirectMatrix() {
    DataMatrix<double, SIZE> directDataMatrix;
    for (int i = 0; i < SIZE; i++) {
      double coefI = i == 0 ? sqrt(1.0 / SIZE) : sqrt(2.0 / SIZE);
      for (int j = 0; j < SIZE; j++) {
        directDataMatrix[i][j] = coefI * cos(M_PI / SIZE * (j + 0.5) * i);
      }
    }

    return directDataMatrix;
  }

  const DataMatrix<double, SIZE> direct_matrix_ = CreateDirectMatrix();
  const DataMatrix<double, SIZE> inverse_matrix_ = direct_matrix_.Transposed();
};

}  // namespace ringli

#endif  // COMMON_DCT_H_
