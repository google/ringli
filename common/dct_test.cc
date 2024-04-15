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

#include <cmath>
#include <functional>

#include "common/data_defs/constants.h"
#include "common/data_defs/data_vector.h"
#include "gtest/gtest.h"

namespace ringli {

template <int SIZE>
class DCT {
 public:
  DataVector<float, SIZE> ApplyDirectDCT(
      const DataVector<float, SIZE>& input) const {
    return direct_matrix_ * input;
  }

  DataVector<float, SIZE> ApplyInverseDCT(
      const DataVector<float, SIZE>& input) const {
    return inverse_matrix_ * input;
  }

 private:
  static DataMatrix<float, SIZE> CreateDirectMatrix() {
    DataMatrix<float, SIZE> directDataMatrix;
    for (int i = 0; i < SIZE; i++) {
      float coefI = i == 0 ? sqrt(1.0 / SIZE) : sqrt(2.0 / SIZE);
      for (int j = 0; j < SIZE; j++) {
        directDataMatrix[i][j] = coefI * cos(M_PI / SIZE * (j + 0.5) * i);
      }
    }

    return directDataMatrix;
  }

  const DataMatrix<float, SIZE> direct_matrix_ = CreateDirectMatrix();
  const DataMatrix<float, SIZE> inverse_matrix_ = direct_matrix_.Transposed();
};

TEST(DCTTest, VerifyInverse) {
  DCT<kDctLength> dct;

  for (int i = 0; i < kDctLength; ++i) {
    DataVector<float, kDctLength> e;
    e[i] = 1.0;
    DataVector<float, kDctLength> e1 = dct.ApplyDirectDCT(e);
    DataVector<float, kDctLength> e2 = dct.ApplyInverseDCT(e1);
    DataVector<float, kDctLength> error = e - e2;
    EXPECT_LT(error.AbsMax(), 1e-7);
  }
}

TEST(DCTTest, VerifyOrthogonal) {
  DCT<kDctLength> dct;

  for (int i = 0; i < kDctLength; ++i) {
    DataVector<float, kDctLength> ei;
    ei[i] = 1.0;
    DataVector<float, kDctLength> dct_ei = dct.ApplyDirectDCT(ei);
    for (int j = 0; j < kDctLength; ++j) {
      DataVector<float, kDctLength> ej;
      ej[j] = 1.0;
      DataVector<float, kDctLength> dct_ej = dct.ApplyDirectDCT(ej);
      float v = dct_ei.dot(dct_ej);
      float expected = (i == j ? 1.0 : 0.0);
      EXPECT_NEAR(v, expected, 5e-7);
    }
  }
}

TEST(DCTTest, TestForwardDCT) {
  DCT<kDctLength> dct;

  for (int i = 0; i < kDctLength; ++i) {
    DataVector<float, kDctLength> e;
    e[i] = 1.0;
    DataVector<float, kDctLength> e1 = dct.ApplyDirectDCT(e);
    DataVector<float, kDctLength> e2 = ForwardDCT(e);
    DataVector<float, kDctLength> error = e1 - e2;
    EXPECT_LT(error.AbsMax(), 3e-7);
  }
  {
    DataVector<float, kDctLength> e;
    for (int i = 0; i < kDctLength; ++i) {
      e[i] = 1.0;
    }
    DataVector<float, kDctLength> e1 = dct.ApplyDirectDCT(e);
    DataVector<float, kDctLength> e2 = ForwardDCT(e);
    DataVector<float, kDctLength> error = e1 - e2;
    EXPECT_LT(error.AbsMax(), 5e-8);
  }
}

TEST(DCTTest, TestInverseDCT) {
  DCT<kDctLength> dct;

  for (int i = 0; i < kDctLength; ++i) {
    DataVector<float, kDctLength> e;
    e[i] = 1.0;
    DataVector<float, kDctLength> e1 = dct.ApplyInverseDCT(e);
    DataVector<float, kDctLength> e2 = InverseDCT(e);
    DataVector<float, kDctLength> error = e1 - e2;
    EXPECT_LT(error.AbsMax(), 5e-7);
  }
  {
    DataVector<float, kDctLength> e;
    for (int i = 0; i < kDctLength; ++i) {
      e[i] = 1.0;
    }
    DataVector<float, kDctLength> e1 = dct.ApplyInverseDCT(e);
    DataVector<float, kDctLength> e2 = InverseDCT(e);
    DataVector<float, kDctLength> error = e1 - e2;
    EXPECT_LT(error.AbsMax(), 5e-6);
  }
}

}  // namespace ringli
