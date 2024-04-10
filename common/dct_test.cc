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

TEST(DCTTest, VerifyInverse) {
  DCT<kDctLength> dct;

  for (int i = 0; i < kDctLength; ++i) {
    DataVector<double, kDctLength> e;
    e[i] = 1.0;
    DataVector<double, kDctLength> e1 = dct.ApplyDirectDCT(e);
    DataVector<double, kDctLength> e2 = dct.ApplyInverseDCT(e1);
    DataVector<double, kDctLength> error = e - e2;
    EXPECT_LT(error.AbsMax(), 1e-14);
  }
}

TEST(DCTTest, VerifyOrthogonal) {
  DCT<kDctLength> dct;

  for (int i = 0; i < kDctLength; ++i) {
    DataVector<double, kDctLength> ei;
    ei[i] = 1.0;
    DataVector<double, kDctLength> dct_ei = dct.ApplyDirectDCT(ei);
    for (int j = 0; j < kDctLength; ++j) {
      DataVector<double, kDctLength> ej;
      ej[j] = 1.0;
      DataVector<double, kDctLength> dct_ej = dct.ApplyDirectDCT(ej);
      float v = dct_ei.dot(dct_ej);
      float expected = (i == j ? 1.0 : 0.0);
      EXPECT_NEAR(v, expected, 1e-14);
    }
  }
}

}  // namespace ringli
