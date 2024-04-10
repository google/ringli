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

#include "common/data_defs/data_matrix.h"

#include "common/data_defs/data_vector.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace ringli {
namespace {

using ::testing::ElementsAre;

TEST(DataMatrixTest, CanCreateFromData) {
  DataMatrix<double, 2> m((double[2][2]){{1, 2}, {3, 4}});

  EXPECT_THAT(m.ToStdVector(),
              ElementsAre(ElementsAre(1, 2), ElementsAre(3, 4)));
}

TEST(DataMatrixTest, ConstructsFromOtherMatrix) {
  DataMatrix<double, 2> m1((double[2][2]){{1, 2}, {3, 4}});

  DataMatrix<double, 2> m2(m1);

  EXPECT_THAT(m2.Row(0).ToStdVector(), ElementsAre(1, 2));
  EXPECT_THAT(m2.Row(1).ToStdVector(), ElementsAre(3, 4));
}

TEST(DataMatrixTest, AssignsFromOtherMatrix) {
  DataMatrix<double, 2> m2;
  {
    DataMatrix<double, 2> m1((double[2][2]){{1, 2}, {3, 4}});
    m2 = m1;
    m1[0][0] = 0;
    // m1 goes out of scope.
  }

  EXPECT_THAT(m2.Row(0).ToStdVector(), ElementsAre(1, 2));
  EXPECT_THAT(m2.Row(1).ToStdVector(), ElementsAre(3, 4));
}

TEST(DataMatrixTest, MultiplyIdentityToVector) {
  DataMatrix<double, 2> m((double[2][2]){{1, 0}, {0, 1}});

  DataVector<double, 2> v((double[]){5, 6});

  EXPECT_THAT((m * v).ToStdVector(), ElementsAre(5, 6));
}

TEST(DataMatrixTest, MultiplyToMatrix) {
  DataMatrix<double, 2> m1((double[2][2]){{0, 1}, {1, 0}});
  DataMatrix<double, 2> m2((double[2][2]){{5, 6}, {7, 8}});

  DataMatrix<double, 2> m3(m1 * m2);

  EXPECT_THAT(m3.Row(0).ToStdVector(), ElementsAre(7, 8));
  EXPECT_THAT(m3.Row(1).ToStdVector(), ElementsAre(5, 6));
}

TEST(DataMatrixTest, Transposes) {
  DataMatrix<double, 2> m1((double[2][2]){{1, 2}, {3, 4}});

  DataMatrix<double, 2> m2 = m1.Transposed();

  EXPECT_THAT(m2.Row(0).ToStdVector(), ElementsAre(1, 3));
  EXPECT_THAT(m2.Row(1).ToStdVector(), ElementsAre(2, 4));
}

}  // namespace
}  // namespace ringli
