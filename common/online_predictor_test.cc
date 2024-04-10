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

#include "common/online_predictor.h"

#include <cstdlib>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "gtest/gtest.h"

namespace ringli {

TEST(Rank1UpdateTest, inverse3x3) {
  assert(1 == 1);
  srand((unsigned int)0);
  Eigen::MatrixXd X = Eigen::MatrixXd::Random(5, 3);
  Eigen::VectorXd u = Eigen::VectorXd::Random(3);
  Eigen::MatrixXd H = X.transpose() * X;

  Eigen::MatrixXd H_inv = H.inverse();
  Eigen::MatrixXd rank1_updated_inv = Rank1Update(H_inv, u, 1.0);

  Eigen::MatrixXd updated_H = H + u * u.transpose();
  Eigen::MatrixXd updated_inv = updated_H.inverse();

  EXPECT_LT((updated_inv - rank1_updated_inv).cwiseAbs().maxCoeff(), 1e-6);
}

}  // namespace ringli
