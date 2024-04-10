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

#ifndef COMMON_ONLINE_PREDICTOR_H_
#define COMMON_ONLINE_PREDICTOR_H_

#include <stddef.h>
#include <stdint.h>

#include <cstring>

#include "Eigen/Core"
#include "absl/log/check.h"
#include "common/data_defs/constants.h"
#include "common/predictor.h"

namespace ringli {

// TODO(antoran): switch third argument to boolean
template <typename Derived1, typename Derived2>
Eigen::MatrixXd Rank1Update(const Eigen::MatrixBase<Derived1>& covariance,
                            const Eigen::MatrixBase<Derived2>& data,
                            double sign) {
  CHECK_EQ(std::abs(sign), 1.0);
  const Eigen::VectorXd cov_u = covariance * data;
  const double predictive_variance = (data.transpose() * cov_u);
  const double denominator = 1.0 + sign * predictive_variance;
  const Eigen::MatrixXd numerator = sign * cov_u * cov_u.transpose();
  return covariance - numerator / denominator;
}

class OnlinePredictor : public Predictor {
 public:
  explicit OnlinePredictor(float regulariser) : regulariser_(regulariser) {
    OnlinePredictor::Reset();
  }

  void Reset() override {
    position_ = 0;
    memset(coeffs_, 0, sizeof(coeffs_));
    memset(data_buffer_, 0, sizeof(data_buffer_));
    gradient_ = Eigen::VectorXd::Zero(kOnlinePredictorOrder);
    covariance_ =
        Eigen::MatrixXd::Identity(kOnlinePredictorOrder, kOnlinePredictorOrder);
    covariance_ = covariance_ * (1.0 / regulariser_);
  }

  float Predict() override;

  void AddNewSample(float sample) override;

  uint32_t circular_index(int i) const {
    return (position_ + i + kOnlinePredictorBufferSize) %
           kOnlinePredictorBufferSize;
  }

 private:
  const float regulariser_;
  uint32_t position_ = 0;
  float coeffs_[kOnlinePredictorOrder + 1] = {0};
  float data_buffer_[kOnlinePredictorBufferSize] = {0};
  Eigen::VectorXd gradient_;
  Eigen::MatrixXd covariance_;
};

}  // namespace ringli

#endif  // COMMON_ONLINE_PREDICTOR_H_
