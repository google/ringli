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

#include "Eigen/Core"
#include "Eigen/Dense"
#include "common/data_defs/constants.h"

namespace ringli {

float OnlinePredictor::Predict() {
  float prediction = 0;
  if (position_ == 0) {
    // do nothing
  } else if (position_ == 1) {
    // predict first sample
    prediction = data_buffer_[0];
  } else if (position_ < 2 * kOnlinePredictorOrder) {
    // use precomputed optimal order 2 predictor
    prediction =
        kOnlineO2PredictorDefaultCoeffs[0] * data_buffer_[position_ - 1] -
        kOnlineO2PredictorDefaultCoeffs[1] * data_buffer_[position_ - 2];
  } else {
    for (int i = 1; i <= kOnlinePredictorOrder; ++i) {
      prediction += data_buffer_[(position_ - i) % kOnlinePredictorBufferSize] *
                    coeffs_[i];
    }
  }
  return prediction;
}

void OnlinePredictor::AddNewSample(float sample) {
  if (position_ >= kOnlinePredictorOrder) {
    // add new feature vector to covariance and gradient
    Eigen::VectorXd x_in = Eigen::VectorXd::Zero(kOnlinePredictorOrder);
    for (int j = 0; j < kOnlinePredictorOrder; ++j) {
      x_in(j) = data_buffer_[circular_index(-j - 1)];
    }
    covariance_ = Rank1Update(covariance_, x_in, 1.0);
    gradient_ = gradient_ + x_in * sample;

    if (position_ >= kOnlinePredictorOrder + kOnlinePredictorBufferSize) {
      // remove old feature vector from covariance and gradient
      Eigen::VectorXd x_out = Eigen::VectorXd::Zero(kOnlinePredictorOrder);
      for (int j = 0; j < kOnlinePredictorOrder; ++j) {
        x_out(j) = data_buffer_[circular_index(kOnlinePredictorOrder - j - 1)];
      }
      covariance_ = Rank1Update(covariance_, x_out, -1.0);
      gradient_ = gradient_ -
                  x_out * data_buffer_[circular_index(kOnlinePredictorOrder)];
    }

    Eigen::VectorXd solution = covariance_ * gradient_;

    for (int i = 1; i <= kOnlinePredictorOrder; ++i) {
      coeffs_[i] = solution[i - 1];
    }
  }

  data_buffer_[position_ % kOnlinePredictorBufferSize] = sample;
  position_++;
}

}  // namespace ringli
