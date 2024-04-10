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

#ifndef COMMON_COVARIANCE_LATTICE_H_
#define COMMON_COVARIANCE_LATTICE_H_

#include <cstddef>
#include <vector>

#include "Eigen/Core"
#include "absl/log/check.h"
namespace ringli {

// Implementation of the "covariance lattice method" for computing optimal
// linear predictor parameters, as described in the following paper:
// J. Makhoul, "New lattice methods for linear prediction", in IEEE Int. Conf.
// Acoust., Speech Signal Process., pp 462-465
template <typename T>
class CovarianceLattice {
 public:
  CovarianceLattice(const T* data, size_t len, size_t max_order,
                    double regulariser)
      : data_(data),
        len_(len),
        max_order_(max_order),
        last_order_(0),
        covariance_(Eigen::MatrixXd::Zero(max_order_ + 1, max_order_ + 1)),
        reflection_coeffs_(max_order_ + 1),
        a_(max_order_ + 1) {
    for (int n = max_order_; n < len_; ++n) {
      covariance_(0, 0) += data_[n] * data_[n];
    }
    covariance_(0, 0) += regulariser * (len_ - max_order_);
  }

  void FitCoeffsCommon(int order) {
    CHECK_GT(order, last_order_);
    for (int i = last_order_ + 1; i <= order; ++i) {
      for (int n = max_order_; n < len_; ++n) {
        covariance_(0, i) += data_[n] * data_[n - i];
      }
    }
    for (int i = 0; i < order; ++i) {
      const int prev_i = max_order_ - 1 - i;
      const int next_i = len_ - 1 - i;
      for (int j = last_order_; j < order; ++j) {
        const int prev_j = max_order_ - 1 - j;
        const int next_j = len_ - 1 - j;
        covariance_(i + 1, j + 1) = covariance_(i, j) +
                                    data_[prev_i] * data_[prev_j] -
                                    data_[next_i] * data_[next_j];
      }
    }
    for (int i = last_order_ + 1; i <= order; ++i) {
      for (int j = 0; j < i; ++j) {
        covariance_(i, j) = covariance_(j, i);
      }
    }
    for (int m = last_order_; m < order; ++m) {
      double num = covariance_(0, m + 1);
      for (int k = 1; k <= m; ++k) {
        num += a_[k] * (covariance_(0, m + 1 - k) + covariance_(k, m + 1));
        num += a_[k] * a_[k] * covariance_(k, m + 1 - k);
      }
      for (int k = 1; k < m; ++k) {
        for (int i = k + 1; i <= m; ++i) {
          num += a_[k] * a_[i] *
                 (covariance_(k, m + 1 - i) + covariance_(i, m + 1 - k));
        }
      }
      double denom = covariance_(0, 0) + covariance_(m + 1, m + 1);
      for (int k = 1; k <= m; ++k) {
        denom +=
            2 * a_[k] * (covariance_(0, k) + covariance_(m + 1, m + 1 - k));
        denom += a_[k] * a_[k] *
                 (covariance_(k, k) + covariance_(m + 1 - k, m + 1 - k));
      }
      for (int k = 1; k < m; ++k) {
        for (int i = k + 1; i <= m; ++i) {
          denom += 2 * a_[k] * a_[i] *
                   (covariance_(k, i) + covariance_(m + 1 - k, m + 1 - i));
        }
      }
      CHECK_GE(denom, 0);
      reflection_coeffs_[m + 1] = denom > 1e-3 ? -2 * num / denom : 0.0;
      std::vector<double> next_a(max_order_ + 1);
      next_a[m + 1] = reflection_coeffs_[m + 1];
      for (int j = 1; j <= m; ++j) {
        next_a[j] = a_[j] + reflection_coeffs_[m + 1] * a_[m + 1 - j];
      }
      a_ = next_a;
    }
    last_order_ = order;
  }

  template <typename TCoef>
  void FitPredictorCoeffs(TCoef* pcoefs, int order) {
    FitCoeffsCommon(order);
    for (int p = 0; p < order; ++p) {
      pcoefs[p] = -a_[p + 1];
    }
  }

  template <typename TCoef>
  void FitReflectionCoeffs(TCoef* k, int order) {
    FitCoeffsCommon(order);
    for (int p = 0; p <= order; ++p) {
      k[p] = reflection_coeffs_[p];
    }
  }

 private:
  const T* const data_;
  const int len_;
  const int max_order_;
  int last_order_;
  Eigen::MatrixXd covariance_;
  std::vector<double> reflection_coeffs_;
  std::vector<double> a_;
};

}  // namespace ringli

#endif  // COMMON_COVARIANCE_LATTICE_H_
