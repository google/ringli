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

#include "common/block_predictor.h"

#include <stdint.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <utility>
#include <vector>

#include "Eigen/Dense"
#include "absl/log/check.h"
#include "absl/types/span.h"
#include "common/data_defs/constants.h"
#include "common/ringli_header.h"

namespace ringli {
namespace {

// Returns true if the synthesis filter defined by the prediction coefficients
// is stable. This condition is sufficient for the existence of the line
// spectral frequencies.
static bool IsInversePredictionFilterStable(const float* pcoefs, int order) {
  // See eqs (50) and (51) in J.Makhoul, "Linear Prediction: A Tutorial
  // Review", Proceedings of the IEEE, vol. 63, pp. 561-580
  const int P = order;
  std::vector<double> reflected_coefs(P + 1);
  reflected_coefs[0] = 1.0;
  for (int p = 1; p <= P; ++p) {
    reflected_coefs[p] = -pcoefs[p - 1];
  }
  for (int i = P; i > 0; --i) {
    double ki = reflected_coefs[i];
    if (std::abs(ki) >= 1.0) {
      return false;
    }
    if (i > 1) {
      double scale = 1.0 / (1.0 - ki * ki);
      std::vector<double> tmp(P + 1);
      for (int j = 1; j < i; ++j) {
        tmp[j] = (reflected_coefs[j] - ki * reflected_coefs[i - j]) * scale;
      }
      reflected_coefs = tmp;
    }
  }
  return true;
}

double EvaluateChebyshevPolynomial(absl::Span<const double> c, double x) {
  double b1 = 0.0;
  double b2 = 0.0;
  for (int k = c.size() - 1; k > 0; --k) {
    double b = 2 * x * b1 - b2 + c[k];
    b2 = b1;
    b1 = b;
  }
  return x * b1 - b2 + c[0];
}

double FindRootOfChebyshevPolynomial(absl::Span<const double> c, double xa,
                                     double xb, double vala, double valb,
                                     absl::Span<const double> vals,
                                     uint16_t* quant_lsf) {
  CHECK_LT(xa, xb);
  CHECK_LE(vala * valb, 0);
  int idx_a = std::upper_bound(vals.begin(), vals.end(), xa) - vals.begin() - 1;
  int idx_b = std::lower_bound(vals.begin(), vals.end(), xb) - vals.begin();
  while (valb != 0 && idx_b > idx_a + 1) {
    const int idx = (idx_a + idx_b) / 2;
    const double x = vals[idx];
    const double val = EvaluateChebyshevPolynomial(c, x);
    if (val * vala < 0) {
      idx_b = idx;
      xb = x;
      valb = val;
    } else {
      idx_a = idx;
      xa = x;
      vala = val;
    }
  }
  *quant_lsf = vals.size() - 1 - idx_b;
  return (valb * xa - vala * xb) / (valb - vala);
}

class LSFQuantizer {
 public:
  LSFQuantizer() {
    for (int p = 0; p < kMaxPredictorOrder; ++p) {
      const float q = kLSFQuant[p];
      // Number of different possible quantized lsf values from 0 to
      // round(q)
      const int num_quant_vals = std::round(q) + 1;
      // Boundaries of the lsf quantization intervals
      std::vector<double> lsf_boundaries(num_quant_vals + 1);
      lsf_boundaries[0] = 0.0;
      for (int i = 0; i < num_quant_vals; ++i) {
        lsf_boundaries[i + 1] = std::min(1.0, (i + 0.5) / q);
      }
      // Boundaries of the quantization intervals after the cos(x/pi)
      // transform
      std::vector<double> cos_boundaries(num_quant_vals + 1);
      for (int i = 0; i <= num_quant_vals; ++i) {
        cos_boundaries[i] = std::cos(lsf_boundaries[num_quant_vals - i] * M_PI);
      }
      all_boundaries_.emplace_back(std::move(cos_boundaries));
    }
  }

  const std::vector<double>& boundaries(int p) const {
    return all_boundaries_[p];
  }

 private:
  std::vector<std::vector<double>> all_boundaries_;
};

const LSFQuantizer& GetLSFQuantizer() {
  static const LSFQuantizer* kLsfQuant = new LSFQuantizer();
  return *kLsfQuant;
}

class LSFDequantizer {
 public:
  LSFDequantizer() {
    for (int p = 0; p < kMaxPredictorOrder; ++p) {
      const float q = kLSFQuant[p];
      const int num_quant_vals = std::round(q) + 1;
      std::vector<double> cos_centers(num_quant_vals);
      for (int i = 0; i < num_quant_vals; ++i) {
        cos_centers[i] = std::cos((i / q) * M_PI);
      }
      all_centers_.emplace_back(std::move(cos_centers));
    }
  }

  double Dequantize(int p, int qval) const {
    CHECK_LT(p, all_centers_.size());
    CHECK_LT(qval, all_centers_[p].size());
    return all_centers_[p][qval];
  }

 private:
  std::vector<std::vector<double>> all_centers_;
};

const LSFDequantizer& GetLSFDequantizer() {
  static const LSFDequantizer* kLsfDeq = new LSFDequantizer();
  return *kLsfDeq;
}

bool ComputeLineSpectralFrequencies(const float* pcoefs, uint16_t* quant_lsf,
                                    int order) {
  // The implementation follows the algorithm in the following paper:
  // P. Kabal, R. P. Ramachandran, "The Computation of Line Spectral
  // Frequencies Using Chebyshev Polynomials", IEEE Trans. Acoustics, Speech,
  // Signal Processing, vol. 34, no. 6, pp. 1419–1426
  double A[kMaxPredictorOrder + 2];
  A[0] = 1.0;
  A[order + 1] = 0;
  for (int i = 1; i <= order; ++i) {
    A[i] = -pcoefs[i - 1];
  }
  double F1[kMaxPredictorOrder + 2];
  double F2[kMaxPredictorOrder + 2];
  for (int i = 0; i <= order + 1; ++i) {
    F1[i] = A[i] + A[order + 1 - i];
    F2[i] = A[i] - A[order + 1 - i];
  }
  CHECK_EQ(order % 2, 0);
  double G1[kMaxPredictorOrder + 2];
  double G2[kMaxPredictorOrder + 2];
  G1[order + 1] = 0;
  G2[order + 1] = 0;
  for (int i = order; i >= 0; --i) {
    G1[i] = F1[i + 1] - G1[i + 1];
    G2[i] = -F2[i + 1] + G2[i + 1];
  }
  const int M = order / 2;
  constexpr static int kMaxM = kMaxPredictorOrder / 2 + 1;
  double cheb[2][kMaxM];
  for (int i = 0; i <= M; ++i) {
    cheb[0][i] = (i == 0 ? 1.0 : 2.0) * G1[M - i];
    cheb[1][i] = (i == 0 ? 1.0 : 2.0) * G2[M - i];
  }
  // The alternating roots of cheb[0] and cheb[1] in decreasing order.
  double roots[kMaxPredictorOrder];
  // First step is to find the M roots of cheb[1]. To do this, we evaluate the
  // polynomial on n + 1 points on the [-1.0, 1.0] interval and count the
  // number of times the sign changes. If we find M sign changes, we know an
  // interval for each root, otherwise we double n and repeat.
  for (int j = 1; j >= 0; --j) {
    int cheb_root_pos[kMaxM];
    int n = 1;
    while (n < M) {
      n *= 2;
    }
    std::vector<std::pair<double, double>> values(n + 1);
    bool found_all_roots = false;
    for (int retry = 0; !found_all_roots && retry < 10; ++retry) {
      const double len = 1.0 / n;
      for (int i = 0; i <= n; ++i) {
        if (retry == 0 || (i % 2 == 1)) {
          values[i].first = std::cos(i * len * M_PI);
          values[i].second = EvaluateChebyshevPolynomial(
              absl::Span<const double>(&cheb[j][0], M + 1), values[i].first);
        }
      }
      int ri = 0;
      for (int i = 0; i < n; ++i) {
        if (values[i].second == 0.0 ||
            values[i].second * values[i + 1].second < 0) {
          cheb_root_pos[ri++] = i;
        }
      }
      if (ri == M) {
        found_all_roots = true;
      } else {
        n *= 2;
        values.resize(n + 1);
        for (int i = n; i > 0; i -= 2) {
          values[i] = values[i / 2];
        }
      }
    }
    if (!found_all_roots) {
      return false;
    }
    const auto& lsfq = GetLSFQuantizer();
    for (int ri = 0; ri < M; ++ri) {
      const int i = cheb_root_pos[ri];
      const int p = 2 * ri + j;
      roots[p] = FindRootOfChebyshevPolynomial(
          absl::Span<const double>(&cheb[j][0], M + 1), values[i + 1].first,
          values[i].first, values[i + 1].second, values[i].second,
          lsfq.boundaries(p), &quant_lsf[p]);
    }
    if (j == 0) break;
    // To find the roots of cheb[0] we use the fact that they are bracketed
    // by the roots of cheb[1]. Here we have to check that this is actually
    // the case because we only computed the roots of cheb[1] approximately.
    bool roots_separated = true;
    for (int i = 0; i <= M; ++i) {
      values[i].first = i == 0 ? 1.0 : roots[2 * i - 1];
      values[i].second = EvaluateChebyshevPolynomial(
          absl::Span<const double>(&cheb[0][0], M + 1), values[i].first);
      if (i > 0 && values[i].second * values[i - 1].second >= 0) {
        roots_separated = false;
      }
    }
    if (!roots_separated) {
      // If the approximate roots of cheb[1] did not separate all the
      // roots of cheb[0], then we compute the roots of cheb[0] with the
      // same method as before.
      continue;
    }
    for (int ri = 0; ri < M; ++ri) {
      const int p = 2 * ri;
      const double xa = roots[p + 1];
      const double xb = ri > 0 ? roots[p - 1] : 1.0;
      const double vala = values[ri + 1].second;
      const double valb = values[ri].second;
      roots[p] = FindRootOfChebyshevPolynomial(
          absl::Span<const double>(&cheb[0][0], M + 1), xa, xb, vala, valb,
          lsfq.boundaries(p), &quant_lsf[p]);
    }
    break;
  }
  return true;
}

}  // namespace

void ComputeLinearPredictorCoeffs(const uint16_t* quant_lsf, float* pcoefs,
                                  int order) {
  CHECK_EQ(order % 2, 0);
  const int M = order / 2;
  float G1[kMaxPredictorOrder + 2] = {1.0};
  float G2[kMaxPredictorOrder + 2] = {1.0};
  const auto& lsfdeq = GetLSFDequantizer();
  for (int i = 0; i < M; ++i) {
    const int d = 2 * i;
    const float v1 = -2 * lsfdeq.Dequantize(d, quant_lsf[d]);
    const float v2 = -2 * lsfdeq.Dequantize(d + 1, quant_lsf[d + 1]);
    G1[d + 2] = G1[d];
    G2[d + 2] = G2[d];
    for (int j = d + 1; j > 1; --j) {
      G1[j] += v1 * G1[j - 1] + G1[j - 2];
      G2[j] += v2 * G2[j - 1] + G2[j - 2];
    }
    G1[1] += v1 * G1[0];
    G2[1] += v2 * G2[0];
  }
  for (int i = 0; i < order; ++i) {
    pcoefs[i] = -0.5 * (G1[i] + G1[i + 1] - G2[i] + G2[i + 1]);
  }
}

void ComputePredictorParams(RingliPredictiveHeader* header, float* pcoefs,
                            int order) {
  CHECK(IsInversePredictionFilterStable(pcoefs, order));
  header->quant_lsf.resize(order);
  if (!ComputeLineSpectralFrequencies(pcoefs, &header->quant_lsf[0], order)) {
    fprintf(stderr, "Could not find LSF, reverting to default predictor\n");
    DefaultLineSpectralFrequencies(&header->quant_lsf[0], order);
  }
  ComputeLinearPredictorCoeffs(&header->quant_lsf[0], pcoefs, order);
}

void DefaultLineSpectralFrequencies(uint16_t* quant_lsf, int order) {
  const float step = 1.0 / order;
  for (int i = 0; i < order; ++i) {
    quant_lsf[i] = std::round(i * step * kLSFQuant[i]);
  }
}

}  // namespace ringli
