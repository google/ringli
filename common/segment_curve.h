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

#ifndef COMMON_SEGMENT_CURVE_H_
#define COMMON_SEGMENT_CURVE_H_

#include <string>
#include <utility>
#include <vector>
namespace ringli {
// A class that contains a set of basis points and returns value f(x) for each
// input following rules:
// - if x <= first_point.x => first_point.y
// - if x >= last_point.x => last_point.y
// Otherwise - linearly extrapolates the value between two points.
class SegmentCurve {
 public:
  SegmentCurve() = default;
  explicit SegmentCurve(const std::vector<std::pair<double, double>>& points);
  explicit SegmentCurve(double const_val);

  // Parse SegmentCurve from input with points. Supported formats:
  // each point is [x, y] or (x, y) or {x, y}.
  // Example: (1.5, 10), (5, 12.5)
  static SegmentCurve ParseFromString(const std::string& input_string);

  std::string ToString() const;

  double GetValue(double x) const;

 private:
  std::vector<std::pair<double, double>> points_;
  friend class SegmentCurveTest;
};

}  // namespace ringli

#endif  // COMMON_SEGMENT_CURVE_H_
