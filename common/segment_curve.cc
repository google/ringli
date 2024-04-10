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

#include "common/segment_curve.h"

#include <algorithm>
#include <regex>  // NOLINT
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace ringli {

SegmentCurve::SegmentCurve(
    const std::vector<std::pair<double, double>>& points) {
  points_ = points;
  std::sort(points_.begin(), points_.end());
  if (points_.size() > 1) {
    auto it = points_.begin();
    double prev_x = it->first;
    it++;
    while (it != points_.end()) {
      double cur_x = it->first;
      if (cur_x == prev_x) {
        it = points_.erase(it);
      } else {
        it++;
      }
      prev_x = cur_x;
    }
  }
}

SegmentCurve::SegmentCurve(double const_val) {
  points_.push_back(std::make_pair(/* does not matter */ 0, const_val));
}

SegmentCurve SegmentCurve::ParseFromString(const std::string& input_string) {
  std::regex point_regex(
      "[\\(\\{\\[]([\\d\\.]+)\\s*;\\s*([\\d\\.]+)[\\)\\}\\]]");
  SegmentCurve curve;
  for (std::sregex_iterator it = std::sregex_iterator(
           input_string.begin(), input_string.end(), point_regex);
       it != std::sregex_iterator(); ++it) {
    std::smatch match = *it;
    curve.points_.push_back(
        std::make_pair(std::stod(match[1]), std::stod(match[2])));
  }
  return curve;
}

std::string SegmentCurve::ToString() const {
  std::stringstream ss;
  for (int i = 0; i < points_.size(); ++i) {
    ss << "(" << points_[i].first << ";" << points_[i].second << ")";
    if (i < points_.size() - 1) {
      ss << ";";
    }
  }
  return ss.str();
}

double SegmentCurve::GetValue(double x) const {
  if (points_.empty()) {
    return x;
  }
  if (points_.size() == 1) {
    return points_[0].second;
  }
  if (x <= points_[0].first) {
    return points_[0].second;
  }
  for (int i = 1; i < points_.size(); ++i) {
    const auto& p1 = points_[i - 1];
    const auto& p2 = points_[i];
    if (x <= p2.first) {
      const double a = (p2.second - p1.second) / (p2.first - p1.first);
      const double b = p1.second - a * p1.first;
      return a * x + b;
    }
  }
  return points_.back().second;
}

}  // namespace ringli
