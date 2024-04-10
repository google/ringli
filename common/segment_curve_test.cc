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
#include <utility>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace ringli {

class SegmentCurveTest : public ::testing::Test {
 public:
  const std::vector<std::pair<double, double>>& GetCurvePoints(
      const SegmentCurve& curve) {
    return curve.points_;
  }
};

TEST_F(SegmentCurveTest, EmptytCurveReturnsValueItselt) {
  SegmentCurve curve;
  EXPECT_EQ(curve.GetValue(2), 2);
  EXPECT_EQ(curve.GetValue(5), 5);
}

TEST_F(SegmentCurveTest, OnePointCurveReturnsValues) {
  SegmentCurve curve({{0, 5}});
  EXPECT_EQ(curve.GetValue(2), 5);
  EXPECT_EQ(curve.GetValue(-2), 5);
}

TEST_F(SegmentCurveTest, OnePointCurveFromSingleValue) {
  SegmentCurve curve(5);
  EXPECT_EQ(curve.GetValue(2), 5);
  EXPECT_EQ(curve.GetValue(-2), 5);
}

TEST_F(SegmentCurveTest, MultiplePointsCurveReturnsSmallestValue) {
  SegmentCurve curve({{1, 2}, {4, 5}});
  EXPECT_EQ(curve.GetValue(0), 2);
  EXPECT_EQ(curve.GetValue(1), 2);
}

TEST_F(SegmentCurveTest, MultiplePointsCurveReturnsLargestValue) {
  SegmentCurve curve({{1, 2}, {4, 5}});
  EXPECT_EQ(curve.GetValue(4), 5);
  EXPECT_EQ(curve.GetValue(5), 5);
}

TEST_F(SegmentCurveTest, MultiplePointsCurveExtrapolatesValueBetweenPoints) {
  SegmentCurve curve({{1, 10}, {5, 12}});
  EXPECT_EQ(curve.GetValue(2), 10.5);
  EXPECT_EQ(curve.GetValue(3), 11);
}

TEST_F(SegmentCurveTest, PointsSortedByX) {
  SegmentCurve curve({{5, 10}, {1, 12}});
  EXPECT_TRUE(std::is_sorted(GetCurvePoints(curve).begin(),
                             GetCurvePoints(curve).end()));
}

TEST_F(SegmentCurveTest, PointsDuplicatesRemoved) {
  SegmentCurve curve({{1, 10}, {1, 12}, {2, 14}, {2, 16}});
  EXPECT_THAT(GetCurvePoints(curve), testing::SizeIs(2));
}

TEST_F(SegmentCurveTest, CanParseFromString) {
  SegmentCurve curve = SegmentCurve::ParseFromString("(1.5; 10); (5; 12.5)");
  SegmentCurve curve2 = SegmentCurve::ParseFromString("[1.5; 10]; [5; 12.5]");
  SegmentCurve curve3 = SegmentCurve::ParseFromString("{1.5; 10}; {5; 12.5}");
  EXPECT_THAT(GetCurvePoints(curve),
              testing::ElementsAre(std::pair<double, double>{1.5, 10},
                                   std::pair<double, double>{5, 12.5}));
  EXPECT_THAT(GetCurvePoints(curve2), GetCurvePoints(curve));
  EXPECT_THAT(GetCurvePoints(curve2), GetCurvePoints(curve3));
}

TEST_F(SegmentCurveTest, CanSerializeToString) {
  SegmentCurve curve({{1, 10}, {2, 14}});
  EXPECT_EQ(curve.ToString(), "(1;10);(2;14)");
}

TEST_F(SegmentCurveTest, CanParseFromSerialized) {
  SegmentCurve curve({{1, 10}, {2, 14}});
  EXPECT_EQ(GetCurvePoints(curve),
            GetCurvePoints(SegmentCurve::ParseFromString(curve.ToString())));
}

}  // namespace ringli
