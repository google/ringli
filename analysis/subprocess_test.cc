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

#include "analysis/subprocess.h"

#include <ctime>
#include <iostream>

#include "gtest/gtest.h"

namespace ringli {

TEST(SubprocessTest, TestExecuteWithoutStdin) {
  SubprocessResult result = Execute("date", {"+%s"}, "");
  EXPECT_EQ(result.status, 0);
  EXPECT_EQ(result.stderr, "");
  size_t got_time;
  sscanf(result.stdout.c_str(), "%zu", &got_time);
  std::time_t expected_time = std::time(0);
  size_t week = 7 * 24 * 3600;
  EXPECT_GT(got_time, expected_time - week);
  EXPECT_LT(got_time, expected_time + week);
}

TEST(SubprocessTest, TestExecuteWithStdin) {
  SubprocessResult result = Execute("cat", {}, "TEST");
  EXPECT_EQ(result.status, 0);
  EXPECT_EQ(result.stderr, "");
  EXPECT_EQ(result.stdout, "TEST");
}

}  // namespace ringli
