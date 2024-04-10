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

#include "analysis/opus_codec.h"

#include "gtest/gtest.h"

namespace ringli {
namespace {

TEST(OpusCodecTest, NameIsRight) {
  StreamingOpusCodec codec;

  EXPECT_EQ(codec.Name(), "opus");
}

TEST(OpusCodecTest, CanParseParams) {
  StreamingOpusCodec codec;

  EXPECT_TRUE(codec.ParseParams({"opus", "br256"}));
  EXPECT_EQ(codec.ToString(), "opus:br256");
}

TEST(OpusCodecTest, IgnoresBadParam) {
  StreamingOpusCodec codec;

  EXPECT_FALSE(codec.ParseParams({"opus", "bad_param"}));
}

}  // namespace
}  // namespace ringli
