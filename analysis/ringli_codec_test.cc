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

#include "analysis/ringli_codec.h"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/log/check.h"
#include "absl/strings/str_split.h"
#include "absl/strings/string_view.h"
#include "analysis/generate_wav.h"
#include "common/error_norm.h"
#include "gtest/gtest.h"

namespace ringli {
namespace {

TEST(RingliCodecTest, NameIsRight) {
  StreamingRingliCodec codec;

  EXPECT_EQ(codec.Name(), "ringli");
}

TEST(RingliCodecTest, DefaultParams) {
  StreamingRingliCodec codec;

  EXPECT_EQ(codec.ParamsToString(), "qc(0;1)");
}

struct RingliTestParams {
  std::string codec_params;
};

class RingliCodecParamTest
    : public ::testing::Test,
      public testing::WithParamInterface<RingliTestParams> {};

INSTANTIATE_TEST_SUITE_P(
    RingliParseParams, RingliCodecParamTest,
    testing::Values(RingliTestParams{"ringli:qc(0;7)"},
                    RingliTestParams{"ringli:pc:o4-28:e7:q7"},
                    RingliTestParams{"ringli:apc:e7:q3"},
                    RingliTestParams{"ringli:aconly:qc(0;7)"}));

TEST_P(RingliCodecParamTest, CanParseParams) {
  StreamingRingliCodec codec;

  const std::string codec_params_string = GetParam().codec_params;
  const std::vector<std::string> codec_params =
      absl::StrSplit(codec_params_string, ':');

  EXPECT_TRUE(codec.ParseParams(codec_params));

  EXPECT_EQ(codec.ToString(), codec_params_string);
}

struct RingliEvaluationTestParams {
  std::string codec_params = "";
  int64_t compressed_size = 0;
  double psnr = 0.0;
};

class RingliCodecEvaluationTest
    : public ::testing::Test,
      public testing::WithParamInterface<RingliEvaluationTestParams> {
 protected:
  explicit RingliCodecEvaluationTest(const std::vector<Waveform>& waveforms,
                                     float noise) {
    const std::string codec_params_string = GetParam().codec_params;
    const std::vector<std::string> codec_params =
        absl::StrSplit(codec_params_string, ':');

    EXPECT_TRUE(codec_.ParseParams(codec_params));

    input_ = GenerateWav(waveforms, 48000.0, 5.0, noise);
  }

  void CheckCompression() {
    EXPECT_TRUE(codec_.Compress(input_, &compressed_));
    EXPECT_NEAR(compressed_.size(), GetParam().compressed_size, 1000)
        << "*** codec_params=" << GetParam().codec_params;
  }

  void CheckDecompression() {
    EXPECT_TRUE(codec_.Decompress(compressed_, &decompressed_));
  }
  void CheckPsnr() {
    const ErrorNorm error_norm = CompareFiles(input_, decompressed_);
    const double expected_psnr = GetParam().psnr;
    if (expected_psnr > 0) {
      EXPECT_NEAR(error_norm.psnr, expected_psnr, 1.0)
          << "*** codec_params=" << GetParam().codec_params;
    } else {
      EXPECT_EQ(error_norm.psnr, std::numeric_limits<double>::infinity());
    }
  }

  StreamingRingliCodec codec_;
  std::string input_;
  std::string compressed_;
  std::string decompressed_;
};

class RingliCodecEvaluationSimpleSineTest : public RingliCodecEvaluationTest {
 protected:
  RingliCodecEvaluationSimpleSineTest()
      : RingliCodecEvaluationTest({{.frequency = 150.0, .amplitude = 0.5}}, 0) {
  }
};

INSTANTIATE_TEST_SUITE_P(
    RingliCompressSimpleSine, RingliCodecEvaluationSimpleSineTest,
    testing::Values(
        RingliEvaluationTestParams{"ringli:qc(0;7)", 90261, 84},
        RingliEvaluationTestParams{"ringli:pc:o2-8:e5:q7", 21447, 84},
        RingliEvaluationTestParams{"ringli:pc:o2-8:e5:q1", 61616, -1},
        RingliEvaluationTestParams{"ringli:apc:e7:q3", 30513, 92},
        RingliEvaluationTestParams{"ringli:apc:aconly:e7:q3", 29690, 92},
        RingliEvaluationTestParams{"ringli:aconly:qc(0;7)", 91029, 84}));

TEST_P(RingliCodecEvaluationSimpleSineTest, CompressedSizeAndPsnrWithinRange) {
  CheckCompression();
  CheckDecompression();
  CheckPsnr();
}

class RingliCodecEvaluationMultipleSineTest : public RingliCodecEvaluationTest {
 protected:
  RingliCodecEvaluationMultipleSineTest()
      : RingliCodecEvaluationTest({{.frequency = 150.0, .amplitude = 0.5},
                                   {.frequency = 5100.0, .amplitude = 0.2}},
                                  0) {}
};

INSTANTIATE_TEST_SUITE_P(
    RingliCompressMultipleSine, RingliCodecEvaluationMultipleSineTest,
    testing::Values(
        RingliEvaluationTestParams{"ringli:qb(0;1);(1000;1000)", 33391, 46},
        RingliEvaluationTestParams{"ringli:pc:o2-8:e5:q7", 133794, 87},
        RingliEvaluationTestParams{"ringli:pc:o2-8:e5:q1", 216750, -1},
        RingliEvaluationTestParams{"ringli:apc:e7:q3", 122955, 94},
        RingliEvaluationTestParams{"ringli:apc:aconly:e7:q3", 115776, 94},
        RingliEvaluationTestParams{"ringli:aconly:qc(0;7)", 227521, 86}));

TEST_P(RingliCodecEvaluationMultipleSineTest,
       CompressedSizeAndPsnrWithinRange) {
  CheckCompression();
  CheckDecompression();
  CheckPsnr();
}

class RingliCodecEvaluationNoiseTest : public RingliCodecEvaluationTest {
 protected:
  RingliCodecEvaluationNoiseTest()
      : RingliCodecEvaluationTest({{.frequency = 150.0, .amplitude = 0.5},
                                   {.frequency = 5100.0, .amplitude = 0.2}},
                                  0.05) {}
};

INSTANTIATE_TEST_SUITE_P(
    RingliCompressNoise, RingliCodecEvaluationNoiseTest,
    testing::Values(
        RingliEvaluationTestParams{"ringli:qb(0;1);(1000;1000)", 83799, 44},
        RingliEvaluationTestParams{"ringli:pc:o2-8:e5:q7", 292405, 87},
        RingliEvaluationTestParams{"ringli:pc:o2-8:e5:q1", 377206, -1},
        RingliEvaluationTestParams{"ringli:apc:e7:q3", 330827, 95},
        RingliEvaluationTestParams{"ringli:apc:aconly:e7:q3", 331928, 95},
        RingliEvaluationTestParams{"ringli:aconly:qc(0;7)", 301957, 87}));
// TODO(antoran): add test cases for noise shaping and adaptive quantizaiton.

TEST_P(RingliCodecEvaluationNoiseTest, CompressedSizeAndPsnrWithinRange) {
  CheckCompression();
  CheckDecompression();
  CheckPsnr();
}

}  // namespace
}  // namespace ringli
