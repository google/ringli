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

#include <filesystem>
#include <fstream>
#include <string>

#include "absl/flags/flag.h"
#include "absl/log/check.h"
#include "absl/strings/substitute.h"
#include "analysis/ffmpeg_codec.h"
#include "analysis/generate_wav.h"
#include "analysis/opus_codec.h"
#include "common/error_norm.h"
#include "gtest/gtest.h"

namespace ringli {
namespace {

class FFMPEGOpusCodec : public FFMPEGCodec {
 public:
  std::string Name() const override { return "ffmpeg_opus"; }

 protected:
  bool ParseParam(const std::string& param) override {
    if (param.compare(0, 2, "br") == 0) {
      bitrate_ = std::stod(param.substr(2));
    } else {
      return false;
    }
    return true;
  }

  std::vector<std::string> Options() const override {
    return {"-c:a", "libopus", "-b:a", absl::StrCat(bitrate_, "k"),
            "-f",   "opus"};
  }

  std::string ParamsToString() const override {
    return absl::Substitute("br$0", bitrate_);
  }

 private:
  int bitrate_ = 128;
};

class OpusCodecCompatibilityTest : public ::testing::Test {
 public:
  OpusCodecCompatibilityTest() {
    ffmpeg_codec_.ParseParams({ffmpeg_codec_.Name(), "br128"});
    libopus_codec_.ParseParams({libopus_codec_.Name(), "br128"});
    input_ =
        GenerateWav({{.frequency = 200, .amplitude = 0.2}}, 48000.0, 5.0, 0.2);
  }

 protected:
  FFMPEGOpusCodec ffmpeg_codec_;
  StreamingOpusCodec libopus_codec_;
  std::string input_;
};

TEST_F(OpusCodecCompatibilityTest, VeryfyOpusEncoder) {
  std::string libopus_data;
  ASSERT_TRUE(libopus_codec_.Compress(input_, &libopus_data));

  // Verify that ffmpeg can decode what the libopus encoder generated.
  std::string roundtrip1;
  ASSERT_TRUE(ffmpeg_codec_.Decompress(libopus_data, &roundtrip1));

  // Verify that the libopus compressed size is about the same as the ffmpeg
  // compressed size.
  std::string ffmpeg_opus_data;
  ASSERT_TRUE(ffmpeg_codec_.Compress(input_, &ffmpeg_opus_data));
  EXPECT_NEAR(libopus_data.size(), ffmpeg_opus_data.size(), 300);

  // Verify that the libopus compressed file has about the same PSNR as the
  // ffmpeg compressed file.
  std::string roundtrip2;
  ASSERT_TRUE(ffmpeg_codec_.Decompress(ffmpeg_opus_data, &roundtrip2));
  ErrorNorm error_norm1 = CompareFiles(input_, roundtrip1);
  ErrorNorm error_norm2 = CompareFiles(input_, roundtrip2);
  EXPECT_NEAR(error_norm1.psnr, error_norm2.psnr, 0.1);
}

TEST_F(OpusCodecCompatibilityTest, VeryfyOpusDecoder) {
  std::string ffmpeg_opus_data;
  CHECK(ffmpeg_codec_.Compress(input_, &ffmpeg_opus_data));

  // Verify that libopus can decode what the ffmpeg encoder generated.
  std::string roundtrip1;
  CHECK(libopus_codec_.Decompress(ffmpeg_opus_data, &roundtrip1));

  // Verify that libopus decoder produces almost the same PSNR as the
  // ffmpeg decoder.
  std::string roundtrip2;
  ASSERT_TRUE(ffmpeg_codec_.Decompress(ffmpeg_opus_data, &roundtrip2));
  ErrorNorm error_norm1 = CompareFiles(input_, roundtrip1);
  ErrorNorm error_norm2 = CompareFiles(input_, roundtrip2);
  EXPECT_NEAR(error_norm1.psnr, error_norm2.psnr, 0.1);
}

}  // namespace
}  // namespace ringli
