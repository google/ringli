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

#ifndef ANALYSIS_FFMPEG_CODEC_H_
#define ANALYSIS_FFMPEG_CODEC_H_

#include <fcntl.h>

#include <cstdlib>
#include <string>

#include "absl/strings/substitute.h"
#include "analysis/audio_codec.h"
#include "analysis/subprocess.h"

namespace ringli {

class FFMPEGCodec : public AudioCodec {
 public:
  virtual ~FFMPEGCodec() = default;
  bool Compress(const std::string& wav_data,
                std::string* compressed_data) override {
    std::vector<std::string> arguments = {"-i", "-", "-"};
    const std::vector<std::string> opts = Options();
    arguments.insert(arguments.end() - 1, opts.begin(), opts.end());
    const SubprocessResult result = Execute("ffmpeg", arguments, wav_data);
    *compressed_data = result.stdout;
    if (result.status != 0) {
      fprintf(stderr, "%s\n", result.stderr.c_str());
    }
    return result.status == 0;
  }
  bool Decompress(const std::string& compressed_data,
                  std::string* wav_data) override {
    const SubprocessResult result =
        Execute("ffmpeg", {"-i", "-", "-c:a", "pcm_s16le", "-f", "wav", "-"},
                compressed_data);
    *wav_data = result.stdout;
    if (result.status != 0) {
      fprintf(stderr, "%s\n", result.stderr.c_str());
    }
    return result.status == 0;
  }

  virtual std::vector<std::string> Options() const = 0;
};

}  // namespace ringli

#endif  // ANALYSIS_FFMPEG_CODEC_H_
