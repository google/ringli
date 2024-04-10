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

#include "analysis/warp_q.h"

#include <cmath>
#include <iostream>
#include <regex>

#include "analysis/subprocess.h"

namespace ringli {

namespace {

std::vector<std::string> Split(const std::string& s, char sep) {
  std::vector<std::string> result;
  std::stringstream lines(s);
  std::string line;
  while (std::getline(lines, line, sep)) {
    result.push_back(line);
  }
  return result;
}

}  // namespace

float WarpQ::MOS(const std::filesystem::path& reference,
                 const std::filesystem::path& degraded) {
  SubprocessResult result = Execute(warp_q_binary, {reference, degraded}, "");
  if (result.status != 0) {
    std::cerr << "Unable to execute '" << warp_q_binary
              << "': " << result.status << std::endl
              << "STDOUT:" << std::endl
              << result.stdout << std::endl
              << "STDERR:" << std::endl
              << result.stderr << std::endl;
  }
  for (const auto& line : Split(result.stdout, '\n')) {
    const size_t pos =
        line.find("Mapped WARP-Q score (higher rating means better quality): ");
    if (pos != std::string::npos) {
      return std::stof(Split(line, ' ').back());
    }
  }
  std::cerr << "failed to find Mapped WARP-Q score in" << std::endl
            << result.stdout << std::endl;
  return std::numeric_limits<float>::quiet_NaN();
}

}  // namespace ringli
