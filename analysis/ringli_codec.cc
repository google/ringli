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

#include <string>
#include <vector>

#include "absl/strings/str_cat.h"
#include "absl/strings/str_join.h"
#include "absl/strings/str_split.h"
#include "absl/strings/substitute.h"
#include "common/ringli_header.h"
#include "common/segment_curve.h"
#include "encode/ringli_encoder.h"

namespace ringli {

namespace {

bool ParseParam(const std::string& param, RingliEncoderConfig& config) {
  if (param[0] == 'q') {
    if (config.dconfig.use_predictive_coding) {
      config.dconfig.pred_quant = std::stoi(param.substr(1));
    } else {
      if (param[1] == 'c') {
        config.quantization_type = QuantizationType::CONST;
      } else if (param[1] == 'p') {
        config.quantization_type = QuantizationType::ADAPTIVE_PER_PACKET;
      } else if (param[1] == 'b') {
        config.quantization_type = QuantizationType::ADAPTIVE_PER_BAND;
      } else {
        return false;
      }
      config.quantization_curve =
          SegmentCurve::ParseFromString(param.substr(2));
    }
  } else if (param[0] == 'o') {
    if (config.dconfig.use_predictive_coding &&
        !config.dconfig.use_online_predictive_coding) {
      const std::vector<std::string>& v = absl::StrSplit(param.substr(1), '-');
      config.pred_order_min = std::stoi(v[0]);
      config.pred_order_max =
          v.size() > 1 ? std::stoi(v[1]) : config.pred_order_min;
    } else {
      return false;
    }
  } else if (param == "aq") {
    config.dconfig.pred_quant = 1;
    config.dconfig.use_adaptive_quantization = true;
  } else if (param[0] == 'e') {
    config.dconfig.effort = std::stoi(param.substr(1));
  } else if (param == "pc") {
    config.dconfig.use_predictive_coding = 1;
  } else if (param == "apc") {
    config.dconfig.use_predictive_coding = 1;
    config.dconfig.use_online_predictive_coding = 1;
  } else if (param == "aconly") {
    config.dconfig.ecparams.arithmetic_only = 1;
  } else if (param == "ns") {
    config.use_noise_shaping = true;
  } else if (param == "nf") {
    config.dconfig.use_noise_filter = true;
  } else {
    return false;
  }
  return true;
}

std::string ConfigToString(const RingliEncoderConfig& config) {
  std::vector<std::string> result;
  if (config.dconfig.use_predictive_coding) {
    if (config.dconfig.use_online_predictive_coding) {
      result.push_back("apc");
    } else {
      result.push_back("pc");
      if (config.pred_order_min == config.pred_order_max) {
        result.push_back(absl::Substitute("o$0", config.pred_order_min));
      } else {
        result.push_back(absl::Substitute("o$0-$1", config.pred_order_min,
                                          config.pred_order_max));
      }
    }
    if (config.dconfig.ecparams.arithmetic_only) {
      result.push_back("aconly");
    }
    result.push_back(absl::Substitute("e$0", config.dconfig.effort));
    result.push_back(absl::Substitute("q$0", config.dconfig.pred_quant));
    if (config.use_noise_shaping) {
      result.push_back("ns");
    }
    if (config.dconfig.use_noise_filter) {
      result.push_back("nf");
    }
    if (config.dconfig.use_adaptive_quantization) {
      result.push_back("aq");
    }
  } else if (config.dconfig.use_predictive_coding &&
             config.dconfig.use_online_predictive_coding) {
    result.push_back(absl::Substitute("q$0", config.dconfig.pred_quant));
  } else {
    if (config.dconfig.ecparams.arithmetic_only) {
      result.push_back("aconly");
    }
    std::string quantization_type;
    switch (config.quantization_type) {
      case QuantizationType::CONST:
        quantization_type = 'c';
        break;
      case QuantizationType::ADAPTIVE_PER_PACKET:
        quantization_type = 'p';
        break;
      case QuantizationType::ADAPTIVE_PER_BAND:
        quantization_type = 'b';
        break;
    }
    result.push_back(absl::StrCat("q", quantization_type,
                                  config.quantization_curve.ToString()));
  }
  return absl::StrJoin(result, ":");
}

}  // namespace

bool StreamingRingliCodec::ParseParam(const std::string& param) {
  return ::ringli::ParseParam(param, config_);
}

std::string StreamingRingliCodec::ParamsToString() const {
  return ConfigToString(config_);
}

}  // namespace ringli
