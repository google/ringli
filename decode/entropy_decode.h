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

#ifndef DECODE_ENTROPY_DECODE_H_
#define DECODE_ENTROPY_DECODE_H_

#include <stddef.h>

#include <vector>

#include "common/context.h"
#include "common/data_defs/constants.h"
#include "common/distributions.h"
#include "common/entropy_coding.h"
#include "common/ringli_header.h"
#include "decode/arith_decode.h"

namespace ringli {

class IntegerArithmeticDecoder {
 public:
  typedef bool (*ProcessOutput)(void* opaque, int val);

  IntegerArithmeticDecoder(int ndirect, int max_sym, void* opaque,
                           ProcessOutput output_cb);

  void set_distribution(Prob* p) { distribution_ = p; }

  bool ProcessInput(uint16_t next_word);

 private:
  bool Output(int value);

  const int ndirect_absval_;
  const int ndirect_symbols_;
  const int max_symbols_;
  void* const opaque_;
  ProcessOutput const output_cb_;
  BinaryArithmeticDecoder ac_;
  enum { SYMBOL_DECODING, EXTRA_BITS_DECODING } state_;
  int val0_;
  int val1_;
  int sign_;
  int msb_;
  int nbits_;
  int bitpos_;
  int extra_bits_val_;
  Prob* distribution_;
};

class EntropyDecoder {
 public:
  typedef bool (*ProcessSamples)(void* opaque, const int* samples);
  EntropyDecoder(size_t num_channels, void* opaque,
                 ProcessSamples process_samples);

  bool ProcessInput(const uint8_t* data, size_t len);

 private:
  bool ProcessOutput(int value);
  static bool ProcessOutputCb(void* opaque, int value) {
    return reinterpret_cast<EntropyDecoder*>(opaque)->ProcessOutput(value);
  }
  void SetContext();

  const size_t num_channels_;
  void* const opaque_;
  ProcessSamples const process_samples_;
  IntegerArithmeticDecoder int_decoder_;
  std::vector<PredictiveContextModel> context_model_;
  std::vector<Prob> symbol_prob_;
  uint16_t next_word_;
  int input_shift_;
  std::vector<int> samples_;
  size_t channel_idx_;
  size_t idx_;
};

bool DecompressPredictiveRingliBlocks(const char* input, size_t input_size,
                                      size_t num_channels, size_t num_blocks,
                                      const RingliDecoderConfig& config,
                                      std::vector<RingliBlock>* ringli_blocks);

bool DecompressCoefficients(const char* input, size_t input_size,
                            size_t num_channels, size_t num_blocks,
                            const RingliDecoderConfig& config,
                            std::vector<RingliBlock>* ringli_blocks);

}  // namespace ringli

#endif  // DECODE_ENTROPY_DECODE_H_
