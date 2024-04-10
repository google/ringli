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

#ifndef ENCODE_ENTROPY_ENCODE_H_
#define ENCODE_ENTROPY_ENCODE_H_

#include <memory>
#include <string>
#include <vector>

#include "common/context.h"
#include "common/data_defs/constants.h"
#include "common/distributions.h"
#include "common/entropy_coding.h"
#include "common/ringli_header.h"
#include "encode/ans_encode.h"
#include "encode/arith_encode.h"

namespace ringli {

struct Histogram {
  Histogram() { Clear(); }

  void Clear() {
    memset(data, 0, sizeof(data));
    total_count = 0;
  }

  void AddHistogram(const Histogram& other) {
    for (int i = 0; i < MAX_SYMBOLS; ++i) {
      data[i] += other.data[i];
    }
    total_count += other.total_count;
  }

  void Add(int val) {
    ++data[val];
    ++total_count;
  }

  int data[MAX_SYMBOLS];
  int total_count;
  // Thise fields are only updated and used by the clustering functions.
  double bit_cost;
  mutable double entropy;
};

// Manages building, clustering and encoding of the histograms of an entropy
// source.
class EntropySource {
 public:
  void Resize(int num_contexts) { histograms_.resize(num_contexts); }

  void AddCode(int code, int histo_ix) { histograms_[histo_ix].Add(code); }

  void ClusterHistograms();

  void EncodeContextMap(size_t* storage_ix, uint8_t* storage) const;

  void BuildAndStoreEntropyCodes(size_t* storage_ix, uint8_t* storage);

  const ANSTable* GetANSTable(int context) const {
    const int entropy_ix = context_map_[context];
    return &ans_tables_[entropy_ix];
  }

  size_t NumHistograms() const { return clustered_.size(); }

  double ClusteredEntropy(int histo_idx);

 private:
  static constexpr int kMaxNumberOfHistograms = 256;
  std::vector<Histogram> histograms_;
  std::vector<Histogram> clustered_;
  std::vector<uint32_t> context_map_;
  std::vector<ANSTable> ans_tables_;
};

// Manages the multiplexing of the ANS-coded and arithmetic coded bits.
class DataStream {
 public:
  explicit DataStream(EntropySource* entropy_source);

  void Resize(int max_num_code_words) {
    code_words_.resize(max_num_code_words);
  }

  void ResizeForBlock();

  void AddCode(int code, int context);

  void AddBits(int nbits, int bits);

  void FlushArithmeticCoder();

  void FlushBitWriter();

  // Encodes the next bit to the bit stream, based on the 8-bit precision
  // probability, i.e. P(bit = 0) = prob / 256.
  void AddBit(uint8_t prob, int bit);

  // Same as above, but tatistics are also updated in 'p'.
  void AddBit(Prob* p, int bit);

  void EncodeCodeWords(const EntropySource& s,
                       const EntropyCodingParams& ecparams, uint8_t* data,
                       size_t* pos, size_t len);

  size_t TotalExtraBits() const { return total_extra_bits_; }

 private:
  struct CodeWord {
    // Add a constructor that does nothing (unlike the default one) to avoid
    // initializing to values that are unused anyway.
    CodeWord() = default;
    uint32_t context;
    uint16_t value;
    uint8_t code;
    uint8_t nbits;
  };
  static constexpr size_t kSlackForOneBlock =
      6 + 2 * (kMaxPredictorOrder + kRingliBlockSize);

  static void WriteUint16(uint16_t val, uint8_t* data, size_t* pos) {
    data[(*pos)++] = val & 0xff;
    data[(*pos)++] = val >> 8;
  }

  int pos_;
  int bw_pos_;
  int ac_pos0_;
  int ac_pos1_;
  uint32_t low_;
  uint32_t high_;
  uint32_t bw_val_;
  int bw_bitpos_;
  std::vector<CodeWord> code_words_;
  EntropySource* entropy_source_;
  size_t total_extra_bits_;
};

class EntropyCoder {
 public:
  EntropyCoder(const EntropyCodingParams& ecparams, uint32_t sampling_freq,
               uint32_t num_channels, bool predictive, bool online);

  void Reset();

  bool ProcessBlock(const RingliBlock& block, std::string* output);

  void ProcessSamples(const int* samples, std::string* output);

  bool Flush(std::string* output);

 private:
  bool ProcessPredictiveBlock(const RingliBlock& block, std::string* output);

  EntropyCodingParams ecparams_;
  uint32_t sampling_freq_;
  uint32_t num_channels_;
  bool predictive_;
  bool online_;
  BinaryArithmeticEncoder arith_encode_;
  std::unique_ptr<EntropySource> entropy_source_;
  std::unique_ptr<DataStream> data_stream_;
  std::vector<PredictiveContextModel> context_model_;
  std::vector<Prob> symbol_prob_;
  int order_histo_[kMaxPredictorOrder + 1];
  int lsf_extra_bits_;
  uint32_t num_samples_;
  size_t idx_;
};

bool CompressCoefficients(const std::vector<RingliBlock>& ringli_blocks,
                          size_t num_channels,
                          const EntropyCodingParams& ecparams,
                          std::string* output);

}  // namespace ringli

#endif  // ENCODE_ENTROPY_ENCODE_H_
