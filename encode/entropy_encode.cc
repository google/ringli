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

#include "encode/entropy_encode.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/log/check.h"
#include "absl/types/span.h"
#include "common/context.h"
#include "common/data_defs/constants.h"
#include "common/distributions.h"
#include "common/entropy_coding.h"
#include "common/log2floor.h"
#include "common/logging.h"
#include "common/ringli_header.h"
#include "encode/ans_encode.h"
#include "encode/arith_encode.h"
#include "encode/cluster.h"
#include "encode/context_map_encode.h"
#include "encode/histogram_encode.h"
#include "encode/write_bits.h"

namespace ringli {

double PopulationCost(const Histogram& h) {
  return PopulationCost(&h.data[0], h.total_count);
}

double HistogramCrossEntropy(const Histogram& a, const Histogram& b) {
  double entropy = 0.0;
  for (int i = 0; i < MAX_SYMBOLS; ++i) {
    if (b.data[i]) {
      entropy -= a.data[i] * std::log(b.data[i]);
    }
  }
  entropy += a.total_count * std::log(b.total_count);
  entropy /= std::log(2.0);
  return entropy;
}

double HistogramEntropy(const Histogram& h) {
  return HistogramCrossEntropy(h, h);
}

double HistogramDistance(const Histogram& a, const Histogram& b) {
  Histogram c;
  c.AddHistogram(a);
  c.AddHistogram(b);
  return HistogramEntropy(c) - a.entropy - b.entropy;
}

void EntropySource::ClusterHistograms() {
  ::ringli::ClusterHistograms(histograms_, 1, histograms_.size(),
                              std::vector<int>(), kMaxNumberOfHistograms,
                              &clustered_, &context_map_);
}

void EntropySource::EncodeContextMap(size_t* storage_ix,
                                     uint8_t* storage) const {
  ::ringli::EncodeContextMap(context_map_, clustered_.size(), storage_ix,
                             storage);
}

void EntropySource::BuildAndStoreEntropyCodes(size_t* storage_ix,
                                              uint8_t* storage) {
  ans_tables_.resize(clustered_.size());
  for (size_t i = 0; i < clustered_.size(); ++i) {
    BuildAndStoreANSEncodingData(&clustered_[i].data[0], &ans_tables_[i],
                                 storage_ix, storage);
  }
}

double EntropySource::ClusteredEntropy(int histo_idx) {
  const Histogram& histo_a = histograms_[histo_idx];
  const Histogram& histo_b = clustered_[context_map_[histo_idx]];
  return HistogramCrossEntropy(histo_a, histo_b);
}

DataStream::DataStream(EntropySource* entropy_source)
    : pos_(3),
      bw_pos_(0),
      ac_pos0_(1),
      ac_pos1_(2),
      low_(0),
      high_(~0),
      bw_val_(0),
      bw_bitpos_(0),
      entropy_source_(entropy_source),
      total_extra_bits_(0) {}

void DataStream::ResizeForBlock() {
  if (pos_ + kSlackForOneBlock > code_words_.size()) {
    static const double kGrowMult = 1.2;
    const size_t new_size =
        kGrowMult * code_words_.capacity() + kSlackForOneBlock;
    code_words_.resize(new_size);
  }
}

void DataStream::AddCode(int code, int context) {
  CodeWord word;
  word.context = context;
  word.code = code;
  word.nbits = 0;
  word.value = 0;
  CHECK(pos_ < code_words_.size());
  code_words_[pos_++] = word;
  entropy_source_->AddCode(code, context);
}

void DataStream::AddBits(int nbits, int bits) {
  bw_val_ |= (bits << bw_bitpos_);
  bw_bitpos_ += nbits;
  total_extra_bits_ += nbits;
  if (bw_bitpos_ > 16) {
    CodeWord word;
    word.context = 0;
    word.code = 0;
    word.nbits = 16;
    word.value = bw_val_ & 0xffff;
    code_words_[bw_pos_] = word;
    bw_pos_ = pos_;
    ++pos_;
    bw_val_ >>= 16;
    bw_bitpos_ -= 16;
  }
}

void DataStream::FlushArithmeticCoder() {
  code_words_[ac_pos0_].value = high_ >> 16;
  code_words_[ac_pos1_].value = high_ & 0xffff;
  code_words_[ac_pos0_].nbits = 16;
  code_words_[ac_pos1_].nbits = 16;
  low_ = 0;
  high_ = ~0;
}

void DataStream::FlushBitWriter() {
  code_words_[bw_pos_].nbits = 16;
  code_words_[bw_pos_].value = bw_val_ & 0xffff;
}

void DataStream::AddBit(Prob* const p, int bit) {
  const uint8_t prob = p->get_proba();
  p->Add(bit);
  AddBit(prob, bit);
}

void DataStream::AddBit(uint8_t prob, int bit) {
  while (((low_ ^ high_) >> 16) == 0) {
    code_words_[ac_pos0_].value = high_ >> 16;
    code_words_[ac_pos0_].nbits = 16;
    ac_pos0_ = ac_pos1_;
    ac_pos1_ = pos_;
    ++pos_;
    low_ <<= 16;
    high_ <<= 16;
    high_ |= 0xffff;
  }
  const uint32_t diff = high_ - low_;
  const uint32_t split = low_ + (((uint64_t)diff * prob) >> 8);
  if (bit) {
    low_ = split + 1;
  } else {
    high_ = split;
  }
}

void DataStream::EncodeCodeWords(const EntropySource& s,
                                 const EntropyCodingParams& ecparams,
                                 uint8_t* data, size_t* pos, size_t len) {
  FlushBitWriter();
  FlushArithmeticCoder();
  if (!ecparams.arithmetic_only) {
    ANSCoder ans;
    for (int i = pos_ - 1; i >= 0; --i) {
      CodeWord* const word = &code_words_[i];
      if (word->nbits == 0) {
        const ANSEncSymbolInfo info =
            s.GetANSTable(word->context)->info_[word->code];
        word->value = ans.PutSymbol(info, &word->nbits);
      }
    }
    const uint32_t state = ans.GetState();
    CHECK(*pos + 4 <= len);
    WriteUint16((state >> 16) & 0xffff, data, pos);
    WriteUint16((state >> 0) & 0xffff, data, pos);
  }
  for (int i = 0; i < pos_; ++i) {
    const CodeWord& word = code_words_[i];
    if (word.nbits) {
      CHECK(*pos + 2 <= len);
      WriteUint16(word.value, data, pos);
    }
  }
}

void EncodeSymbol(int val, int alphabet_size, Prob* p,
                  DataStream* data_stream) {
  int val0 = 0;
  int val1 = alphabet_size;
  while (val0 + 1 < val1) {
    int mid = (val0 + val1) >> 1;
    int bit = (val >= mid);
    data_stream->AddBit(&p[mid - 1], bit);
    if (bit) {
      val0 = mid;
    } else {
      val1 = mid;
    }
  }
}

size_t Base128Size(size_t val) {
  size_t size = 1;
  for (; val >= 128; val >>= 7) ++size;
  return size;
}

void EncodeBase128Fix(size_t val, size_t len, uint8_t* data) {
  for (size_t i = 0; i < len; ++i) {
    *data++ = (val & 0x7f) | (i + 1 < len ? 0x80 : 0);
    val >>= 7;
  }
}

uint16_t ConvertToUnsigned(int16_t val) {
  int x = val;
  int ux = x >= 0 ? 2 * x : -2 * x - 1;
  DCHECK_GE(ux, 0);
  DCHECK_LE(ux, 0xffff);
  return ux;
}

int EncodeValue(int value, int ndirect, int* nbits, int* extra_bits) {
  const int absval = std::abs(value);
  if (absval < ndirect) {
    return ConvertToUnsigned(value);
  } else {
    const int n = Log2FloorNonZero(absval - ndirect + 2);
    const int m = absval - ndirect + 2 - (1 << n);
    *nbits = n - 1;
    const int msb = m >> (*nbits);
    *extra_bits = m & ((1 << (*nbits)) - 1);
    const int sign = value >= 0 ? 0 : 1;
    return 2 * ndirect - 1 + 4 * (*nbits) + 2 * msb + sign;
  }
}

void ProcessCoefficients(absl::Span<const RingliBlock> ringli_blocks,
                         size_t num_channels,
                         const EntropyCodingParams& ecparams,
                         EntropySource* entropy_source,
                         DataStream* data_stream) {
  const size_t num_blocks = ringli_blocks.size();
  const size_t num_contexts = 2 + kNumZeroDensityContexts;
  entropy_source->Resize(num_contexts);

  std::vector<Prob> last_nz_prob(kDctLength - 1);
  std::vector<Prob> is_zero_prob(kNumZeronessContexts);
  std::vector<Prob> sign_prob(kDctLength);
  std::vector<Prob> symbol_prob;
  if (ecparams.arithmetic_only) {
    symbol_prob.resize(num_contexts * (MAX_SYMBOLS - 1));
  }

  for (size_t i = 0; i < num_blocks; ++i) {
    data_stream->ResizeForBlock();
    const auto& header = ringli_blocks[i].header.dct;
    for (int band = 0; band < kNumDctBands; ++band) {
      uint8_t quant_msb = header.quant[band] >> 8;
      uint8_t quant_lsb = header.quant[band] & 0xff;
      if (ecparams.arithmetic_only) {
        EncodeSymbol(quant_msb, MAX_SYMBOLS, &symbol_prob[0], data_stream);
        EncodeSymbol(quant_lsb, MAX_SYMBOLS, &symbol_prob[MAX_SYMBOLS - 1],
                     data_stream);
      } else {
        data_stream->AddCode(quant_msb, 0);
        data_stream->AddCode(quant_lsb, 1);
      }
    }
    for (size_t c = 0; c < num_channels; ++c) {
      for (size_t b = 0; b < kDctNumber; ++b) {
        const size_t offset = b * kDctLength;
        const int32_t* block = &ringli_blocks[i].channels[c][offset];
        data_stream->ResizeForBlock();

        int last_nz = 0;
        for (int k = 1; k < kDctLength; ++k) {
          if (block[k]) last_nz = k;
        }
        EncodeSymbol(last_nz, kDctLength, &last_nz_prob[0], data_stream);
        int num_nzeros = 0;
        for (int k = last_nz; k >= 0; --k) {
          const int coeff = block[k];
          const int is_zero = (coeff == 0);
          if (k == 0 || k < last_nz) {
            const int is_zero_ctx = ZeronessContext(num_nzeros, k);
            Prob* const p = &is_zero_prob[is_zero_ctx];
            data_stream->AddBit(p, is_zero);
          }
          if (!is_zero) {
            const int sign = (coeff > 0 ? 0 : 1);
            const size_t sign_ctx = k;
            Prob* const p = &sign_prob[sign_ctx];
            data_stream->AddBit(p, sign);
            const int absval = sign ? -coeff : coeff;
            const size_t absval_ctx = 2 + ZeroDensityContext(num_nzeros, k);
            int symbol;
            int nbits = 0;
            int extra_bits = 0;
            if (absval <= NUM_DIRECT_CODES) {
              symbol = absval - 1;
            } else {
              nbits = Log2FloorNonZero(absval - NUM_DIRECT_CODES + 1);
              symbol = NUM_DIRECT_CODES + nbits - 1;
              extra_bits = absval - (NUM_DIRECT_CODES - 1 + (1 << nbits));
              extra_bits &= (1 << nbits) - 1;
            }
            if (ecparams.arithmetic_only) {
              EncodeSymbol(symbol, MAX_SYMBOLS,
                           &symbol_prob[absval_ctx * (MAX_SYMBOLS - 1)],
                           data_stream);
            } else {
              data_stream->AddCode(symbol, absval_ctx);
            }
            if (nbits > 0) {
              data_stream->AddBits(nbits, extra_bits);
            }
            ++num_nzeros;
          }
        }
      }
    }
  }
  if (!ecparams.arithmetic_only) {
    entropy_source->ClusterHistograms();
  }
}

bool CompressCoefficients(const std::vector<RingliBlock>& ringli_blocks,
                          size_t num_channels,
                          const EntropyCodingParams& ecparams,
                          std::string* output) {
  EntropySource entropy_source;
  DataStream data_stream(&entropy_source);
  ProcessCoefficients(ringli_blocks, num_channels, ecparams, &entropy_source,
                      &data_stream);
  const size_t num_coeffs =
      ringli_blocks.size() * num_channels * kRingliBlockSize;

  // TODO(szabadka) Improve on this bound.
  const size_t max_compressed_size = num_coeffs * 4 + (1 << 12);
  const size_t size_bytes = Base128Size(max_compressed_size);
  const size_t data_start = output->size();
  output->resize(data_start + max_compressed_size);
  size_t pos = data_start;

  if (!ecparams.arithmetic_only) {
    // Skip some bytes for the size of the histogram data.
    pos += size_bytes;
    uint8_t* storage = reinterpret_cast<uint8_t*>(&(*output)[pos]);
    size_t storage_ix = 0;
    WriteBitsPrepareStorage(storage_ix, storage);
    entropy_source.EncodeContextMap(&storage_ix, storage);
    entropy_source.BuildAndStoreEntropyCodes(&storage_ix, storage);
    // Encode histogram data size before the histogram data.
    const size_t histograms_size = (storage_ix + 7) >> 3;
    EncodeBase128Fix(histograms_size, size_bytes,
                     reinterpret_cast<uint8_t*>(&(*output)[data_start]));
    pos += histograms_size;
  }

  const size_t coeff_data_start = pos;
  // Skip some bytes for the size of the coefficient data.
  pos += size_bytes;
  data_stream.EncodeCodeWords(entropy_source, ecparams,
                              reinterpret_cast<uint8_t*>(&(*output)[0]), &pos,
                              output->size());
  // Encode coeff data size before the coeff data.
  const size_t coeff_data_size = pos - coeff_data_start - size_bytes;
  EncodeBase128Fix(coeff_data_size, size_bytes,
                   reinterpret_cast<uint8_t*>(&(*output)[coeff_data_start]));

  output->resize(pos);
  return true;
}

void AppendUint16ToString(void* opaque, uint16_t val) {
  std::string* s = reinterpret_cast<std::string*>(opaque);
  s->push_back(val & 0xff);
  s->push_back(val >> 8);
}

void WriteSymbol(int val, int alphabet_size, Prob* probs,
                 BinaryArithmeticEncoder* ac, std::string* output) {
  int val0 = 0;
  int val1 = alphabet_size;
  while (val0 + 1 < val1) {
    const int mid = (val0 + val1) >> 1;
    const int bit = (val >= mid);
    Prob& p = probs[mid - 1];
    const uint8_t prob = p.get_proba();
    p.Add(bit);
    ac->AddBit(prob, bit, output, AppendUint16ToString);
    if (bit) {
      val0 = mid;
    } else {
      val1 = mid;
    }
  }
}

EntropyCoder::EntropyCoder(const EntropyCodingParams& ecparams,
                           uint32_t sampling_freq, uint32_t num_channels,
                           bool predictive, bool online)
    : ecparams_(ecparams),
      sampling_freq_(sampling_freq),
      num_channels_(num_channels),
      predictive_(predictive),
      online_(online) {
  Reset();
}

void EntropyCoder::Reset() {
  if (predictive_) {
    if (online_ && ecparams_.arithmetic_only) {
      context_model_.resize(num_channels_);
      symbol_prob_.resize(context_model_[0].NumContexts() * (MAX_SYMBOLS - 1));
    } else {
      context_model_.resize(1);
      entropy_source_ = std::make_unique<EntropySource>();
      data_stream_ = std::make_unique<DataStream>(entropy_source_.get());
      const size_t num_contexts =
          3 + kNumLSFContexts + context_model_[0].NumContexts();
      entropy_source_->Resize(num_contexts);
      if (ecparams_.arithmetic_only) {
        symbol_prob_.resize(num_contexts * (MAX_SYMBOLS - 1));
      }
      memset(order_histo_, 0, sizeof(order_histo_));
      lsf_extra_bits_ = 0;
    }
  }
  num_samples_ = 0;
  arith_encode_.Reset();
  idx_ = 0;
}

void EntropyCoder::ProcessSamples(const int* samples, std::string* output) {
  for (uint32_t ci = 0; ci < num_channels_; ++ci) {
    const int val = samples[ci];
    int nbits = 0;
    int extra_bits = 0;
    const int symbol =
        EncodeValue(val, kPredNumDirectAbsval, &nbits, &extra_bits);
    const int ctx = context_model_[ci].Context();
    WriteSymbol(symbol, MAX_SYMBOLS, &symbol_prob_[ctx * (MAX_SYMBOLS - 1)],
                &arith_encode_, output);
    for (int b = 0; b < nbits; ++b) {
      arith_encode_.AddBit(128, (extra_bits >> b) & 1, output,
                           AppendUint16ToString);
    }
    context_model_[ci].Add(val);
  }
  ++idx_;
  if (idx_ == kRingliBlockSize) {
    for (uint32_t ci = 0; ci < num_channels_; ++ci) {
      context_model_[ci].Reset();
    }
    idx_ = 0;
  }
}

bool EntropyCoder::ProcessPredictiveBlock(const RingliBlock& block,
                                          std::string* output) {
  for (uint32_t ci = 0; ci < num_channels_; ++ci) {
    data_stream_->ResizeForBlock();
    const auto& header = block.header.pred[ci];
    if (!online_) {
      const int order = header.quant_lsf.size();
      if (ecparams_.arithmetic_only) {
        WriteSymbol(order, MAX_SYMBOLS, &symbol_prob_[2 * (MAX_SYMBOLS - 1)],
                    &arith_encode_, output);
      } else {
        data_stream_->AddCode(order, 2);
      }
      ++order_histo_[order];
      for (int p = 0; p < order; ++p) {
        const int pred_lsf = p * (kLSFQuant[p] / order);
        const int residual = header.quant_lsf[p] - pred_lsf;
        int nbits = 0;
        int extra_bits = 0;
        const int symbol = EncodeValue(residual, 16, &nbits, &extra_bits);
        const int ctx = 3 + LSFContext(p, order);
        if (ecparams_.arithmetic_only) {
          WriteSymbol(symbol, MAX_SYMBOLS,
                      &symbol_prob_[ctx * (MAX_SYMBOLS - 1)], &arith_encode_,
                      output);
          for (int b = 0; b < nbits; ++b) {
            arith_encode_.AddBit(128, (extra_bits >> b) & 1, output,
                                 AppendUint16ToString);
          }
        } else {
          data_stream_->AddCode(symbol, ctx);
          if (nbits > 0) {
            data_stream_->AddBits(nbits, extra_bits);
          }
        }
        lsf_extra_bits_ += nbits;
      }
    }
    const auto& channel = block.channels[ci];
    context_model_[0].Reset();
    for (int i = 0; i < kRingliBlockSize; ++i) {
      const int val = channel[i];
      int nbits = 0;
      int extra_bits = 0;
      const int symbol =
          EncodeValue(val, kPredNumDirectAbsval, &nbits, &extra_bits);
      const int ctx = 3 + kNumLSFContexts + context_model_[0].Context();
      if (ecparams_.arithmetic_only) {
        WriteSymbol(symbol, MAX_SYMBOLS, &symbol_prob_[ctx * (MAX_SYMBOLS - 1)],
                    &arith_encode_, output);
        for (int b = 0; b < nbits; ++b) {
          arith_encode_.AddBit(128, (extra_bits >> b) & 1, output,
                               AppendUint16ToString);
        }
      } else {
        data_stream_->AddCode(symbol, ctx);
        if (nbits > 0) {
          data_stream_->AddBits(nbits, extra_bits);
        }
      }
      context_model_[0].Add(val);
    }
  }
  return true;
}

bool EntropyCoder::ProcessBlock(const RingliBlock& block, std::string* output) {
  num_samples_ += kRingliBlockSize * num_channels_;
  if (predictive_) {
    return ProcessPredictiveBlock(block, output);
  }
  return false;
}

bool EntropyCoder::Flush(std::string* output) {
  if (ecparams_.arithmetic_only) {
    arith_encode_.Flush(output, AppendUint16ToString);
    return true;
  }
  const double duration = 1.0 * num_samples_ / num_channels_ / sampling_freq_;
  if (!ecparams_.arithmetic_only) {
    entropy_source_->ClusterHistograms();
  }
  PrintHistogram("predictor order", order_histo_, kMaxPredictorOrder + 1);

  // TODO(szabadka) Improve on this bound.
  const size_t max_compressed_size =
      (num_samples_ * sizeof(int32_t)) + (1 << 12);
  const size_t size_bytes = Base128Size(max_compressed_size);
  const size_t data_start = output->size();
  output->resize(data_start + max_compressed_size);
  size_t pos = data_start;

  if (!ecparams_.arithmetic_only) {
    // Skip some bytes for the size of the histogram data.
    pos += size_bytes;
    uint8_t* storage = reinterpret_cast<uint8_t*>(&(*output)[pos]);
    size_t storage_ix = 0;
    WriteBitsPrepareStorage(storage_ix, storage);
    entropy_source_->EncodeContextMap(&storage_ix, storage);
    entropy_source_->BuildAndStoreEntropyCodes(&storage_ix, storage);
    // Encode histogram data size before the histogram data.
    const size_t histograms_size = (storage_ix + 7) >> 3;
    EncodeBase128Fix(histograms_size, size_bytes,
                     reinterpret_cast<uint8_t*>(&(*output)[data_start]));
    pos += histograms_size;
    PrintSize("histograms", histograms_size, duration);
  }

  const size_t compressed_data_start = pos;
  data_stream_->EncodeCodeWords(*entropy_source_, ecparams_,
                                reinterpret_cast<uint8_t*>(&(*output)[0]), &pos,
                                output->size());
  const int total_extra_bits = data_stream_->TotalExtraBits();
  PrintSize("extra bits", total_extra_bits / 8, duration);
  PrintSize("entropy coded", pos - compressed_data_start - total_extra_bits / 8,
            duration);
  if (absl::GetFlag(FLAGS_log_level) >= 1) {
    printf("Num histograms: %zu\n", entropy_source_->NumHistograms());
  }
  if (absl::GetFlag(FLAGS_log_level) >= 2 && !ecparams_.arithmetic_only) {
    double quantizer_entropy = entropy_source_->ClusteredEntropy(0) +
                               entropy_source_->ClusteredEntropy(1);
    PrintSize("Quantizer entropy", quantizer_entropy / 8, duration);
    double pred_order_entropy = entropy_source_->ClusteredEntropy(2);
    PrintSize("Pred order entropy", pred_order_entropy / 8, duration);
    double lsf_entropy = 0.0;
    for (int ctx = 0; ctx < kNumLSFContexts; ++ctx) {
      lsf_entropy += entropy_source_->ClusteredEntropy(3 + ctx);
    }
    PrintSize("LSF entropy", (lsf_extra_bits_ + lsf_entropy) / 8, duration);
    double residual_entropy = 0.0;
    for (int ctx = 0; ctx < context_model_[0].NumContexts(); ++ctx) {
      residual_entropy +=
          entropy_source_->ClusteredEntropy(3 + kNumLSFContexts + ctx);
    }
    PrintSize("Residual entropy", residual_entropy / 8, duration);
  }
  output->resize(pos);
  return true;
}

}  // namespace ringli
