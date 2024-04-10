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

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "absl/log/check.h"
#include "common/wav_header.h"
#include "common/wav_reader.h"
#include "common/wav_writer.h"
#include "opus.h"
#include "opus_defines.h"
#include "opus_multistream.h"

namespace ringli {
namespace {

uint32_t UpdateCRC(uint32_t crc, const void* p, size_t len) {
  // TODO(szabadka) This is slow. Figure out how to use util/hash or some
  // other library to produce this result.
  static const uint32_t kCRCTable[] = {
      0x00000000, 0xB71DC104, 0x6E3B8209, 0xD926430D, 0xDC760413, 0x6B6BC517,
      0xB24D861A, 0x0550471E, 0xB8ED0826, 0x0FF0C922, 0xD6D68A2F, 0x61CB4B2B,
      0x649B0C35, 0xD386CD31, 0x0AA08E3C, 0xBDBD4F38, 0x70DB114C, 0xC7C6D048,
      0x1EE09345, 0xA9FD5241, 0xACAD155F, 0x1BB0D45B, 0xC2969756, 0x758B5652,
      0xC836196A, 0x7F2BD86E, 0xA60D9B63, 0x11105A67, 0x14401D79, 0xA35DDC7D,
      0x7A7B9F70, 0xCD665E74, 0xE0B62398, 0x57ABE29C, 0x8E8DA191, 0x39906095,
      0x3CC0278B, 0x8BDDE68F, 0x52FBA582, 0xE5E66486, 0x585B2BBE, 0xEF46EABA,
      0x3660A9B7, 0x817D68B3, 0x842D2FAD, 0x3330EEA9, 0xEA16ADA4, 0x5D0B6CA0,
      0x906D32D4, 0x2770F3D0, 0xFE56B0DD, 0x494B71D9, 0x4C1B36C7, 0xFB06F7C3,
      0x2220B4CE, 0x953D75CA, 0x28803AF2, 0x9F9DFBF6, 0x46BBB8FB, 0xF1A679FF,
      0xF4F63EE1, 0x43EBFFE5, 0x9ACDBCE8, 0x2DD07DEC, 0x77708634, 0xC06D4730,
      0x194B043D, 0xAE56C539, 0xAB068227, 0x1C1B4323, 0xC53D002E, 0x7220C12A,
      0xCF9D8E12, 0x78804F16, 0xA1A60C1B, 0x16BBCD1F, 0x13EB8A01, 0xA4F64B05,
      0x7DD00808, 0xCACDC90C, 0x07AB9778, 0xB0B6567C, 0x69901571, 0xDE8DD475,
      0xDBDD936B, 0x6CC0526F, 0xB5E61162, 0x02FBD066, 0xBF469F5E, 0x085B5E5A,
      0xD17D1D57, 0x6660DC53, 0x63309B4D, 0xD42D5A49, 0x0D0B1944, 0xBA16D840,
      0x97C6A5AC, 0x20DB64A8, 0xF9FD27A5, 0x4EE0E6A1, 0x4BB0A1BF, 0xFCAD60BB,
      0x258B23B6, 0x9296E2B2, 0x2F2BAD8A, 0x98366C8E, 0x41102F83, 0xF60DEE87,
      0xF35DA999, 0x4440689D, 0x9D662B90, 0x2A7BEA94, 0xE71DB4E0, 0x500075E4,
      0x892636E9, 0x3E3BF7ED, 0x3B6BB0F3, 0x8C7671F7, 0x555032FA, 0xE24DF3FE,
      0x5FF0BCC6, 0xE8ED7DC2, 0x31CB3ECF, 0x86D6FFCB, 0x8386B8D5, 0x349B79D1,
      0xEDBD3ADC, 0x5AA0FBD8, 0xEEE00C69, 0x59FDCD6D, 0x80DB8E60, 0x37C64F64,
      0x3296087A, 0x858BC97E, 0x5CAD8A73, 0xEBB04B77, 0x560D044F, 0xE110C54B,
      0x38368646, 0x8F2B4742, 0x8A7B005C, 0x3D66C158, 0xE4408255, 0x535D4351,
      0x9E3B1D25, 0x2926DC21, 0xF0009F2C, 0x471D5E28, 0x424D1936, 0xF550D832,
      0x2C769B3F, 0x9B6B5A3B, 0x26D61503, 0x91CBD407, 0x48ED970A, 0xFFF0560E,
      0xFAA01110, 0x4DBDD014, 0x949B9319, 0x2386521D, 0x0E562FF1, 0xB94BEEF5,
      0x606DADF8, 0xD7706CFC, 0xD2202BE2, 0x653DEAE6, 0xBC1BA9EB, 0x0B0668EF,
      0xB6BB27D7, 0x01A6E6D3, 0xD880A5DE, 0x6F9D64DA, 0x6ACD23C4, 0xDDD0E2C0,
      0x04F6A1CD, 0xB3EB60C9, 0x7E8D3EBD, 0xC990FFB9, 0x10B6BCB4, 0xA7AB7DB0,
      0xA2FB3AAE, 0x15E6FBAA, 0xCCC0B8A7, 0x7BDD79A3, 0xC660369B, 0x717DF79F,
      0xA85BB492, 0x1F467596, 0x1A163288, 0xAD0BF38C, 0x742DB081, 0xC3307185,
      0x99908A5D, 0x2E8D4B59, 0xF7AB0854, 0x40B6C950, 0x45E68E4E, 0xF2FB4F4A,
      0x2BDD0C47, 0x9CC0CD43, 0x217D827B, 0x9660437F, 0x4F460072, 0xF85BC176,
      0xFD0B8668, 0x4A16476C, 0x93300461, 0x242DC565, 0xE94B9B11, 0x5E565A15,
      0x87701918, 0x306DD81C, 0x353D9F02, 0x82205E06, 0x5B061D0B, 0xEC1BDC0F,
      0x51A69337, 0xE6BB5233, 0x3F9D113E, 0x8880D03A, 0x8DD09724, 0x3ACD5620,
      0xE3EB152D, 0x54F6D429, 0x7926A9C5, 0xCE3B68C1, 0x171D2BCC, 0xA000EAC8,
      0xA550ADD6, 0x124D6CD2, 0xCB6B2FDF, 0x7C76EEDB, 0xC1CBA1E3, 0x76D660E7,
      0xAFF023EA, 0x18EDE2EE, 0x1DBDA5F0, 0xAAA064F4, 0x738627F9, 0xC49BE6FD,
      0x09FDB889, 0xBEE0798D, 0x67C63A80, 0xD0DBFB84, 0xD58BBC9A, 0x62967D9E,
      0xBBB03E93, 0x0CADFF97, 0xB110B0AF, 0x060D71AB, 0xDF2B32A6, 0x6836F3A2,
      0x6D66B4BC, 0xDA7B75B8, 0x035D36B5, 0xB440F7B1,
  };
  const uint8_t* data = reinterpret_cast<const uint8_t*>(p);
  for (int i = 0; i < len; ++i) {
    crc = kCRCTable[(crc & 0xff) ^ data[i]] ^ (crc >> 8);
  }
  return crc;
}

uint32_t BSwap32(uint32_t x) {
  return ((((x >> 24) & 0xff) << 0) | (((x >> 16) & 0xff) << 8) |
          (((x >> 8) & 0xff) << 16) | (((x >> 0) & 0xff) << 24));
}

struct OggPageHeader {
  uint8_t version;
  uint8_t header_type;
  uint64_t granule_position;
  uint32_t serial_number;
  uint32_t sequence_number;
  uint32_t crc;
  std::vector<uint16_t> packet_lengths;
  uint16_t data_length;
};

struct OpusHeader {
  uint8_t version;
  uint8_t num_channels;
  uint16_t pre_skip;
  uint32_t input_sampling_frequency;
  uint16_t output_gain;
  uint8_t mapping_family;
  std::vector<uint8_t> mapping_table;
};

struct OpusTags {
  std::string vendor;
  std::vector<std::string> comments;
};

template <typename T>
void Write(const T& t, std::string* output) {
  const size_t pos = output->size();
  output->resize(pos + sizeof(T));
  memcpy(&(*output)[pos], &t, sizeof(T));
}

// Writes the given data packets to an ogg page as specified
// in https://datatracker.ietf.org/doc/html/rfc3533
void WriteOggPage(const std::vector<uint16_t> packet_lengths,
                  uint32_t sequence_number, uint64_t granule_position,
                  uint32_t serial_number, uint8_t is_last, const uint8_t* data,
                  std::string* output) {
  const size_t start_pos = output->size();
  output->append("OggS", 4);
  const uint8_t version = 0;
  Write(version, output);
  const uint8_t is_first = sequence_number == 0;
  const uint8_t header_type = (is_first << 1) | (is_last << 2);
  Write(header_type, output);
  Write(granule_position, output);
  Write(serial_number, output);
  Write(sequence_number, output);
  const size_t crc_pos = output->size();
  uint32_t crc32 = 0;  // Will be filled in later.
  Write(crc32, output);
  uint8_t num_segments = 0;
  uint16_t total_len = 0;
  for (uint16_t packet_len : packet_lengths) {
    int new_segments = 1 + packet_len / 255;
    CHECK_LE(new_segments + num_segments, 255);
    num_segments += new_segments;
    total_len += packet_len;
  }
  Write(num_segments, output);
  for (uint16_t packet_len : packet_lengths) {
    uint8_t segment_len;
    do {
      segment_len = std::min<uint16_t>(255, packet_len);
      Write(segment_len, output);
      packet_len -= segment_len;
    } while (segment_len == 255);
  }
  const size_t pos = output->size();
  output->resize(pos + total_len);
  memcpy(&(*output)[pos], reinterpret_cast<const char*>(data), total_len);
  // Fill in CRC checksum.
  crc32 = UpdateCRC(crc32, &(*output)[start_pos], output->size() - start_pos);
  crc32 = BSwap32(crc32);
  memcpy(&(*output)[crc_pos], &crc32, 4);
}

void FillOpusHeader(const FormatChunk& format, uint16_t pre_skip,
                    OpusHeader* opus_header) {
  opus_header->version = 1;
  opus_header->num_channels = format.number_of_channels;
  opus_header->pre_skip = pre_skip;
  opus_header->input_sampling_frequency = format.sampling_frequency;
  opus_header->output_gain = 0;
  opus_header->mapping_family = 0;
}

void WriteOpusHeader(const OpusHeader& opus_header, uint32_t serial_number,
                     uint32_t* sequence_number, std::string* output) {
  std::string opus_packet;
  opus_packet.append("OpusHead", 8);
  Write(opus_header.version, &opus_packet);
  Write(opus_header.num_channels, &opus_packet);
  Write(opus_header.pre_skip, &opus_packet);
  Write(opus_header.input_sampling_frequency, &opus_packet);
  Write(opus_header.output_gain, &opus_packet);
  Write(opus_header.mapping_family, &opus_packet);
  if (opus_header.mapping_family != 0) {
    for (uint8_t val : opus_header.mapping_table) {
      Write(val, &opus_packet);
    }
  }
  WriteOggPage(std::vector<uint16_t>(1, opus_packet.size()),
               (*sequence_number)++, /*granule_position=*/0, serial_number,
               /*is_last=*/false,
               reinterpret_cast<const uint8_t*>(opus_packet.data()), output);
  // Write placeholder OpusTags packet
  opus_packet.clear();
  opus_packet.append("OpusTags", 8);
  const std::string kVendor("Google");
  uint32_t len = kVendor.size();
  Write(len, &opus_packet);
  opus_packet.append(kVendor);
  len = 0;
  Write(len, &opus_packet);
  WriteOggPage(std::vector<uint16_t>(1, opus_packet.size()),
               (*sequence_number)++, /*granule_position=*/0, serial_number,
               /*is_last=*/false,
               reinterpret_cast<const uint8_t*>(opus_packet.data()), output);
}

template <typename T>
bool Read(const uint8_t* data, size_t len, size_t* pos, T* t) {
  if (*pos + sizeof(T) > len) return false;
  memcpy(t, data + (*pos), sizeof(T));
  *pos += sizeof(T);
  return true;
}

bool Read(const uint8_t* data, size_t len, size_t* pos, char* id) {
  if (*pos + 4 > len) return false;
  memcpy(id, data + (*pos), 4);
  *pos += 4;
  return true;
}

bool Read(const uint8_t* data, size_t len, size_t* pos, std::string* out) {
  uint32_t length;
  if (!Read(data, len, pos, &length)) {
    return false;
  }
  out->resize(length);
  if (*pos + length > len) return false;
  memcpy(reinterpret_cast<uint8_t*>(out->data()), data + (*pos), length);
  *pos += length;
  return true;
}

bool ParseOggHeader(const uint8_t* data, size_t len, size_t* pos,
                    uint32_t* sequence_number, OggPageHeader* header) {
  const size_t start_pos = *pos;
  char magic[4];
  uint8_t num_segments;
  if (!Read(data, len, pos, magic) || (memcmp(magic, "OggS", 4) != 0) ||
      !Read(data, len, pos, &header->version) ||
      !Read(data, len, pos, &header->header_type) ||
      !Read(data, len, pos, &header->granule_position) ||
      !Read(data, len, pos, &header->serial_number) ||
      !Read(data, len, pos, &header->sequence_number)) {
    return false;
  }
  const size_t crc_pos = *pos;
  if (!Read(data, len, pos, &header->crc) ||
      !Read(data, len, pos, &num_segments)) {
    return false;
  }
  header->crc = BSwap32(header->crc);
  if (header->version != 0) {
    fprintf(stderr, "Invalid page version %d\n", header->version);
    return false;
  }
  if (header->header_type & 0x01) {
    fprintf(stderr, "Continuation pages are not supported\n");
    return false;
  }
  const bool bos = (header->sequence_number == 0);
  if (!!(header->header_type & 0x02) != bos) {
    fprintf(stderr, "Invalid bos page header flag\n");
    return false;
  }
  if (header->sequence_number != *sequence_number) {
    fprintf(stderr, "Out of oder page sequence number %d, expected %d\n",
            header->sequence_number, *sequence_number);
    return false;
  }
  ++(*sequence_number);
  header->data_length = 0;
  uint16_t packet_length = 0;
  for (size_t i = 0; i < num_segments; ++i) {
    uint8_t segment_len;
    if (!Read(data, len, pos, &segment_len)) {
      return false;
    }
    packet_length += segment_len;
    if (segment_len < 255) {
      header->packet_lengths.push_back(packet_length);
      header->data_length += packet_length;
      packet_length = 0;
    } else if (i + 1 == num_segments) {
      fprintf(stderr, "Data packet spanning multiple pages is not supported\n");
      return false;
    }
  }
  const size_t end_pos = *pos + header->data_length;
  if (end_pos > len) {
    fprintf(stderr, "Invalid page data length %d (remaining: %zu)\n",
            header->data_length, len - *pos);
    return false;
  }
  const bool eos = (end_pos == len);
  if (!!(header->header_type & 0x04) != eos) {
    fprintf(stderr, "Invalid eos page header flag\n");
    return false;
  }
  // Verify the CRC checksum;
  uint32_t crc32 = 0;
  crc32 = UpdateCRC(crc32, &data[start_pos], crc_pos - start_pos);
  uint32_t zero = 0;
  crc32 = UpdateCRC(crc32, &zero, 4);
  crc32 = UpdateCRC(crc32, &data[crc_pos + 4], end_pos - crc_pos - 4);
  if (crc32 != header->crc) {
    fprintf(stderr, "CRC checksum failure: 0x%08x expected: 0x%08x\n", crc32,
            header->crc);
    return false;
  }
  return true;
}

bool ParseOpusHeader(const uint8_t* data, size_t len, size_t* pos,
                     OpusHeader* header) {
  char magic[8];
  if (!Read(data, len, pos, magic) || !Read(data, len, pos, magic + 4) ||
      (memcmp(magic, "OpusHead", 8) != 0) ||
      !Read(data, len, pos, &header->version) ||
      !Read(data, len, pos, &header->num_channels) ||
      !Read(data, len, pos, &header->pre_skip) ||
      !Read(data, len, pos, &header->input_sampling_frequency) ||
      !Read(data, len, pos, &header->output_gain) ||
      !Read(data, len, pos, &header->mapping_family)) {
    return false;
  }
  if (header->version != 1) {
    fprintf(stderr, "Invalid opus version %d\n", header->version);
    return false;
  }
  if (header->mapping_family == 0) {
    return true;
  }
  header->mapping_table.resize(header->num_channels);
  for (size_t c = 0; c < header->num_channels; ++c) {
    if (!Read(data, len, pos, &header->mapping_table[c])) {
      return false;
    }
  }
  return true;
}

bool ParseOpusTags(const uint8_t* data, size_t len, size_t* pos,
                   OpusTags* tags) {
  char magic[8];
  if (!Read(data, len, pos, magic) || !Read(data, len, pos, magic + 4) ||
      (memcmp(magic, "OpusTags", 8) != 0) ||
      !Read(data, len, pos, &tags->vendor)) {
    return false;
  }
  uint32_t num_comments;
  if (!Read(data, len, pos, &num_comments)) {
    return false;
  }
  tags->comments.resize(num_comments);
  for (uint32_t i = 0; i < num_comments; ++i) {
    if (!Read(data, len, pos, &tags->comments[i])) {
      return false;
    }
  }
  return true;
}

void FillWavHeader(const OpusHeader& opus_header, size_t wav_len,
                   WavHeader* wav_header) {
  memcpy(wav_header->riff_chunk.riff_chunk_id, "RIFF", 4);
  wav_header->riff_chunk.riff_chunk_size = wav_len - 8;
  memcpy(wav_header->riff_chunk.wave_format, "WAVE", 4);
  memcpy(wav_header->format_chunk.format_chunk_id, "fmt ", 4);
  wav_header->format_chunk.format_chunk_size = 16;
  wav_header->format_chunk.audio_format = 1;
  wav_header->format_chunk.number_of_channels = opus_header.num_channels;
  // TODO(szabadka): Support decoding to input frequency;
  wav_header->format_chunk.sampling_frequency = 48000;
  wav_header->format_chunk.bits_per_sample = 16;
  wav_header->format_chunk.block_align =
      wav_header->format_chunk.number_of_channels *
      (wav_header->format_chunk.bits_per_sample / 8);
  wav_header->format_chunk.byte_rate =
      wav_header->format_chunk.sampling_frequency *
      wav_header->format_chunk.number_of_channels *
      (wav_header->format_chunk.bits_per_sample / 8);
  memcpy(wav_header->channel_id, "data", 4);
  wav_header->channel_data_length = wav_len - kWavHeaderSize;
}

bool DecodeOpusPackets(const OpusHeader& opus_header,
                       const OggPageHeader& page_header,
                       const uint8_t* page_data, uint64_t* total_samples,
                       OpusDecoder* dec, std::string* wav_data) {
  constexpr int kMaxFrameSize = 96000;
  std::vector<int16_t> output(kMaxFrameSize * opus_header.num_channels);
  if (page_header.granule_position < *total_samples) {
    fprintf(stderr, "Invalid granule position %lu, expected at least %lu\n",
            page_header.granule_position, *total_samples);
    return false;
  }
  const int num_channels = opus_header.num_channels;
  const uint64_t input_samples = page_header.granule_position - *total_samples;
  wav_data->reserve(wav_data->size() +
                    input_samples * num_channels * sizeof(int16_t));
  size_t pos = 0;
  for (size_t i = 0; i < page_header.packet_lengths.size(); ++i) {
    const size_t packet_len = page_header.packet_lengths[i];
    int output_samples = opus_decode(dec, &page_data[pos], packet_len,
                                     output.data(), kMaxFrameSize, 0);
    if (output_samples < 0) {
      fprintf(stderr, "Error decoding opus packet: %s\n",
              opus_strerror(output_samples));
      return false;
    }
    if (*total_samples + output_samples >= opus_header.pre_skip) {
      const uint16_t skip = *total_samples < opus_header.pre_skip
                                ? opus_header.pre_skip - *total_samples
                                : 0;
      output_samples = std::min<uint64_t>(
          output_samples, page_header.granule_position - *total_samples);
      WriteSamples(output.data() + skip * num_channels,
                   (output_samples - skip) * num_channels, wav_data);
    }
    *total_samples += output_samples;
    pos += packet_len;
  }
  if (*total_samples != page_header.granule_position) {
    fprintf(stderr, "Unexpected total samples: %lu vs %lu\n", *total_samples,
            page_header.granule_position);
    return false;
  }
  return true;
}

bool OpusDecompress(const std::string& opus_data, std::string* wav_data) {
  // Make room for the wav header, which will be generated later based on
  // the total number of output_samples.
  wav_data->append(kWavHeaderSize, 0);
  const uint8_t* data = reinterpret_cast<const uint8_t*>(opus_data.data());
  size_t pos = 0;
  uint32_t sequence_number = 0;
  OpusHeader opus_header;
  std::unique_ptr<OpusDecoder, void (*)(OpusDecoder*)> dec(
      nullptr, opus_decoder_destroy);
  size_t total_samples = 0;
  while (pos < opus_data.size()) {
    OggPageHeader page_header;
    if (!ParseOggHeader(data, opus_data.size(), &pos, &sequence_number,
                        &page_header)) {
      fprintf(stderr, "Failed to parse ogg page header\n");
      return false;
    }
    const size_t page_end = pos + page_header.data_length;
    if (page_header.sequence_number == 0) {
      if (!ParseOpusHeader(data, page_end, &pos, &opus_header)) {
        fprintf(stderr, "Failed to parse Opus header\n");
        return false;
      }
      int error;
      dec.reset(opus_decoder_create(48000, opus_header.num_channels, &error));
      if (error != OPUS_OK) {
        fprintf(stderr, "Opus error: %s\n", opus_strerror(error));
        return false;
      }
    } else if (page_header.sequence_number == 1) {
      OpusTags opus_tags;
      if (!ParseOpusTags(data, page_end, &pos, &opus_tags)) {
        fprintf(stderr, "Failed to parse Opus tags\n");
        return false;
      }
    } else {
      if (!DecodeOpusPackets(opus_header, page_header, &data[pos],
                             &total_samples, dec.get(), wav_data)) {
        return false;
      }
    }
    pos = page_end;
  }
  // Fill in the wav header at the start of the stream with the actual data
  // sizes.
  WavHeader wav_header;
  FillWavHeader(opus_header, wav_data->size(), &wav_header);
  std::string wav_header_str;
  WriteWavHeader(wav_header, &wav_header_str);
  CHECK_EQ(wav_header_str.size(), kWavHeaderSize);
  memcpy(wav_data->data(), wav_header_str.data(), kWavHeaderSize);
  return true;
}

}  // namespace

StreamingOpusEncoder::StreamingOpusEncoder(int bitrate)
    : bitrate_(bitrate), enc_(nullptr, opus_multistream_encoder_destroy) {
  format_.format_chunk_size = 0;
  wav_reader_.RegisterCallback("fmt ", this, ParseFormatCb,
                               sizeof(FormatChunk));
}

void StreamingOpusEncoder::Reset() {
  opus_data_.clear();
  wav_reader_.Reset();
  format_.format_chunk_size = 0;
  output_pos_ = 0;
  sequence_number_ = 0;
  total_encoded_ = 0;
  page_pos_ = 0;
  num_segments_ = 0;
  packet_lengths_.clear();
}

bool StreamingOpusEncoder::ParseFormatCb(void* opaque, const uint8_t* data,
                                         size_t len, size_t chunk_pos,
                                         size_t chunk_size) {
  auto enc = reinterpret_cast<StreamingOpusEncoder*>(opaque);
  if (chunk_size != len) {
    fprintf(stderr, "Invalid format chunk size %zu\n", chunk_size);
    return false;
  }
  enc->format_.format_chunk_size = chunk_size;
  size_t pos = 0;
  if (!ParseFormatChunk(data, len, &pos, &enc->format_)) {
    return false;
  }
  return enc->InitForFormat();
}

bool StreamingOpusEncoder::InitForFormat() {
  const size_t num_channels = format_.number_of_channels;
  const size_t bytes_per_sample = format_.bits_per_sample / 8;
  const size_t fs = format_.sampling_frequency;
  framesize_ = fs / 50;
  const size_t block_size = framesize_ * num_channels * bytes_per_sample;
  wav_reader_.RegisterCallback("data", this, ProcessDataCb, block_size);

  if (fs != 48000) {
    fprintf(stderr, "Only 48 kHz sampling rate is supported.");
    return false;
  }
  int coupled_streams = num_channels / 2;
  int streams = coupled_streams + (num_channels % 2);
  std::vector<uint8_t> channel_map(num_channels);
  for (int i = 0; i < num_channels; ++i) {
    channel_map[i] = i;
  }
  int error;
  enc_.reset(opus_multistream_encoder_create(
      fs, num_channels, streams, coupled_streams, channel_map.data(),
      OPUS_APPLICATION_AUDIO, &error));
  if (error != OPUS_OK) {
    fprintf(stderr, "Opus error: %s\n", opus_strerror(error));
    return false;
  }
  opus_multistream_encoder_ctl(enc_.get(), OPUS_SET_BITRATE(bitrate_ * 1000));
  opus_multistream_encoder_ctl(enc_.get(), OPUS_SET_COMPLEXITY(10));
  opus_multistream_encoder_ctl(enc_.get(), OPUS_SET_VBR(1));
  opus_multistream_encoder_ctl(enc_.get(), OPUS_SET_VBR_CONSTRAINT(0));

  opus_multistream_encoder_ctl(enc_.get(), OPUS_GET_LOOKAHEAD(&pre_skip_));
  if (pre_skip_ > 0xffff) {
    fprintf(stderr, "Pre-skip too high: %d\n", pre_skip_);
    return false;
  }
  block_in_.resize(framesize_ * num_channels);
  page_data_.resize(1 << 17);
  return true;
}

bool StreamingOpusEncoder::ProcessDataCb(void* opaque, const uint8_t* data,
                                         size_t len, size_t chunk_pos,
                                         size_t chunk_size) {
  auto enc = reinterpret_cast<StreamingOpusEncoder*>(opaque);
  if (enc->format_.format_chunk_size == 0) {
    fprintf(stderr, "Missing fmt chunk.\n");
    return false;
  }
  return enc->ProcessData(data, len, chunk_pos, chunk_size);
}

bool StreamingOpusEncoder::ProcessData(const uint8_t* data, size_t len,
                                       size_t chunk_pos, size_t chunk_size) {
  if (chunk_pos == 0) {
    WriteHeader(chunk_size);
  }
  CopyBlock(data, len);
  return ProcessBlock();
}

void StreamingOpusEncoder::WriteHeader(size_t chunk_size) {
  input_samples_ = chunk_size / sizeof(int16_t) / format_.number_of_channels;
  const float duration = input_samples_ * 1.0 / format_.sampling_frequency;
  const int estimated_output_size = duration * bitrate_ * 128;
  opus_data_.reserve(estimated_output_size);
  OpusHeader opus_header;
  FillOpusHeader(format_, pre_skip_, &opus_header);
  WriteOpusHeader(opus_header, serial_number_, &sequence_number_, &opus_data_);
}

void StreamingOpusEncoder::CopyBlock(const uint8_t* data, size_t len) {
  const size_t num_channels = format_.number_of_channels;
  size_t pos = 0;
  for (size_t i = 0; i < framesize_; ++i) {
    for (size_t c = 0; c < num_channels; ++c) {
      int16_t val = 0;
      if (pos + sizeof(val) <= len) {
        memcpy(&val, data + pos, sizeof(val));
        pos += sizeof(val);
      }
      block_in_[i * num_channels + c] = val;
    }
  }
}

bool StreamingOpusEncoder::ProcessBlock() {
  const int kMaxPacket = 1276;  // recommended by opus API doc
  int nb = opus_multistream_encode(enc_.get(), block_in_.data(), framesize_,
                                   &page_data_[page_pos_], kMaxPacket);
  if (nb < 0) {
    fprintf(stderr, "Error encoding opus packet: %s\n", opus_strerror(nb));
    return false;
  }
  if (nb > 0xffff) {
    fprintf(stderr, "Packet too long, nb = %d\n", nb);
    return false;
  }
  int encoded = opus_packet_get_samples_per_frame(&page_data_[page_pos_],
                                                  format_.sampling_frequency) *
                opus_packet_get_nb_frames(&page_data_[page_pos_], nb);
  if (encoded != framesize_) {
    fprintf(stderr,
            "Could not fit %zu sample frame into packet "
            "(encoded = %d)\n",
            framesize_, encoded);
  }
  int new_segments = 1 + nb / 255;
  if (num_segments_ + new_segments > 255) {
    WriteOggPage(packet_lengths_, sequence_number_++, total_encoded_,
                 serial_number_, /*is_last=*/false, page_data_.data(),
                 &opus_data_);
    for (int i = 0; i < nb; ++i) {
      page_data_[i] = page_data_[i + page_pos_];
    }
    page_pos_ = 0;
    num_segments_ = 0;
    packet_lengths_.clear();
  }
  page_pos_ += nb;
  num_segments_ += new_segments;
  packet_lengths_.push_back(nb);
  total_encoded_ += encoded;
  return true;
}

bool StreamingOpusEncoder::ProcessInput(const uint8_t* data, size_t len) {
  return wav_reader_.ProcessInput(data, len);
}

bool StreamingOpusEncoder::Flush() {
  while (total_encoded_ < input_samples_ + pre_skip_) {
    CopyBlock(nullptr, 0);
    if (!ProcessBlock()) {
      return false;
    }
  }
  total_encoded_ =
      std::min<uint64_t>(total_encoded_, input_samples_ + pre_skip_);
  WriteOggPage(packet_lengths_, sequence_number_++, total_encoded_,
               serial_number_, /*is_last=*/true, page_data_.data(),
               &opus_data_);
  return true;
}

size_t StreamingOpusEncoder::OutputSize() const {
  return opus_data_.size() - output_pos_;
}

size_t StreamingOpusEncoder::CopyOutput(uint8_t* buffer, size_t len) {
  size_t nbytes = std::min(len, OutputSize());
  memcpy(buffer, reinterpret_cast<const uint8_t*>(&opus_data_[output_pos_]),
         nbytes);
  output_pos_ += nbytes;
  return nbytes;
}

void StreamingOpusDecoder::Reset() {
  opus_data_.clear();
  wav_data_.clear();
  output_pos_ = 0;
}

bool StreamingOpusDecoder::ProcessInput(const uint8_t* data, size_t len) {
  opus_data_.append(reinterpret_cast<const char*>(data), len);
  return true;
}

bool StreamingOpusDecoder::Flush() {
  return OpusDecompress(opus_data_, &wav_data_);
}

size_t StreamingOpusDecoder::OutputSize() const {
  return wav_data_.size() - output_pos_;
}

size_t StreamingOpusDecoder::CopyOutput(uint8_t* buffer, size_t len) {
  size_t nbytes = std::min(len, OutputSize());
  memcpy(buffer, reinterpret_cast<const uint8_t*>(&wav_data_[output_pos_]),
         nbytes);
  output_pos_ += nbytes;
  return nbytes;
}

}  // namespace ringli
