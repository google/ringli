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

#include "analysis/flac_codec.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "FLAC/format.h"
#include "FLAC/ordinals.h"
#include "FLAC/stream_decoder.h"
#include "FLAC/stream_encoder.h"
#include "common/wav_header.h"
#include "common/wav_reader.h"
#include "common/wav_writer.h"

namespace ringli {

namespace {

FLAC__StreamEncoderWriteStatus AppendToString(
    const FLAC__StreamEncoder* encoder, const FLAC__byte buffer[], size_t bytes,
    uint32_t samples, uint32_t current_frame, void* client_data) {
  std::string* output = reinterpret_cast<std::string*>(client_data);
  output->append(reinterpret_cast<const char*>(&buffer[0]), bytes);
  return FLAC__STREAM_ENCODER_WRITE_STATUS_OK;
}

bool FlacCompress(const std::string& wav_data, int compression_level,
                  std::string* flac_data) {
  const uint8_t* data = reinterpret_cast<const uint8_t*>(wav_data.data());
  WavHeader wav_header;
  size_t pos = 0;
  if (!ParseWavHeader(data, wav_data.size(), &pos, &wav_header)) {
    return false;
  }
  const int num_channels = wav_header.format_chunk.number_of_channels;
  const int bits_per_sample = wav_header.format_chunk.bits_per_sample;
  const int fs = wav_header.format_chunk.sampling_frequency;
  std::unique_ptr<FLAC__StreamEncoder, void (*)(FLAC__StreamEncoder*)> enc(
      FLAC__stream_encoder_new(), FLAC__stream_encoder_delete);
  FLAC__stream_encoder_set_channels(enc.get(), num_channels);
  FLAC__stream_encoder_set_bits_per_sample(enc.get(), bits_per_sample);
  FLAC__stream_encoder_set_sample_rate(enc.get(), fs);
  FLAC__stream_encoder_set_compression_level(enc.get(), compression_level);
  FLAC__StreamEncoderInitStatus status = FLAC__stream_encoder_init_stream(
      enc.get(), AppendToString, nullptr, nullptr, nullptr, flac_data);
  if (status != FLAC__STREAM_ENCODER_INIT_STATUS_OK) {
    fprintf(stderr, "Failed to initialize flac stream: %s\n",
            FLAC__StreamEncoderInitStatusString[status]);
  }
  size_t remaining_samples =
      (wav_data.size() - pos) / sizeof(int16_t) / num_channels;
  static const size_t kBufferSize = 1024;
  std::vector<FLAC__int32> buffer(kBufferSize * num_channels);
  while (remaining_samples > 0) {
    size_t samples = std::min(remaining_samples, kBufferSize);
    for (size_t si = 0, i = 0; si < samples; ++si) {
      for (size_t ci = 0; ci < num_channels; ++ci, ++i) {
        int16_t value;
        memcpy(&value, data + pos, sizeof(value));
        pos += sizeof(value);
        buffer[i] = value;
      }
    }
    if (!FLAC__stream_encoder_process_interleaved(enc.get(), buffer.data(),
                                                  samples)) {
      fprintf(stderr, "Failed to encode samples: %s\n",
              FLAC__StreamEncoderStateString[FLAC__stream_encoder_get_state(
                  enc.get())]);
      return false;
    }
    remaining_samples -= samples;
  }
  if (!FLAC__stream_encoder_finish(enc.get())) {
    fprintf(stderr, "Failed to flush stream: %s\n",
            FLAC__StreamEncoderStateString[FLAC__stream_encoder_get_state(
                enc.get())]);
    return false;
  }
  return true;
}

struct FlacIo {
  const std::string& flac_data;
  size_t input_pos = 0;
  std::string* wav_data;
  uint32_t sample_rate;
  uint32_t channels;
  uint32_t bits_per_sample;
};

void FlacLogError(const FLAC__StreamDecoder* decoder,
                  FLAC__StreamDecoderErrorStatus status, void* client_data) {
  fprintf(stderr, "decode error: %s\n",
          FLAC__StreamDecoderErrorStatusString[status]);
}

FLAC__StreamDecoderReadStatus FlacReadFromStream(
    const FLAC__StreamDecoder* decoder, FLAC__byte buffer[], size_t* bytes,
    void* client_data) {
  FlacIo* io = reinterpret_cast<FlacIo*>(client_data);
  if (*bytes == 0) {
    return FLAC__STREAM_DECODER_READ_STATUS_ABORT;
  }
  if (io->input_pos >= io->flac_data.size()) {
    *bytes = 0;
    return FLAC__STREAM_DECODER_READ_STATUS_END_OF_STREAM;
  }
  const size_t remaining = io->flac_data.size() - io->input_pos;
  *bytes = std::min(*bytes, remaining);
  memcpy(&buffer[0], &io->flac_data[io->input_pos], *bytes);
  io->input_pos += *bytes;
  return FLAC__STREAM_DECODER_READ_STATUS_CONTINUE;
}

void FillWavHeader(const FLAC__StreamMetadata_StreamInfo& streaminfo,
                   WavHeader* wav_header) {
  memcpy(wav_header->riff_chunk.riff_chunk_id, "RIFF", 4);
  wav_header->riff_chunk.riff_chunk_size = 0;  // Filled in later.
  memcpy(wav_header->riff_chunk.wave_format, "WAVE", 4);
  memcpy(wav_header->format_chunk.format_chunk_id, "fmt ", 4);
  wav_header->format_chunk.format_chunk_size = 16;
  wav_header->format_chunk.audio_format = 1;
  wav_header->format_chunk.number_of_channels = streaminfo.channels;
  wav_header->format_chunk.sampling_frequency = streaminfo.sample_rate;
  wav_header->format_chunk.bits_per_sample = streaminfo.bits_per_sample;
  wav_header->format_chunk.block_align =
      wav_header->format_chunk.number_of_channels *
      (wav_header->format_chunk.bits_per_sample / 8);
  wav_header->format_chunk.byte_rate =
      wav_header->format_chunk.sampling_frequency *
      wav_header->format_chunk.number_of_channels *
      (wav_header->format_chunk.bits_per_sample / 8);
  memcpy(wav_header->channel_id, "data", 4);
  wav_header->channel_data_length = 0;  // Filled in later.
}

void FlacWriteMetadata(const FLAC__StreamDecoder* decoder,
                       const FLAC__StreamMetadata* metadata,
                       void* client_data) {
  FlacIo* io = reinterpret_cast<FlacIo*>(client_data);
  if (metadata->type == FLAC__METADATA_TYPE_STREAMINFO) {
    const auto& si = metadata->data.stream_info;
    io->sample_rate = si.sample_rate;
    io->channels = si.channels;
    io->bits_per_sample = si.bits_per_sample;
    WavHeader wav_header;
    FillWavHeader(si, &wav_header);
    WriteWavHeader(wav_header, io->wav_data);
  }
}

FLAC__StreamDecoderWriteStatus FlacWriteSamples(
    const FLAC__StreamDecoder* decoder, const FLAC__Frame* frame,
    const FLAC__int32* const buffer[], void* client_data) {
  FlacIo* io = reinterpret_cast<FlacIo*>(client_data);
  if (io->sample_rate != frame->header.sample_rate ||
      io->channels != frame->header.channels ||
      io->bits_per_sample != frame->header.bits_per_sample) {
    fprintf(stderr, "Frame header mismatch\n");
    return FLAC__STREAM_DECODER_WRITE_STATUS_ABORT;
  }
  uint32_t samples = frame->header.blocksize;
  std::vector<int32_t> output(samples * io->channels);
  io->wav_data->reserve(io->wav_data->size() + output.size() * sizeof(int16_t));
  for (uint32_t i = 0; i < samples; ++i) {
    for (int c = 0; c < io->channels; ++c) {
      output[i * io->channels + c] = buffer[c][i];
    }
  }
  WriteSamples(output.data(), output.size(), io->wav_data);
  return FLAC__STREAM_DECODER_WRITE_STATUS_CONTINUE;
}

bool FlacDecompress(const std::string& flac_data, std::string* wav_data) {
  std::unique_ptr<FLAC__StreamDecoder, void (*)(FLAC__StreamDecoder*)> dec(
      FLAC__stream_decoder_new(), FLAC__stream_decoder_delete);
  FlacIo io = {flac_data, 0, wav_data};
  FLAC__StreamDecoderInitStatus status = FLAC__stream_decoder_init_stream(
      dec.get(), FlacReadFromStream, nullptr, nullptr, nullptr, nullptr,
      FlacWriteSamples, FlacWriteMetadata, FlacLogError, &io);
  if (status != FLAC__STREAM_DECODER_INIT_STATUS_OK) {
    fprintf(stderr, "Failed to init flac decoder: %s\n",
            FLAC__StreamDecoderInitStatusString[status]);
    return false;
  }
  bool done = false;
  FLAC__StreamDecoderState state = FLAC__stream_decoder_get_state(dec.get());
  while (!done) {
    if (!FLAC__stream_decoder_process_single(dec.get())) {
      done = true;
    }
    state = FLAC__stream_decoder_get_state(dec.get());
    if (state == FLAC__STREAM_DECODER_END_OF_STREAM) {
      done = true;
    }
  }
  if (state != FLAC__STREAM_DECODER_END_OF_STREAM) {
    fprintf(stderr, "Error decoding flac stream: %s\n",
            FLAC__StreamDecoderStateString[state]);
    return false;
  }
  if (!FLAC__stream_decoder_finish(dec.get())) {
    fprintf(stderr, "MD5 checksum mismatch\n");
    return false;
  }
  // Fill in the length data in the wav header.
  uint32_t riff_len = wav_data->size() - 8;
  uint8_t* wav_out = reinterpret_cast<uint8_t*>(wav_data->data());
  memcpy(&wav_out[4], &riff_len, 4);
  uint32_t data_len = wav_data->size() - kWavHeaderSize;
  memcpy(&wav_out[kWavHeaderSize - 4], &data_len, 4);
  return true;
}

}  // namespace

void StreamingFlacEncoder::Reset() {
  wav_data_.clear();
  flac_data_.clear();
  output_pos_ = 0;
}

bool StreamingFlacEncoder::ProcessInput(const uint8_t* data, size_t len) {
  wav_data_.append(reinterpret_cast<const char*>(data), len);
  return true;
}

bool StreamingFlacEncoder::Flush() {
  return FlacCompress(wav_data_, compression_level_, &flac_data_);
}

size_t StreamingFlacEncoder::OutputSize() const {
  return flac_data_.size() - output_pos_;
}

size_t StreamingFlacEncoder::CopyOutput(uint8_t* buffer, size_t len) {
  size_t nbytes = std::min(len, OutputSize());
  memcpy(buffer, reinterpret_cast<const uint8_t*>(&flac_data_[output_pos_]),
         nbytes);
  output_pos_ += nbytes;
  return nbytes;
}

void StreamingFlacDecoder::Reset() {
  flac_data_.clear();
  wav_data_.clear();
  output_pos_ = 0;
}

bool StreamingFlacDecoder::ProcessInput(const uint8_t* data, size_t len) {
  flac_data_.append(reinterpret_cast<const char*>(data), len);
  return true;
}

bool StreamingFlacDecoder::Flush() {
  return FlacDecompress(flac_data_, &wav_data_);
}

size_t StreamingFlacDecoder::OutputSize() const {
  return wav_data_.size() - output_pos_;
}

size_t StreamingFlacDecoder::CopyOutput(uint8_t* buffer, size_t len) {
  size_t nbytes = std::min(len, OutputSize());
  memcpy(buffer, reinterpret_cast<const uint8_t*>(&wav_data_[output_pos_]),
         nbytes);
  output_pos_ += nbytes;
  return nbytes;
}

}  // namespace ringli
