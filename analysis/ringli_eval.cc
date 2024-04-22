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

#include <sys/stat.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/log/check.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_split.h"
#include "absl/strings/string_view.h"
#include "absl/strings/substitute.h"
#include "absl/time/time.h"
#include "absl/types/span.h"
#include "analysis/audio_codec.h"
#include "analysis/flac_codec.h"
#include "analysis/opus_codec.h"
#include "analysis/ringli_codec.h"
#include "analysis/visqol.h"
#include "analysis/warp_q.h"
#include "common/data_defs/constants.h"
#include "common/error_norm.h"
#include "common/logging.h"
#include "common/streaming.h"
#include "common/wav_header.h"
#include "common/wav_reader.h"
#include "hwy/aligned_allocator.h"
#include "zimt/mos.h"
#include "zimt/zimtohrli.h"

ABSL_FLAG(int, trim_track_seconds, 0,
          "The track is first trimmed to this length. If 0 the track will not "
          "be trimmed.");
ABSL_FLAG(std::string, input_file, "", "Input file.");
ABSL_FLAG(
    std::string, output_file, "",
    "Output file. If not empty, will be populated by the roundtrip-compressed "
    "input file. If more than one codec is provided the name prefix "
    "will be appended with codec name. For example, "
    "'--codecs=ringli:apc:aconly:e6:q1,ringli:apc:aconly:e5:q6 "
    "--output_file=output.wav' will produce the files "
    "'output.ringli:apc:aconly:e6:q1.wav' "
    "and 'output.ringli:apc:aconly:e5:q6.wav'");
ABSL_FLAG(std::vector<std::string>, codecs, std::vector<std::string>(),
          "Comma-separated list of codecs to compare.");
ABSL_FLAG(bool, skip_comparison, false,
          "Whether to skip audio metrics comparison.");
ABSL_FLAG(int, encode_reps, 1, "Encoder repetitions for benchmarking");
ABSL_FLAG(int, decode_reps, 1, "Decoder repetitions for benchmarking");
ABSL_FLAG(bool, report_speed, false,
          "Whether to include speed measurements in the evaluation report");
ABSL_FLAG(std::optional<std::string>, warp_q_binary, std::nullopt,
          "Path to a WARP-Q binary to use for evaluating codecs. Run "
          "analysis/install_warp_q.sh to install the WARP-Q binary.");
ABSL_FLAG(bool, run_visqol, false,
          "Whether to run ViSQOL when analyzing results.");

namespace ringli {
namespace {

template <typename S>
int CheckFileOk(const S& fs, absl::string_view filename) {
  if (!fs.is_open() || fs.fail()) {
    std::cout << "Could not open file: " << filename << std::endl;
    return -1;
  }
  return 0;
}

std::unique_ptr<StreamingAudioCodec> CreateCodec(
    const std::string& codec_name) {
  const std::vector<std::string> params = absl::StrSplit(codec_name, ':');
  std::vector<std::unique_ptr<StreamingAudioCodec>> available_codecs;
  available_codecs.push_back(std::make_unique<StreamingRingliCodec>());
  available_codecs.push_back(std::make_unique<StreamingOpusCodec>());
  available_codecs.push_back(std::make_unique<StreamingFlacCodec>());

  const std::string& name = params[0];
  for (auto& codec : available_codecs) {
    if (codec->Name() == name) {
      CHECK(codec->ParseParams(params));
      return std::move(codec);
    }
  }
  CHECK(false) << "Unknown codec " << name;
}

std::pair<std::string, std::string> SplitFilename(const std::string& filename) {
  if (const size_t pos = filename.find_last_of('.'); pos != std::string::npos) {
    return {filename.substr(0, pos), filename.substr(pos)};
  }
  return {std::string(filename), ""};
}

std::string CodecFilename(const std::string& filename,
                          const std::string& codec_name) {
  std::pair<std::string, std::string> filename_splits = SplitFilename(filename);
  return absl::StrCat(filename_splits.first, ".", codec_name,
                      filename_splits.second);
}

WavHeader ReadWavHeaderAndTrim(std::string& input, int max_length_seconds) {
  uint8_t* data = reinterpret_cast<uint8_t*>(input.data());
  size_t pos = 0;
  WavHeader header;
  CHECK(ParseWavHeader(data, input.size(), &pos, &header));
  const FormatChunk& format = header.format_chunk;
  const size_t num_channels = format.number_of_channels;
  const size_t bytes_per_tick = num_channels * format.bits_per_sample / 8;
  const uint32_t max_data_length =
      bytes_per_tick * format.sampling_frequency * max_length_seconds;
  if (max_data_length && header.channel_data_length > max_data_length) {
    input.resize(pos + max_data_length);
    // Change the RIFF chunk size and data chunk size fields.
    const uint32_t riff_length = input.size() - 8;
    memcpy(data + 4, &riff_length, sizeof(riff_length));
    memcpy(data + pos - 4, &max_data_length, sizeof(max_data_length));
    header.riff_chunk.riff_chunk_size = riff_length;
    header.channel_data_length = max_data_length;
  }
  return header;
}

class ResultsReport {
 public:
  std::string encoding_;
  double bitrate_;
  double max_delay_;
  double encode_speed_;
  double decode_speed_;
  ErrorNorm error_norm_;
  double visqol_ = 0.0;
  double warp_q_ = 0.0;
  double zimtohrli_ = 0.0;

  static void PrintHeader(int encoding_max_length = 20) {
    PrintColumn("Encoding", encoding_max_length);
    PrintColumn("kbps");
    PrintColumn("PSNR");
    PrintColumn("P4");
    PrintColumn("Delay");
    if (absl::GetFlag(FLAGS_report_speed)) {
      PrintColumn("Enc speed");
      PrintColumn("Dec speed");
    }
    if (!absl::GetFlag(FLAGS_skip_comparison)) {
      PrintColumn("Zimtohrli");
      if (absl::GetFlag(FLAGS_run_visqol)) {
        PrintColumn("Visqol");
      }
      if (absl::GetFlag(FLAGS_warp_q_binary).has_value()) {
        PrintColumn("Warp-q");
      }
    }
    std::cout << std::endl;
  }

  void Print(int encoding_max_length = 20) const {
    PrintColumn(encoding_, encoding_max_length);
    std::cout << std::setprecision(4) << bitrate_ << "\t" << (error_norm_.psnr)
              << "\t" << (error_norm_.p4_diff) << "\t";
    if (max_delay_ < 1e-3) {
      std::cout << max_delay_ * 1e6 << "us\t";
    } else if (max_delay_ < 1) {
      std::cout << max_delay_ * 1e3 << "ms\t";
    } else {
      std::cout << max_delay_ << "s\t";
    }
    if (absl::GetFlag(FLAGS_report_speed)) {
      std::cout << encode_speed_ << "x\t\t" << decode_speed_ << "x\t\t";
    }
    if (!absl::GetFlag(FLAGS_skip_comparison)) {
      std::cout << std::setprecision(4) << zimtohrli_ << "\t\t";
      if (absl::GetFlag(FLAGS_run_visqol)) {
        std::cout << visqol_ << "\t";
      }
      if (absl::GetFlag(FLAGS_warp_q_binary).has_value()) {
        std::cout << warp_q_;
      }
    }
    std::cout << std::endl;
  }

  static void PrintReports(absl::Span<const ResultsReport> reports) {
    std::vector<int> encoding_names_length = {8};
    std::transform(
        reports.begin(), reports.end(),
        std::back_inserter(encoding_names_length),
        [](const ResultsReport& report) { return report.encoding_.size(); });
    const int max_encoding_name_length = *std::max_element(
        encoding_names_length.begin(), encoding_names_length.end());
    PrintHeader(max_encoding_name_length);
    for (const ResultsReport& report : reports) {
      report.Print(max_encoding_name_length);
    }
  }

 private:
  static void PrintColumn(absl::string_view column, int min_chars = 0) {
    std::cout << column;
    if (column.size() < min_chars) {
      for (int i = 0; i < min_chars - column.size(); ++i) {
        std::cout << " ";
      }
    }
    std::cout << "\t";
  }

  static void PrintColumn(int column, int min_chars = 0) {
    PrintColumn(std::to_string(column), min_chars);
  }
};

float ComputeSpeed(double length, const std::function<void()> f) {
  std::chrono::time_point start = std::chrono::high_resolution_clock::now();
  f();
  std::chrono::time_point stop = std::chrono::high_resolution_clock::now();
  int64_t microseconds =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start)
          .count();
  return length * 1e6 / static_cast<double>(microseconds);
}

template <typename T>
double ReadSample(const T* data) {
  int16_t sample;
  memcpy(&sample, data, sizeof(int16_t));
  return sample;
}

template <typename T>
hwy::AlignedNDArray<float, 2> BytesToSamples(const T* data, size_t num_channels,
                                             size_t num_frames) {
  const double scale = 1.0f / std::numeric_limits<int16_t>::max();
  hwy::AlignedNDArray<float, 2> result({num_channels, num_frames});
  for (size_t channel_index = 0; channel_index < num_channels;
       ++channel_index) {
    for (size_t step_index = 0; step_index < num_frames; ++step_index) {
      const double d =
          ReadSample(data + ((step_index * num_channels + channel_index)));
      result[{channel_index}][step_index] = d * scale;
    }
  }
  return result;
}

void CompressStreaming(StreamingInterface* encoder, StreamingInterface* decoder,
                       const WavHeader& wav_header, const uint8_t* input_data,
                       size_t input_size, uint8_t* output_data,
                       size_t* output_size, size_t* compressed_size = nullptr,
                       size_t* max_delay_tick = nullptr) {
  const FormatChunk& format = wav_header.format_chunk;
  const size_t num_channels = format.number_of_channels;
  const size_t bytes_per_tick = num_channels * format.bits_per_sample / 8;
  const size_t kChannelBufferSize = 1024;
  std::vector<uint8_t> channel_buffer(kChannelBufferSize);
  size_t input_tick = 0;
  size_t output_tick = 0;
  size_t input_pos = 0;
  size_t output_pos = 0;
  bool encoder_header_done = !encoder;
  bool decoder_header_done = !decoder;
  bool encoder_flushed = !encoder;
  bool decoder_flushed = !decoder;
  bool done = false;
  if (encoder && decoder) {
    CHECK(max_delay_tick);
    CHECK(compressed_size);
    *max_delay_tick = 0;
    *compressed_size = 0;
  }
  if (encoder) encoder->Reset();
  if (decoder) decoder->Reset();
  while (!done) {
    if (decoder && decoder->OutputSize() > 0) {
      CHECK_LT(output_pos, *output_size);
      size_t nbytes = decoder->CopyOutput(output_data + output_pos,
                                          *output_size - output_pos);
      output_pos += nbytes;
      if (output_pos < kWavHeaderSize) {
        continue;
      } else if (!decoder_header_done) {
        VerifyWavHeader(wav_header, output_data, kWavHeaderSize);
        decoder_header_done = true;
      } else {
        output_tick = (output_pos - kWavHeaderSize) / bytes_per_tick;
      }
    } else if (encoder && encoder->OutputSize() > 0) {
      if (decoder) {
        size_t nbytes =
            encoder->CopyOutput(channel_buffer.data(), channel_buffer.size());
        *compressed_size += nbytes;
        CHECK(decoder->ProcessInput(channel_buffer.data(), nbytes));
      } else {
        CHECK_LT(output_pos, *output_size);
        size_t nbytes = encoder->CopyOutput(output_data + output_pos,
                                            *output_size - output_pos);
        output_pos += nbytes;
      }
    } else if (!encoder_header_done) {
      CHECK(encoder->ProcessInput(input_data, wav_header.header_size));
      input_pos += wav_header.header_size;
      encoder_header_done = true;
    } else if (input_pos < input_size) {
      if (encoder) {
        CHECK(encoder->ProcessInput(input_data + input_pos, bytes_per_tick));
        input_pos += bytes_per_tick;
        ++input_tick;
        if (decoder) {
          *max_delay_tick = std::max(*max_delay_tick, input_tick - output_tick);
        }
      } else {
        size_t len = std::min(input_size - input_pos, kChannelBufferSize);
        CHECK(decoder->ProcessInput(input_data + input_pos, len));
        input_pos += len;
      }
    } else if (!encoder_flushed) {
      CHECK(encoder->Flush());
      encoder_flushed = true;
    } else if (!decoder_flushed) {
      CHECK(decoder->Flush());
      decoder_flushed = true;
    } else {
      done = true;
    }
  }
  *output_size = output_pos;
}

void EvaluateCodec(const std::string& codec_name, const WavHeader& wav_header,
                   const std::string& input, std::string& output,
                   ResultsReport& report) {
  const FormatChunk& format = wav_header.format_chunk;
  const size_t num_channels = format.number_of_channels;
  const size_t bytes_per_sample = format.bits_per_sample / 8;
  const size_t bytes_per_tick = num_channels * bytes_per_sample;
  CHECK_EQ(wav_header.channel_data_length % bytes_per_tick, 0);
  const size_t num_ticks = wav_header.channel_data_length / bytes_per_tick;
  const double input_length_s =
      wav_header.channel_data_length * 1.0 / wav_header.format_chunk.byte_rate;

  const auto codec = CreateCodec(codec_name);
  report.encoding_ = codec->ToString();

  CHECK_GE(input.size(), wav_header.header_size + num_ticks * bytes_per_tick);
  const uint8_t* const input_data = reinterpret_cast<const uint8_t*>(&input[0]);
  output.resize(input.size(), 0);
  uint8_t* const output_data = reinterpret_cast<uint8_t*>(&output[0]);
  {
    size_t output_size = output.size();
    size_t compressed_size, max_delay_tick;
    CompressStreaming(codec->encoder(), codec->decoder(), wav_header,
                      input_data, input.size(), output_data, &output_size,
                      &compressed_size, &max_delay_tick);
    output.resize(output_size);
    CHECK_EQ(output_size, kWavHeaderSize + num_ticks * bytes_per_tick);
    report.bitrate_ = compressed_size * 8.0 / 1024.0 / input_length_s;
    report.max_delay_ = max_delay_tick * 1.0 / format.sampling_frequency;
  }

  StreamComparator comparator(num_channels);
  size_t input_pos = wav_header.header_size;
  size_t output_pos = kWavHeaderSize;
  for (size_t t = 0; t < num_ticks; ++t) {
    for (int c = 0; c < num_channels; ++c) {
      double a = ReadSample(input_data + input_pos);
      double b = ReadSample(output_data + output_pos);
      comparator.AddSamples(c, a, b);
      input_pos += bytes_per_sample;
      output_pos += bytes_per_sample;
    }
  }
  report.error_norm_.psnr = comparator.PSNR();
  report.error_norm_.p4_diff = comparator.P4Diff();
  if (absl::GetFlag(FLAGS_log_level) >= 1) {
    size_t channel, pos;
    const double max_diff = comparator.MaxDiff(&channel, &pos);
    printf("Max diff: %.0f  [channel: %zu block: %zu idx: %zu]\n", max_diff,
           channel, pos / kRingliBlockSize, pos % kRingliBlockSize);
  }
  if (absl::GetFlag(FLAGS_report_speed)) {
    int log_level = absl::GetFlag(FLAGS_log_level);
    absl::SetFlag(&FLAGS_log_level, 0);
    std::vector<uint8_t> compressed(input.size());
    {
      const int reps = absl::GetFlag(FLAGS_encode_reps);
      const auto run = [&]() {
        for (int i = 0; i < reps; ++i) {
          size_t compressed_size = compressed.size();
          CompressStreaming(codec->encoder(), nullptr, wav_header, input_data,
                            input.size(), compressed.data(), &compressed_size);
          compressed.resize(compressed_size);
        }
      };
      report.encode_speed_ = ComputeSpeed(reps * input_length_s, run);
    }
    {
      const int reps = absl::GetFlag(FLAGS_decode_reps);
      const auto run = [&]() {
        for (int i = 0; i < reps; ++i) {
          size_t output_size = output.size();
          CompressStreaming(nullptr, codec->decoder(), wav_header,
                            compressed.data(), compressed.size(), output_data,
                            &output_size);
        }
      };
      report.decode_speed_ = ComputeSpeed(reps * input_length_s, run);
    }
    absl::SetFlag(&FLAGS_log_level, log_level);
  }
}

void EvaluateCodecs(const std::string& input_file,
                    const std::string& output_file, int max_length_seconds,
                    const std::vector<std::string>& codec_names) {
  std::ifstream input_stream(input_file, std::ios_base::binary);
  CHECK(input_stream.good());
  std::string input((std::istreambuf_iterator<char>(input_stream)),
                    std::istreambuf_iterator<char>());
  input_stream.close();
  const WavHeader wav_header = ReadWavHeaderAndTrim(input, max_length_seconds);
  std::vector<ResultsReport> results;
  std::unique_ptr<ViSQOL> visqol;
  if (!absl::GetFlag(FLAGS_skip_comparison) &&
      absl::GetFlag(FLAGS_run_visqol)) {
    visqol = std::make_unique<ViSQOL>();
  }
  for (const std::string& codec_name : codec_names) {
    std::string output;
    ResultsReport report;
    EvaluateCodec(codec_name, wav_header, input, output, report);
    if (!output_file.empty()) {
      const std::string output_fn = codec_names.size() > 1
                                        ? CodecFilename(output_file, codec_name)
                                        : output_file;
      std::ofstream output_stream(output_fn);
      CHECK(output_stream.good());
      output_stream.write(output.data(), output.size());
      CHECK(output_stream.good());
      output_stream.close();
      CHECK(output_stream.good());
      if (!absl::GetFlag(FLAGS_skip_comparison)) {
        hwy::AlignedNDArray<float, 2> input_samples = BytesToSamples(
            input.data(), wav_header.NumChannels(), wav_header.NumSamples());
        hwy::AlignedNDArray<float, 2> output_samples = BytesToSamples(
            output.data(), wav_header.NumChannels(), wav_header.NumSamples());
        if (visqol.get() != nullptr) {
          report.visqol_ = visqol->MOS(input_samples, output_samples,
                                       wav_header.SampleRate());
        }
        if (absl::GetFlag(FLAGS_warp_q_binary).has_value()) {
          WarpQ warp_q{.warp_q_binary =
                           absl::GetFlag(FLAGS_warp_q_binary).value()};
          report.warp_q_ = warp_q.MOS(input_file, output_file);
        }
        {
          float sum_of_squares = 0.0f;
          zimtohrli::Cam cam;
          cam.minimum_bandwidth_hz = 4;
          zimtohrli::Zimtohrli z{
              .cam_filterbank = cam.CreateFilterbank(wav_header.SampleRate()),
              .full_scale_sine_db = 80};
          const size_t num_bands = z.cam_filterbank->filter.Size();
          const size_t perceptual_sample_rate = 100;
          const size_t num_energy_frames = wav_header.NumSamples() *
                                           perceptual_sample_rate /
                                           wav_header.SampleRate();
          hwy::AlignedNDArray<float, 2> bands(
              {wav_header.NumSamples(), num_bands});
          hwy::AlignedNDArray<float, 2> energy_bands_db(
              {num_energy_frames, num_bands});
          hwy::AlignedNDArray<float, 2> in_spectrogram(
              {num_energy_frames, num_bands});
          hwy::AlignedNDArray<float, 2> out_spectrogram(
              {num_energy_frames, num_bands});
          for (size_t c = 0; c < wav_header.NumChannels(); ++c) {
            hwy::Span<float> input_signal = input_samples[{c}];
            hwy::Span<float> output_signal = output_samples[{c}];
            z.Spectrogram(input_signal, bands, energy_bands_db, in_spectrogram,
                          in_spectrogram);
            z.Spectrogram(output_signal, bands, energy_bands_db,
                          out_spectrogram, out_spectrogram);
            float distance = z.Distance(/* verbose */ false, in_spectrogram,
                                        out_spectrogram, perceptual_sample_rate)
                                 .value;
            sum_of_squares += distance * distance;
          }
          report.zimtohrli_ = zimtohrli::MOSFromZimtohrli(
              std::sqrt(sum_of_squares / wav_header.NumChannels()));
        }
      }
    }
    results.push_back(report);
  }
  std::cout << input_file << std::endl;
  ResultsReport::PrintReports(results);
}

int Main(int argc, char* argv[]) {
  absl::ParseCommandLine(argc, argv);

  const std::filesystem::path input_file = absl::GetFlag(FLAGS_input_file);
  const std::filesystem::path output_file = absl::GetFlag(FLAGS_output_file);
  if (input_file.empty()) {
    std::cout << "input_file must be specified." << std::endl;
    return 0;
  }
  const bool save_output =
      !output_file.empty() || !absl::GetFlag(FLAGS_skip_comparison);
  const std::filesystem::path actual_output_file =
      output_file.empty() && save_output
          ? std::filesystem::path(absl::StrCat(std::tmpnam(nullptr), ".wav"))
          : output_file;
  if (save_output) {
    std::filesystem::create_directories(actual_output_file.parent_path());
  }

  EvaluateCodecs(input_file, actual_output_file,
                 absl::GetFlag(FLAGS_trim_track_seconds),
                 absl::GetFlag(FLAGS_codecs));

  if (actual_output_file != output_file) {
    std::filesystem::remove_all(actual_output_file);
  }
  return 0;
}

}  // namespace
}  // namespace ringli

/* clang-format off
How to run:

 blaze run -c opt //util/compression/ringli/analysis:ringli_eval -- \
 --input_file=/google/data/rw/users/sz/szabadka/test_music/48kHz16bit/Laura.wav \
 --trim_track_seconds=5 --codecs="ringli:qb(0;1);(1000;10)"

 To resample a 44.1 kHz track to 48 kHz for comparison with opus codec, run:

 ffmpeg -y -i <input.wav> -c:a pcm_s16le -ar 48000 <output.wav>
clang-format on

Input dir options:
experimental/amik/new_codec/input
util/compression/tabuli/testdata/silph31.wav
/google/data/rw/users/jy/jyrki/hq_test_music/
*/
int main(int argc, char* argv[]) { return ringli::Main(argc, argv); }
