add_library(analysis STATIC
    audio_codec.h
    ffmpeg_codec.h
    flac_codec.cc
    flac_codec.h
    opus_codec.cc
    opus_codec.h
    ringli_codec.cc
    ringli_codec.h
    subprocess.cc
    subprocess.h
)

target_link_libraries(analysis common encode decode absl::span absl::strings absl::log hwy opus FLAC samplerate visqol)

add_executable(ringli_eval
    ringli_eval.cc
    visqol.cc
    visqol.h
    warp_q.cc
    warp_q.h
    resample.h
)

target_link_libraries(ringli_eval analysis absl::flags_parse visqol zimtohrli_base)

add_executable(ringli_analysis_test
    subprocess_test.cc
    opus_codec_compatibility_test.cc
    opus_codec_test.cc
    ringli_codec_test.cc
    generate_wav.h
    generate_wav.cc
)

target_link_libraries(ringli_analysis_test common gtest gmock_main analysis)
target_compile_definitions(ringli_analysis_test PRIVATE CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})

gtest_discover_tests(ringli_analysis_test)