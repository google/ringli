add_library(common STATIC
    adaptive_quant.cc
    adaptive_quant.h
    ans_params.h
    block_predictor.cc
    block_predictor.h
    context.h
    convolve.h
    covariance_lattice.h
    dct.cc
    dct.h
    distributions.h
    entropy_coding.h
    error_norm.cc
    error_norm.h
    fast_online_predictor.cc
    fast_online_predictor.h
    log2floor.h
    logging.cc
    logging.h
    online_predictor.cc
    online_predictor.h
    predictor.h
    ringli_header.h
    segment_curve.cc
    segment_curve.h
    streaming.h
    wav_header.h
    wav_reader.cc
    wav_reader.h
    wav_writer.cc
    wav_writer.h
    data_defs/constants.h
    data_defs/data_matrix.h
    data_defs/data_vector.h
)

target_link_libraries(common Eigen3::Eigen absl::log hwy absl::log_internal_check_impl)

add_executable(ringli_common_test
    adaptive_quant_test.cc
    block_predictor_test.cc
    dct_test.cc
    online_predictor_test.cc
    segment_curve_test.cc
    data_defs/data_matrix_test.cc
    data_defs/data_vector_test.cc
)

target_link_libraries(ringli_common_test common absl::span gtest gmock_main Eigen3::Eigen)

gtest_discover_tests(ringli_common_test)
