add_library(decode STATIC
    ans_decode.cc
    ans_decode.h
    arith_decode.h
    bit_reader.cc
    bit_reader.h
    context_map_decode.cc
    context_map_decode.h
    entropy_decode.cc
    entropy_decode.h
    histogram_decode.cc
    histogram_decode.h
    huffman_decode.cc
    huffman_decode.h
    huffman_table.cc
    huffman_table.h
    noise_filtering.cc
    noise_filtering.h
    ringli_decoder.cc
    ringli_decoder.h
    ringli_input.h
)

target_link_libraries(decode PRIVATE absl::log)
target_link_libraries(decode PUBLIC Eigen3::Eigen)
