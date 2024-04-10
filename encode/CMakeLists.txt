add_library(encode STATIC
    ans_encode.cc
    ans_encode.h
    arith_encode.h
    cluster.h
    context_map_encode.cc
    context_map_encode.h
    entropy_encode.cc
    entropy_encode.h
    fast_log.h
    histogram_encode.cc
    histogram_encode.h
    huffman_tree.cc
    huffman_tree.h
    noise_shaping.cc
    noise_shaping.h
    ringli_encoder.cc
    ringli_encoder.h
    write_bits.h
)

target_link_libraries(encode absl::log common)

