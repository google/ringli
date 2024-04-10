FetchContent_Declare(googletest
     EXCLUDE_FROM_ALL
     GIT_REPOSITORY https://github.com/google/googletest.git
     GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(googletest)

file(GLOB_RECURSE googletest_files ${CMAKE_CURRENT_BINARY_DIR}/_deps/googletest-src *.cc *.c *.h)
set_source_files_properties(
    ${googletest_files}
    TARGET_DIRECTORY gmock gtest
    PROPERTIES SKIP_LINTING ON
)
