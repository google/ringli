FetchContent_Declare(highway
    EXCLUDE_FROM_ALL
    GIT_REPOSITORY https://github.com/google/highway.git
    GIT_TAG 1.1.0
)
set(HWY_ENABLE_TESTS OFF CACHE INTERNAL "")
set(BUILD_TESTING OFF CACHE INTERNAL "")
FetchContent_MakeAvailable(highway)

file(GLOB_RECURSE hwy_files ${CMAKE_CURRENT_BINARY_DIR}/_deps/highway-src *.cc *.c *.h)
set_source_files_properties(
    ${hwy_files}
    TARGET_DIRECTORY hwy
    PROPERTIES SKIP_LINTING ON
)
