FetchContent_Declare(eigen
    EXCLUDE_FROM_ALL
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
)
set(BUILD_TESTING OFF CACHE INTERNAL "")
FetchContent_MakeAvailable(eigen)