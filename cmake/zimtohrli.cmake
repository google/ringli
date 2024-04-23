FetchContent_Declare(zimtohrli
    EXCLUDE_FROM_ALL
    GIT_REPOSITORY https://github.com/google/zimtohrli.git
    GIT_TAG v0.1.10
)
set(BUILD_ZIMTOHRLI_TESTS OFF CACHE INTERNAL "")
FetchContent_MakeAvailable(zimtohrli)
