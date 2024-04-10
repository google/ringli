FetchContent_Declare(flac
    EXCLUDE_FROM_ALL
    GIT_REPOSITORY https://github.com/xiph/flac.git
    GIT_TAG 1.4.3
)
set(INSTALL_MANPAGES OFF  CACHE INTERNAL "")
FetchContent_MakeAvailable(flac)