FetchContent_Declare(opus
    EXCLUDE_FROM_ALL
    GIT_REPOSITORY https://github.com/xiph/opus.git 
    GIT_TAG v1.5.1
)
FetchContent_MakeAvailable(opus)