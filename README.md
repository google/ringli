[![Tests](https://github.com/google/ringli/workflows/Test%20Ringli/badge.svg)](https://github.com/google/ringli/actions)

# Ringli music compression codec.

A high quality media compression codec.

## Compatibility

Ringli is a project under development, and is built and tested in a Debian-like environment.

## Build

Many dependencies for Ringli are downloaded and managed by the build script, but a few
need to be installed beforehand:

- ffmpeg
- libogg-dev
- cmake
- ninja-build
- xxd
- clang
- clang-tidy

To install these in a Debian-like system:

```
sudo apt install -y ffmpeg libogg-dev cmake ninja-build xxd clang clang-tidy
```

Once they are installed, configure the project:

```
./configure.sh
```

Build the project:
```
(cd build && ninja)
```

### Address sanitizer build

To build with address sanitizer, configure a new build directory with asan configured:


```
./configure.sh asan
```

Build the project:
```
(cd asan_build && ninja)
```

### Debug build

To build with debug symbols, configure a new build directory with debugging configured:


```
./configure.sh debug
```

Build the project:
```
(cd debug_build && ninja)
```

### Testing

```
(cd build && ninja && ninja test)
```

## Analysis

To analyze Ringli performance, build  and run the `ringli_eval` binary with WAV files:

```
(cd build && ninja ringli_eval && \
analysis/ringli_eval \
  --input_file=../analysis/input_files/01_symphonic.wav \
  --codecs="ringli:pc:o2-16:q5")
```

### Metrics

The [ViSQOL](https://github.com/google/visqol) perceptual metric is built into ringli_eval.

To run the [WARP-Q](https://github.com/wjassim/WARP-Q) metric, use the install script
`analysis/install_warp_q.sh`:

```
analysis/install_warp_q.sh $HOME/tmp/warp_q/
```

This will install a script that runs WARP-Q inside $HOME/tmp/warp_q that can be provided to the `ringli_eval` binary:

```
(cd build && ninja ringli_eval && \
analysis/ringli_eval \
  --warp_q_binary=$HOME/tmp/warp_q/warp_q.sh \
  --input_file=$HOME/tmp/01_symphonic.wav \
  --codecs="ringli:pc:o2-16:q5")
```