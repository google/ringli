#!/bin/sh

# Copyright 2024 The Ringli Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

if [ "${1}" = "debug" ]; then
    mkdir -p debug_build
    (cd debug_build && cmake -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ ..)
elif [ "${1}" = "asan" ]; then
    mkdir -p asan_build
    (cd asan_build && cmake -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_C_FLAGS='-fsanitize=address' -DCMAKE_CXX_FLAGS='-fsanitize=address' -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ ..)
else
    mkdir -p build
    (cd build && cmake -G Ninja -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ ..)
fi
