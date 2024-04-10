// Copyright 2024 The Ringli Authors. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef COMMON_DATA_DEFS_DATA_VECTOR_H_
#define COMMON_DATA_DEFS_DATA_VECTOR_H_

#include <cmath>
#include <cstring>
#include <vector>

namespace ringli {

template <typename T, int SIZE>
class DataVector {
 public:
  explicit DataVector(const T* data) {
    AllocateData();
    CopyData(data);
  }

  explicit DataVector(const std::vector<T>& data) {
    AllocateData();
    CopyData(data.data());
  }

  DataVector(const DataVector<T, SIZE>& other_vector) {
    AllocateData();
    CopyData(other_vector.data_);
  }

  DataVector(DataVector<T, SIZE>&& other_vector) {
    AllocateData();
    CopyData(other_vector.data_);
  }

  ~DataVector() { delete[] data_; }

  DataVector() {
    AllocateData();
    memset(data_, 0, SIZE * sizeof(T));
  }

  DataVector<T, SIZE>& operator=(const DataVector<T, SIZE>& other_vector) {
    CopyData(other_vector.data_);
    return *this;
  }

  bool operator==(const DataVector<T, SIZE>& other_vector) const {
    return 0 == std::memcmp(data_, other_vector.data_, SIZE * sizeof(T));
  }

  T* begin() { return data_; }
  T* end() { return data_ + SIZE; }

  const T* begin() const { return data_; }
  const T* end() const { return data_ + SIZE; }

  T* Data() { return data_; }
  const T* Data() const { return data_; }
  std::vector<T> ToStdVector() const {
    return std::vector<T>(data_, data_ + SIZE);
  }

  const T& operator[](int i) const { return data_[i]; }
  T& operator[](int i) { return data_[i]; }

  DataVector<T, SIZE> operator+(const DataVector<T, SIZE>& other) const {
    DataVector<T, SIZE> res;
    for (int i = 0; i < SIZE; ++i) {
      res[i] = data_[i] + other[i];
    }
    return res;
  }

  DataVector<T, SIZE> operator-(const DataVector<T, SIZE>& other) const {
    DataVector<T, SIZE> res;
    for (int i = 0; i < SIZE; ++i) {
      res[i] = data_[i] - other[i];
    }
    return res;
  }

  T dot(const DataVector<T, SIZE>& other) const {
    T res = 0;
    for (int i = 0; i < SIZE; ++i) {
      res += data_[i] * other[i];
    }
    return res;
  }

  T AbsMax() const {
    T res = 0;
    for (int i = 0; i < SIZE; ++i) {
      res = std::max(res, std::abs(data_[i]));
    }
    return res;
  }

 private:
  void AllocateData() { data_ = new T[SIZE]; }
  void CopyData(const T* data) { std::copy(data, data + SIZE, data_); }
  T* data_ = nullptr;
};

template <typename T, int SIZE>
class DataVectorPack {
 public:
  explicit DataVectorPack(std::vector<ringli::DataVector<T, SIZE>>& channels)
      : channels_(channels) {}
  explicit DataVectorPack(int n_channels) : channels_(n_channels) {}

  const std::vector<ringli::DataVector<T, SIZE>>& GetChannels() const {
    return channels_;
  }
  std::vector<ringli::DataVector<T, SIZE>>& GetChannels() { return channels_; }

  ringli::DataVector<T, SIZE>& operator[](int i) { return channels_[i]; }
  const ringli::DataVector<T, SIZE>& operator[](int i) const {
    return channels_[i];
  }

 private:
  std::vector<ringli::DataVector<T, SIZE>> channels_;
};

}  // namespace ringli

#endif  // COMMON_DATA_DEFS_DATA_VECTOR_H_
