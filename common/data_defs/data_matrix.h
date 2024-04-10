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

#ifndef COMMON_DATA_DEFS_DATA_MATRIX_H_
#define COMMON_DATA_DEFS_DATA_MATRIX_H_

#include <cmath>
#include <functional>
#include <iomanip>
#include <ios>
#include <iostream>
#include <ostream>
#include <type_traits>
#include <vector>

#include "common/data_defs/data_vector.h"

namespace ringli {

template <typename T, int SIZE>
class DataMatrix {
 public:
  explicit DataMatrix(const T data[SIZE][SIZE]) {
    AllocateData();
    CopyData((const T*)data);
  }

  DataMatrix(const DataMatrix<T, SIZE>& other_matrix) {
    AllocateData();
    CopyData(other_matrix.data_);
  }

  DataMatrix(DataMatrix<T, SIZE>&& other_matrix) {
    AllocateData();
    CopyData(other_matrix.data_);
  }

  DataMatrix() {
    AllocateData();
    std::fill(data_, data_ + SIZE * SIZE, 0);
  }

  ~DataMatrix() { delete[] data_; }

  DataMatrix<T, SIZE>& operator=(const DataMatrix<T, SIZE>& other_matrix) {
    CopyData(other_matrix.data_);
    return *this;
  }

  DataVector<T, SIZE> operator*(const DataVector<T, SIZE>& other_vector) const {
    DataVector<T, SIZE> res;
    for (int i = 0; i < SIZE; i++) {
      double sum = 0.0;
      for (int j = 0; j < SIZE; j++) {
        sum += (*this)[i][j] * other_vector[j];
      }
      res[i] = sum;
    }
    return res;
  }

  static DataMatrix<T, SIZE> Identity() {
    DataMatrix<T, SIZE> res;
    for (int i = 0; i < SIZE; i++) {
      for (int j = 0; j < SIZE; j++) {
        res[i][j] = i == j ? 1 : 0;
      }
    }
    return res;
  }

  DataMatrix<T, SIZE> operator*(const DataMatrix<T, SIZE>& other_matrix) const {
    DataMatrix<T, SIZE> res;
    for (int i = 0; i < SIZE; i++) {
      for (int j = 0; j < SIZE; j++) {
        T sum = 0;
        for (int k = 0; k < SIZE; k++) {
          sum += (*this)[i][k] * other_matrix[k][j];
        }
        res[i][j] = sum;
      }
    }
    return res;
  }

  DataMatrix<T, SIZE> Transposed() const {
    DataMatrix<T, SIZE> res;
    for (int i = 0; i < SIZE; i++) {
      for (int j = 0; j < SIZE; j++) {
        res[i][j] = (*this)[j][i];
      }
    }
    return res;
  }

  void PrettyPrint() const {
    std::cout << std::fixed;
    std::cout << std::setprecision(2);
    for (int i = 0; i < SIZE; i++) {
      for (int j = 0; j < SIZE; j++) {
        std::cout << (*this)[i][j] << "\t";
      }
      std::cout << std::endl;
    }
  }

  std::vector<std::vector<T>> ToStdVector() const {
    std::vector<std::vector<T>> res;
    for (int i = 0; i < SIZE; i++) {
      res.push_back(std::vector<T>((*this)[i], (*this)[i] + SIZE));
    }
    return res;
  }

  DataVector<T, SIZE> Row(int i) { return DataVector<T, SIZE>((*this)[i]); }

  const T* operator[](int i) const { return data_ + i * SIZE; }
  T* operator[](int i) { return data_ + i * SIZE; }

 private:
  void AllocateData() { data_ = new T[SIZE * SIZE]; }
  void CopyData(const T* data) { std::copy(data, data + SIZE * SIZE, data_); }
  T* data_ = nullptr;
};

}  // namespace ringli

#endif  // COMMON_DATA_DEFS_DATA_MATRIX_H_
