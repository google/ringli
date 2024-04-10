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

#ifndef COMMON_DATA_DEFS_CONSTANTS_H_
#define COMMON_DATA_DEFS_CONSTANTS_H_

#include <cstddef>
#include <cstdint>
#include <limits>

#include "common/data_defs/data_vector.h"
namespace ringli {

constexpr char kRingliId[8] = "RINGLI ";

constexpr size_t kRingliBlockSize = 1024;
constexpr size_t kOnlinePredictorBufferSize = 512;
constexpr size_t kMaxPredictorOrder = 16;
constexpr size_t kOnlinePredictorOrder = 8;
constexpr float kLSFQuant[kMaxPredictorOrder] = {
    160.0, 128.0, 128.0, 128.0, 96.0, 96.0, 96.0, 96.0,
    96.0,  96.0,  96.0,  96.0,  96.0, 96.0, 96.0, 96.0,
};
constexpr float kOnlinePredictorRegulariser = 1000.0;
// Optimal order 8 coefficients precomputed over 40 tracks
constexpr float kOnlinePredictorDefaultCoeffs[kOnlinePredictorOrder] = {
    2.3136797f, -2.549989f,   2.1421785f,  -1.5938007f,
    1.1497451f, -0.76307625f, 0.38299596f, -0.10532121f};
// Optimal order 2 coefficients precomputed over 40 tracks
constexpr float kOnlineO2PredictorDefaultCoeffs[2] = {1.7549847f, 0.77578706f};

constexpr size_t kShapingFilterOrder = 10;
// Coefficients obtained via least squares to match a modified ISO 226 loudness
// curve at 40 dB SPL. These will be convolved with previous samples'
// quantization noise.
constexpr float kShapingO1Coeffs[1] = {-0.55f};
constexpr float kShapingCoeffs[kShapingFilterOrder] = {
    -2.4323456f, 3.523993f, -4.6548705f, 5.494728f,   -5.34814f,
    4.5529437f,  -3.44325f, 2.1963468f,  -1.1031115f, 0.32658848f};

constexpr size_t kNoiseFilterOrder = 17;
// Coefficients obtained via least squares for a linear phase FIR lowpass
// design with cutoff frequency at 15500 Hz and  1500Hz transition width.
// 50 times more loss weight was put on reducing pass-band ripple than
// rejection band attenuation.
constexpr float kNoiseFilterCoeffs[kNoiseFilterOrder] = {
    0.00131554f,  -0.02018712f, 0.0468852f,  -0.05967534f, 0.03253596f,
    0.04501965f,  -0.15389215f, 0.25015644f, 0.71139493f,  0.25015644f,
    -0.15389215f, 0.04501965f,  0.03253596f, -0.05967534f, 0.0468852f,
    -0.02018712f, 0.00131554f};

// Our IIR filters are going to deal with 4 filters of order 3
constexpr int kNumIIRFilters = 2;
constexpr int kIIROrder = 3;

typedef float FilterCoeffArray[kNumIIRFilters][kIIROrder];

// IIR filter coefficients for adaptive quantization band power masking
// estimation
// Filter passabnd centers: [10400, 14000]
// Passband centers in barks: [ 23.7, 24.75]
constexpr FilterCoeffArray kBwdCoeffs = {{0.07061465f, 0, -0.07061465f},
                                         {0.05607456f, 0, -0.05607456f}};

constexpr FilterCoeffArray kFwdCoeffs = {{1.0f, -0.38886467f, 0.8587707f},
                                         {1.0f, 0.4892825f, 0.8878509f}};

constexpr float kEmaMomentum = 0.97f;

// coefficient that relates the power measured by each filter and the maximum
// unshaped noise power acceptable for the respective band
// These values are computed heuristically taking into account the shapes of our
// noise shaping filter and masking power estimation filters
constexpr float kMaskingFilterPowerSpreadRatio[2] = {0.3758374f, 1.4125376f};

// amount we expect original signal to mask in-band
constexpr float kMaskingRatio = 3.1623f;  // std::pow(10, 5 dB / 10)

constexpr float kMaskingFilterTargetNSR[kNumIIRFilters] = {
    kMaskingFilterPowerSpreadRatio[0] / kMaskingRatio,
    kMaskingFilterPowerSpreadRatio[1] / kMaskingRatio,
};

constexpr int kHQSteps = 6;
constexpr int kAdaptiveQuantRampUpSteps = 16;
constexpr int kAdaptiveQuantRampDownSteps = 8;

constexpr float kQuantMax =
    650;                         //  corresponds to an unshaped PSNR of 50.36 dB
constexpr float kQuantMin = 30;  // corresponds to an unshaped PSNR of 77 dB
constexpr float kQuantStart = 35;  // corresponds to an unshaped PSNR of 75.6 dB

constexpr size_t kDctLength = 64;
constexpr size_t kDctNumber = kRingliBlockSize / kDctLength;

constexpr int kNumDctBands = 3;
constexpr int kDctBandEnd[kNumDctBands] = {8, 20, kDctLength};

typedef DataVector<int32_t, kRingliBlockSize> RingliVector;

typedef DataVectorPack<int32_t, kRingliBlockSize> AudioBlock;

constexpr int kNumACPredictionSteps = 3;
constexpr int kACPredictionStart = 4;
constexpr double kACPredictionSigma = 16.0f;
constexpr size_t kACPredictionBorder = kNumACPredictionSteps * kDctLength;
constexpr size_t kACPredictionWindowSize =
    kRingliBlockSize + 2 * kACPredictionBorder;

}  // namespace ringli

#endif  // COMMON_DATA_DEFS_CONSTANTS_H_
