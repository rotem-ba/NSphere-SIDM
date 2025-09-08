/*
 * Copyright 2025 Kris Sigurdson
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

 #ifndef SIGNAL_PROCESSING_H
 #define SIGNAL_PROCESSING_H

 #include <fftw3.h>

 void fft_gaussian_convolution(const double *density_grid, int grid_size, const double *log_r_grid, double sigma_log, double *result);
 void direct_gaussian_convolution(const double *density_grid, int grid_size, const double *log_r_grid, double sigma_log, double *result);
 void gaussian_convolution(const double *density_grid, int grid_size, const double *log_r_grid, double sigma_log, double *result);

 #endif // SIGNAL_PROCESSING_H
