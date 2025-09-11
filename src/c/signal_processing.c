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

#include <math.h>
#include <fftw3.h>
#include "globals.h"

// =========================================================================
// SIGNAL PROCESSING AND FILTERING UTILITIES
// =========================================================================
//
// Advanced numerical processing utilities for density field handling including:
// - FFT-based convolution for density smoothing
// - Direct Gaussian convolution for smaller datasets
// - Signal filtering and processing functions
// =========================================================================
// FFT METHODS AND CONVOLUTION IMPLEMENTATIONS
// =========================================================================
/**
 * @brief Applies Gaussian smoothing using FFT-based convolution (thread-safe via critical section).
 * @details Smooths a density field defined on a potentially non-uniform grid
 *          (`log_r_grid`) using FFT convolution with a Gaussian kernel of width
 *          `sigma_log` (defined in log-space). Uses zero-padding to avoid
 *          wrap-around artifacts. This is generally faster than direct convolution
 *          for large `grid_size`. Assumes log_r_grid is uniformly spaced.
 *
 * Parameters
 * ----------
 * density_grid : const double*
 *     Input density grid array (values corresponding to `log_r_grid`).
 * grid_size : int
 *     Number of points in the input grid and density arrays.
 * log_r_grid : const double*
 *     Array of logarithmic radial grid coordinates (log10(r)). Must be uniformly spaced.
 * sigma_log : double
 *     Width (standard deviation) of the Gaussian kernel in log10-space.
 * result : double*
 *     Output array (pre-allocated, size `grid_size`) for the smoothed density field.
 *
 * Returns
 * -------
 * None (populates the `result` array).
 *
 * @note Uses FFTW library for Fast Fourier Transforms (`fftw_malloc`, `fftw_plan_dft_r2c_1d`, etc.).
 *       Requires FFTW to be installed. FFTW operations are protected by `omp critical(fftw)`.
 * @note Resulting smoothed density is clamped to a minimum value of 1e-10.
 * @warning Prints errors to stderr and returns early on memory allocation failures or
 *          FFTW plan creation failures.
 * @see direct_gaussian_convolution
 * @see gaussian_convolution
 */
void fft_gaussian_convolution(const double *density_grid, int grid_size, const double *log_r_grid, double sigma_log, double *result) {
#pragma omp critical(fftw)
    {
        // Create zero-padded arrays (step 1)
        int padded_size = 2 * grid_size;

        // Allocate memory for padded input array
        double *padded_input = (double *)fftw_malloc(sizeof(double) * padded_size);
        if (!padded_input) {
            fprintf(stderr, "Error: Failed to allocate memory for padded_input\n");
            return;
        }

        double *padded_kernel = (double *)fftw_malloc(sizeof(double) * padded_size);
        if (!padded_kernel) {
            fprintf(stderr, "Error: Failed to allocate memory for padded_kernel\n");
            fftw_free(padded_input);
            return;
        }

        double *padded_output = (double *)fftw_malloc(sizeof(double) * padded_size);
        if (!padded_output) {
            fprintf(stderr, "Error: Failed to allocate memory for padded_output\n");
            fftw_free(padded_input);
            fftw_free(padded_kernel);
            return;
        }

        // Initialize padded arrays with zeros
        for (int i = 0; i < padded_size; i++) {
            padded_input[i] = 0.0;
            padded_kernel[i] = 0.0;
        }

        // Copy input data to first half of padded array
        for (int i = 0; i < grid_size; i++)
            padded_input[i] = density_grid[i];

        // Create Gaussian kernel in spatial domain
        double dlog = log_r_grid[1] - log_r_grid[0]; // Grid spacing in log space.
        double norm = 0.0;

        for (int i = 0; i < grid_size; i++) {
            // Distance in log space.
            double x = i * dlog;

            // Gaussian kernel centered at 0.
            double kernel_val = (1.0 / (sigma_log * sqrt(2.0 * M_PI))) * exp(-0.5 * (x / sigma_log) * (x / sigma_log));

            padded_kernel[i] = kernel_val;
            norm += kernel_val;
        }

        // Normalize the kernel for unit sum
        for (int i = 0; i < grid_size; i++)
            padded_kernel[i] /= norm;

        // Prepare for FFT computation (step 2)
        fftw_complex *fft_input = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (padded_size / 2 + 1));
        if (!fft_input) {
            fprintf(stderr, "Error: Failed to allocate memory for fft_input\n");
            fftw_free(padded_input);
            fftw_free(padded_kernel);
            fftw_free(padded_output);
            return;
        }

        fftw_complex *fft_kernel = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (padded_size / 2 + 1));
        if (!fft_kernel) {
            fprintf(stderr, "Error: Failed to allocate memory for fft_kernel\n");
            fftw_free(padded_input);
            fftw_free(padded_kernel);
            fftw_free(padded_output);
            fftw_free(fft_input);
            return;
        }

        fftw_complex *fft_output = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (padded_size / 2 + 1));
        if (!fft_output) {
            fprintf(stderr, "Error: Failed to allocate memory for fft_output\n");
            fftw_free(padded_input);
            fftw_free(padded_kernel);
            fftw_free(padded_output);
            fftw_free(fft_input);
            fftw_free(fft_kernel);
            return;
        }

        fftw_plan plan_forward_input = fftw_plan_dft_r2c_1d(padded_size, padded_input, fft_input, FFTW_ESTIMATE);
        if (!plan_forward_input) {
            fprintf(stderr, "Error: Failed to create forward FFTW plan for input\n");
            fftw_free(padded_input);
            fftw_free(padded_kernel);
            fftw_free(padded_output);
            fftw_free(fft_input);
            fftw_free(fft_kernel);
            fftw_free(fft_output);
            return;
        }

        fftw_plan plan_forward_kernel = fftw_plan_dft_r2c_1d(padded_size, padded_kernel, fft_kernel, FFTW_ESTIMATE);
        if (!plan_forward_kernel) {
            fprintf(stderr, "Error: Failed to create forward FFTW plan for kernel\n");
            fftw_destroy_plan(plan_forward_input);
            fftw_free(padded_input);
            fftw_free(padded_kernel);
            fftw_free(padded_output);
            fftw_free(fft_input);
            fftw_free(fft_kernel);
            fftw_free(fft_output);
            return;
        }

        fftw_plan plan_backward = fftw_plan_dft_c2r_1d(padded_size, fft_output, padded_output, FFTW_ESTIMATE);
        if (!plan_backward) {
            fprintf(stderr, "Error: Failed to create backward FFTW plan\n");
            fftw_destroy_plan(plan_forward_input);
            fftw_destroy_plan(plan_forward_kernel);
            fftw_free(padded_input);
            fftw_free(padded_kernel);
            fftw_free(padded_output);
            fftw_free(fft_input);
            fftw_free(fft_kernel);
            fftw_free(fft_output);
            return;
        }

        // Execute forward FFTs.
        fftw_execute(plan_forward_input);
        fftw_execute(plan_forward_kernel);

        // Perform complex multiplication in frequency domain (convolution in spatial domain)
        // For each frequency component, multiply signal and kernel transforms
        for (int i = 0; i < padded_size / 2 + 1; i++) {
            double re_in = fft_input[i][0];   // Real part of input transform
            double im_in = fft_input[i][1];   // Imaginary part of input transform
            double re_ker = fft_kernel[i][0]; // Real part of kernel transform
            double im_ker = fft_kernel[i][1]; // Imaginary part of kernel transform

            // Complex multiplication: (a+bi)(c+di) = (ac-bd) + (ad+bc)i
            fft_output[i][0] = re_in * re_ker - im_in * im_ker; // Real component
            fft_output[i][1] = re_in * im_ker + im_in * re_ker; // Imaginary component
        }

        // Execute inverse FFT.
        fftw_execute(plan_backward);

        // Apply normalization to compensate for FFTW's unnormalized inverse transform
        // FFTW's implementation requires division by array length to get properly scaled result
        for (int i = 0; i < padded_size; i++)
            padded_output[i] /= padded_size; // Scale by 1/N to get normalized values

        // Extract the valid portion of the convolution result
        // Only the first grid_size elements contain the actual result (rest is padding)
        for (int i = 0; i < grid_size; i++) {
            result[i] = padded_output[i];

            // Enforce minimum density threshold to avoid numerical instability
            // Consistent with minimum threshold in direct convolution method
            if (result[i] < 1e-10)
                result[i] = 1e-10;
        }

        // Release FFTW resources
        fftw_destroy_plan(plan_forward_input);
        fftw_destroy_plan(plan_forward_kernel);
        fftw_destroy_plan(plan_backward);

        fftw_free(padded_input);
        fftw_free(padded_kernel);
        fftw_free(padded_output);
        fftw_free(fft_input);
        fftw_free(fft_kernel);
        fftw_free(fft_output);
    } // End of critical section.
}

/**
 * @brief Applies Gaussian smoothing using direct convolution.
 * @details Smooths a density field defined on a potentially non-uniform grid
 *          (`log_r_grid`) using direct spatial convolution with a Gaussian kernel
 *          of width `sigma_log` (defined in log-space).
 *          For each point `i` in the output `result` array, it computes a weighted
 *          sum of the input `density_grid` values:
 *          `result[i] = sum(density_grid[j] * kernel(log_r_grid[i] - log_r_grid[j])) / sum(kernel(...))`
 *          where the `kernel` is a Gaussian function `G(x) = (1/(σ√2π)) * exp(-0.5*(x/σ)²)`,
 *          with `x` being the distance `log_r_grid[i] - log_r_grid[j]` and `σ = sigma_log`.
 *          This method is generally more accurate than FFT-based convolution, especially
 *          for non-uniform grids or near boundaries, but has a higher computational
 *          cost (O(N²)) which makes it slower for large `grid_size`.
 *
 * Parameters
 * ----------
 * density_grid : const double*
 *     Input density grid array (values corresponding to `log_r_grid`).
 * grid_size : int
 *     Number of points in the input grid and density arrays.
 * log_r_grid : const double*
 *     Array of logarithmic radial grid coordinates (log10(r)). Can be non-uniformly spaced.
 * sigma_log : double
 *     Width (standard deviation) of the Gaussian kernel in log10-space.
 * result : double*
 *     Output array (pre-allocated, size `grid_size`) for the smoothed density field.
 *
 * Returns
 * -------
 * None (populates the `result` array).
 *
 * @note Resulting smoothed density is clamped to a minimum value of 1e-10.
 *       This function is inherently thread-safe as it only reads inputs and writes
 *       to distinct elements of the output array without shared intermediate state.
 * @see fft_gaussian_convolution
 * @see gaussian_convolution
 */
void direct_gaussian_convolution(
    const double *density_grid,
    int grid_size,
    const double *log_r_grid,
    double sigma_log,
    double *result)
{
    // Precompute normalization factor for Gaussian kernel if sigma is valid
    double kernel_norm_factor = (sigma_log > 1e-15) ? (1.0 / (sigma_log * sqrt(2.0 * M_PI))) : 1.0;
    double sig_sq_inv = (sigma_log > 1e-15) ? (1.0 / (sigma_log * sigma_log)) : 0.0; // Avoid division by zero

    // Perform direct convolution
    for (int i = 0; i < grid_size; i++) {
        double sum = 0.0;
        double norm = 0.0;

        for (int j = 0; j < grid_size; j++) {
            double dlog_r = log_r_grid[i] - log_r_grid[j];
            double kernel = (sigma_log > 1e-15) ?
                (kernel_norm_factor * exp(-0.5 * dlog_r * dlog_r * sig_sq_inv)) :
                ((i == j) ? 1.0 : 0.0); // Delta function if sigma=0

            sum += density_grid[j] * kernel;
            norm += kernel;
        }

        // Normalize the result, handle potential division by zero if norm is too small
        result[i] = (norm > 1e-15) ? (sum / norm) : density_grid[i]; // Fallback to original value if norm is zero

        // Clamp to minimum density
        if (result[i] < 1e-10)
            result[i] = 1e-10;
    }
}

/**
 * @brief Performs Gaussian smoothing by selecting the appropriate convolution method.
 * @details Acts as a routing function that delegates density smoothing to either
 *          `direct_gaussian_convolution` or `fft_gaussian_convolution` based on the
 *          value of the global `debug_direct_convolution` flag.
 *          - If `debug_direct_convolution` is non-zero, `direct_gaussian_convolution`
 *            is called (more accurate, O(N²) complexity, suitable for smaller or
 *            non-uniform grids).
 *          - If `debug_direct_convolution` is zero (default), `fft_gaussian_convolution`
 *            is called (faster for large grids, O(N log N) complexity, requires
 *            uniformly spaced logarithmic grid).
 *
 * Parameters
 * ----------
 * density_grid : const double*
 *     Input density grid array (values corresponding to `log_r_grid`).
 * grid_size : int
 *     Number of points in the input grid and density arrays.
 * log_r_grid : const double*
 *     Array of logarithmic radial grid coordinates (log10(r)). Must be uniform if FFT is used.
 * sigma_log : double
 *     Width (standard deviation) of the Gaussian kernel in log10-space.
 * result : double*
 *     Output array (pre-allocated, size `grid_size`) for the smoothed density field.
 *
 * Returns
 * -------
 * None (populates the `result` array).
 *
 * @see direct_gaussian_convolution
 * @see fft_gaussian_convolution
 * @see debug_direct_convolution
 */
void gaussian_convolution(const double *density_grid, int grid_size, const double *log_r_grid, double sigma_log, double *result) {
    // Select convolution method based on global configuration flag.
    if (debug_direct_convolution) // Use direct spatial-domain convolution.
        direct_gaussian_convolution(density_grid, grid_size, log_r_grid, sigma_log, result);
    else // Use FFT-based frequency-domain convolution.
        fft_gaussian_convolution(density_grid, grid_size, log_r_grid, sigma_log, result);
}
