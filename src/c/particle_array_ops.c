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

 #include "globals.h"
 #include "particle_array_ops.h"
 #include "exit.h"
 #include "nsphere_sort.h"
 #include "globals.h"
 #include <string.h>

// =========================================================================
// PARTICLE DATA STRUCTURES AND OPERATIONS
// =========================================================================

/**
 * @brief Check if an array is strictly monotonically increasing.
 *
 * @param arr Array to check
 * @param n Number of elements
 * @param name Name of the array for debug messages
 * @return 1 if strictly monotonic, 0 otherwise
 */
int check_strict_monotonicity(const double *arr, int n, const char *name) {
    int i;
    for (i = 1; i < n; i++) {
        if (arr[i] <= arr[i-1]) {
            fprintf(stderr, "MONOTONICITY_CHECK FAILED for '%s': arr[%d]=%.17e <= arr[%d]=%.17e\n",
                       name, i, arr[i], i-1, arr[i-1]);
            fflush(stderr);
            // Print a few surrounding values for context
            for (int k = (i > 2 ? i - 2 : 0); k < (i + 3 < n ? i + 3 : n); k++) {
                fprintf(stderr, "  Context: %s[%d] = %.17e\n", name, k, arr[k]);
                fflush(stderr);
            }
            return 0;
        }
    }
    return 1;
}

/**
 * @brief Comparison function for qsort to sort PartData structures by radius.
 * @details Compares two PartData structures based on their `rad` (radius) member
 *          for sorting in ascending order. Handles NaN values by placing them
 *          consistently (e.g., at the beginning or end, behavior might depend on qsort NaN handling).
 *
 * @param a Pointer to the first PartData structure.
 * @param b Pointer to the second PartData structure.
 * @return int -1 if pa->rad < pb->rad, 1 if pa->rad > pb->rad, 0 otherwise.
 *             Specific return for NaNs ensures consistent ordering.
 */
int compare_partdata_by_rad(const void *a, const void *b)
{
    if (!a && !b) return 0;
    if (!a) return -1;
    if (!b) return 1;

    const struct PartData *pa = (const struct PartData *)a;
    const struct PartData *pb = (const struct PartData *)b;

    int pa_is_nan = (pa->rad != pa->rad); // Check for NaN (safe with fast-math)
    int pb_is_nan = (pb->rad != pb->rad); // Check for NaN (safe with fast-math)

    if (pa_is_nan && pb_is_nan) return 0;
    if (pa_is_nan) return -1;
    if (pb_is_nan) return 1;

    if (pa->rad < pb->rad) return -1;
    if (pa->rad > pb->rad) return 1;
    return 0;
}

/**
 * @brief Comparison function for qsort to sort RrPsiPair structures by the 'rr' member.
 * @details Used to sort an array of RrPsiPair structures in ascending order
 *          based on their radial (`rr`) component. This is primarily used when
 *          preparing data for splines where the x-axis (e.g., radius or -Psi)
 *          must be strictly monotonic.
 *
 * @param a Pointer to the first RrPsiPair structure.
 * @param b Pointer to the second RrPsiPair structure.
 * @return int -1 if pa->rr < pb->rr, 1 if pa->rr > pb->rr, 0 otherwise.
 */
int compare_by_rr(const void *a, const void *b)
{
    const struct RrPsiPair *pa = (const struct RrPsiPair *)a;
    const struct RrPsiPair *pb = (const struct RrPsiPair *)b;
    if (pa->rr < pb->rr)
        return -1;
    if (pa->rr > pb->rr)
        return 1;
    return 0;
}

/**
 * @brief Sorts an array of PartData structures by their radial position (`rad`).
 * @details Uses the standard library `qsort` function with `compare_partdata_by_rad`
 *          as the comparison function. Includes basic safety checks for NULL array
 *          or non-positive `npts`. The sort is performed in-place.
 *
 * @param array [in,out] Array of PartData structures to be sorted.
 * @param npts  [in] Number of elements in the array.
 */
void sort_by_rad(struct PartData *array, int npts)
{
    if (!array)
    {
        fprintf(stderr, "ERROR: sort_by_rad called with NULL array\n");
        return;
    }
    if (npts <= 0)
    {
        // Sorting an empty or negatively sized array is meaningless or an error.
        // fprintf(stderr, "Warning: sort_by_rad called with npts <= 0: %d\n", npts);
        return; // Nothing to sort
    }

    qsort(array, (size_t)npts, sizeof(struct PartData), compare_partdata_by_rad);
}

/**
 * @brief Comparison function for sorting double values in ascending order.
 * @details This function is designed to be used with `qsort` or other
 *          standard library sorting functions that require a comparator.
 *          It takes two void pointers, casts them to `const double*`,
 *          dereferences them, and compares their values.
 *
 * @param a [in] Pointer to the first double value.
 * @param b [in] Pointer to the second double value.
 * @return int -1 if the first double is less than the second,
 *              1 if the first double is greater than the second,
 *              0 if they are equal.
 */
int double_cmp(const void *a, const void *b)
{
    double da = *(const double *)a;
    double db = *(const double *)b;

    if (da < db)
        return -1;
    if (da > db)
        return 1;
    return 0;
}

// =========================================================================
// PARALLEL SORTING ALGORITHM FUNCTIONS
// =========================================================================
//
// Functions implementing parallel sorting algorithms.
// The sorting implementation uses a parallel chunk-based approach with overlapping
// regions between chunks to ensure correct ordering at chunk boundaries.
// Constants controlling behavior are defined at the top of this file.

/**
 * @brief Compares two particle entries for sorting based on the first column value.
 *
 * @details Used as a comparison function for qsort, quadsort, and other sorting algorithms.
 *          Particles with smaller values in column 0 will be sorted before those with larger values.
 *
 * @param a Pointer to the first particle entry (as void pointer, expected double**).
 * @param b Pointer to the second particle entry (as void pointer, expected double**).
 * @return -1 if a<b, 1 if a>b, 0 if equal based on the first column value (radius).
 *
 * @note Assumes the particle data is structured as [particle_index][component_index]
 *       when passed via `columns` array in sorting functions, and compares `columns[i][0]`.
 */
int compare_particles(const void *a, const void *b)
{
    double *col_a = *(double **)a;
    double *col_b = *(double **)b;
    if (col_a[0] < col_b[0])
        return -1;
    if (col_a[0] > col_b[0])
        return 1;
    return 0;
}

/**
 * @brief Implementation of classic insertion sort algorithm for particle data.
 *
 * @details Sorts an array of particle data using the insertion sort algorithm.
 *          While not the fastest algorithm for large datasets, it is stable and works well
 *          for small arrays or nearly sorted data.
 *
 * @param columns 2D array of particle data to be sorted
 * @param n Number of elements to sort
 */
void insertion_sort(double **columns, int n)
{
    for (int i = 1; i < n; i++)
    {
        double *temp = columns[i];
        int j = i - 1;
        while (j >= 0 && compare_particles(&columns[j], &temp) > 0)
        {
            columns[j + 1] = columns[j];
            j--;
        }
        columns[j + 1] = temp;
    }
}

/**
 * @brief Wrapper for the standard C library quicksort function.
 *
 * @details Provides a consistent interface to the standard library `qsort` function
 *          using the `compare_particles` function as the comparison callback.
 *
 * @param columns 2D array of particle data to be sorted (passed as `double**`).
 * @param n Number of elements (particles) to sort.
 */
void stdlib_qsort_wrapper(double **columns, int n)
{
    qsort(columns, n, sizeof(double *), compare_particles);
}

/**
 * Wrapper for the external quadsort algorithm.
 *
 * Provides a consistent interface to the external quadsort function,
 * which is typically faster than standard quicksort for many distributions.
 *
 * @param columns 2D array of particle data to be sorted
 * @param n Number of elements to sort
 */
void quadsort_wrapper(double **columns, int n)
{
    quadsort(columns, n, sizeof(double *), compare_particles);
}

/**
 * @brief Helper function that performs insertion sort on a subarray of particle data.
 *
 * @details Similar to `insertion_sort` but operates on a specific range `[start..end]` inclusive.
 *          Used by the parallel sorting algorithms to sort individual chunks and seam regions.
 *
 * @param columns 2D array of particle data containing the target subarray (`double**`).
 * @param start Starting index of the subarray (inclusive).
 * @param end Ending index of the subarray (inclusive).
 */
static void insertion_sort_sub(double **columns, int start, int end)
{
    for (int i = start + 1; i <= end; i++)
    {
        double *temp = columns[i];
        int j = i - 1;
        while (j >= start && compare_particles(&columns[j], &temp) > 0)
        {
            columns[j + 1] = columns[j];
            j--;
        }
        columns[j + 1] = temp;
    }
}

/**
 * Parallel insertion sort implementation using chunk-based approach with overlap.
 *
 * Divides the data into num_sort_sections chunks, sorts each chunk in parallel,
 * then fixes the boundaries between chunks by re-sorting overlap regions.
 * This approach balances parallelism with the need to ensure properly sorted output.
 *
 * @param columns 2D array of particle data to be sorted
 * @param n Number of elements to sort
 */
/**
 * @brief Parallel implementation of insertion sort algorithm optimized for particle sorting.
 *
 * @details Uses multiple OpenMP threads to sort sections of particle data in parallel,
 * followed by seam-fixing operations to ensure global ordering. Includes dynamic section
 * calculation and optimized overlap sizing based on chunk characteristics.
 *
 * @param columns Column-major data array [particle][component]
 * @param n Number of particles to sort
 */
void insertion_parallel_sort(double **columns, int n)
{
    // Dynamically determine number of sections based on runtime threads and constants.
    int active_num_sort_sections;
    #ifdef _OPENMP
        int n_runtime_threads = omp_get_max_threads();
        if (n_runtime_threads <= 0) n_runtime_threads = 1;
        active_num_sort_sections = n_runtime_threads * PARALLEL_SORT_SECTIONS_PER_THREAD;
        active_num_sort_sections = n_runtime_threads * PARALLEL_SORT_SECTIONS_PER_THREAD;
        if (active_num_sort_sections <= 0) active_num_sort_sections = PARALLEL_SORT_DEFAULT_SECTIONS;
    #else
        active_num_sort_sections = 1; // Force serial behavior if OpenMP is not compiled in
    #endif

    // Ensure a reasonable number of sections
    if (active_num_sort_sections < 1) active_num_sort_sections = 1;
    if (n > 0 && active_num_sort_sections > n) active_num_sort_sections = n;
    // Optional: Add a hard cap for maximum sections if desired, e.g.:
    // if (active_num_sort_sections > 96) active_num_sort_sections = 96;

    // Fallback to serial sort for small N or if chunks would be too small
    int estimated_avg_chunk_size = (n > 0 && active_num_sort_sections > 0) ? (n / active_num_sort_sections) : n;
    if (n < PARALLEL_SORT_MIN_CHUNK_SIZE_THRESHOLD || \
        active_num_sort_sections <= 1 || \
        estimated_avg_chunk_size < PARALLEL_SORT_MIN_CHUNK_SIZE_THRESHOLD) {
        // If PARALLEL_SORT_MIN_CHUNK_SIZE_THRESHOLD is set carefully (e.g. >= 2 * PARALLEL_SORT_MIN_CORRECTNESS_OVERLAP),
        // this also helps ensure chunks are large enough for meaningful overlap.
        insertion_sort(columns, n); // Call serial insertion sort
        return;
    }

    // Determine chunk boundaries, distributing remainder
    int base_chunk_size = n / active_num_sort_sections;
    int remainder = n % active_num_sort_sections;
    int *startIdx = (int *)malloc(active_num_sort_sections * sizeof(int));
    int *endIdx = (int *)malloc(active_num_sort_sections * sizeof(int));
    if (!startIdx || !endIdx) { /* Handle error */ CLEAN_EXIT(1); }

    int offset = 0;
    for (int c = 0; c < active_num_sort_sections; c++)
    {
        int size_c = base_chunk_size + (c < remainder ? 1 : 0);
        startIdx[c] = offset;
        endIdx[c] = offset + size_c - 1;
        offset += size_c;
    }

    // Sort each chunk in parallel using insertion sort
#pragma omp parallel for schedule(dynamic)
    for (int c = 0; c < active_num_sort_sections; c++)
    {
        insertion_sort_sub(columns, startIdx[c], endIdx[c]);
    }

    // Calculate minChunkSize based on actual chunk distribution using active_num_sort_sections
    int minChunkSize = n;
    if (active_num_sort_sections > 0 && n > 0 && endIdx && startIdx) { // Check endIdx/startIdx validity
        minChunkSize = (endIdx[0] - startIdx[0] + 1);
        for (int c = 1; c < active_num_sort_sections; c++) {
            int csize = endIdx[c] - startIdx[c] + 1;
            if (csize < minChunkSize) minChunkSize = csize;
        }
    }
    if (minChunkSize <= 0 && n > 0) minChunkSize = 1; // Safety for valid n

    int overlapSize;
    if (n <= 1 || active_num_sort_sections <= 1 || minChunkSize <= 0) {
        overlapSize = 0;
    } else {
        int proportional_overlap = minChunkSize / PARALLEL_SORT_OVERLAP_DIVISOR;
        if (proportional_overlap == 0 && minChunkSize > 0) {
            proportional_overlap = 1;
        }

        // Ensure overlap is at least the minimum required for correctness,
        // but only if that minimum isn't itself making the overlap too large for the chunk.
        if (PARALLEL_SORT_MIN_CORRECTNESS_OVERLAP > 0 && proportional_overlap < PARALLEL_SORT_MIN_CORRECTNESS_OVERLAP) {
            overlapSize = PARALLEL_SORT_MIN_CORRECTNESS_OVERLAP;
        } else {
            overlapSize = proportional_overlap;
        }

        // Cap the overlap: It should not be an excessive fraction of the smallest chunk.
        // This also handles cases where MIN_CORRECTNESS_OVERLAP might be too large for a small chunk.
        int max_permissible_relative_overlap = minChunkSize / 2; // Example: Cap at 50% of chunk
        if (max_permissible_relative_overlap < 1 && minChunkSize > 0) max_permissible_relative_overlap = 1; // Ensure cap is at least 1 if chunk exists

        if (overlapSize > max_permissible_relative_overlap && minChunkSize > 1) {
            overlapSize = max_permissible_relative_overlap;
        }

        // If, after all logic, overlap is 0 but we have multiple sections and data, ensure minimal overlap.
        if (overlapSize == 0 && minChunkSize > 0 && active_num_sort_sections > 1) {
             overlapSize = 1;
        }
    }
    if (overlapSize < 0) overlapSize = 0; // Final safety check
    // Additional absolute cap based on total N, mostly for sanity with very few sections.
    if (n > 1 && overlapSize > n / 2) overlapSize = n / 2;


    // Optional debug print (controlled by -DDEBUG_SORT_PARAMS compile flag)
    #ifdef DEBUG_SORT_PARAMS
    #ifdef _OPENMP
    if (omp_get_thread_num() == 0) // Print only from one thread
    #endif
    {
        printf("[NSPHERE_IS_PARALLEL_DEBUG] N=%d, Sections=%d, MinChunkSz=%d, OverlapSize=%d (Using DIV:%d, MIN_CORRECT:%d)\n",
               n, active_num_sort_sections, minChunkSize, overlapSize,
               PARALLEL_SORT_OVERLAP_DIVISOR, PARALLEL_SORT_MIN_CORRECTNESS_OVERLAP);
        fflush(stdout);
    }
    #endif

    // Merge/fix the seams in parallel
    int nSeams = active_num_sort_sections - 1;
    if (overlapSize > 0 && nSeams > 0) {
#pragma omp parallel for schedule(dynamic)
        for (int s = 0; s < nSeams; s++)
        {
            int c_left = s;
            int c_right = s + 1;

            // Robust boundary calculations for seam sorting
            int seam_sort_start = endIdx[c_left] - overlapSize + 1;
            if (seam_sort_start < startIdx[c_left]) seam_sort_start = startIdx[c_left];

            int seam_sort_end = startIdx[c_right] + overlapSize - 1;
            if (seam_sort_end > endIdx[c_right]) seam_sort_end = endIdx[c_right];

            // Sort the combined overlap region
            if (seam_sort_start <= seam_sort_end)
            {
                insertion_sort_sub(columns, seam_sort_start, seam_sort_end);
            }
        }
    }

    free(startIdx);
    free(endIdx);
}

/**
 * Parallel quadsort implementation using chunk-based approach with overlap.
 *
 * Similar to insertion_parallel_sort but uses the quadsort algorithm for both
 * the initial chunk sorting and the seam fixing. Quadsort is generally faster
 * than insertion sort for larger datasets while maintaining stability.
 *
 * @param columns 2D array of particle data to be sorted
 * @param n Number of elements to sort
 */
/**
 * @brief Parallel implementation of quadsort algorithm optimized for particle sorting.
 *
 * @details Uses multiple OpenMP threads to sort sections of particle data in parallel,
 * followed by seam-fixing operations to ensure global ordering. Includes dynamic section
 * calculation and optimized overlap sizing based on chunk characteristics.
 *
 * @param columns Column-major data array [particle][component]
 * @param n Number of particles to sort
 */
void quadsort_parallel_sort(double **columns, int n)
{
    // Dynamically determine number of sections based on runtime threads and constants.
    int active_num_sort_sections;
    #ifdef _OPENMP
        int n_runtime_threads = omp_get_max_threads();
        if (n_runtime_threads <= 0) n_runtime_threads = 1;
        active_num_sort_sections = n_runtime_threads * PARALLEL_SORT_SECTIONS_PER_THREAD;
        active_num_sort_sections = n_runtime_threads * PARALLEL_SORT_SECTIONS_PER_THREAD;
        if (active_num_sort_sections <= 0) active_num_sort_sections = PARALLEL_SORT_DEFAULT_SECTIONS;
    #else
        active_num_sort_sections = 1; // Force serial behavior if OpenMP is not compiled in
    #endif

    // Ensure a reasonable number of sections
    if (active_num_sort_sections < 1) active_num_sort_sections = 1;
    if (n > 0 && active_num_sort_sections > n) active_num_sort_sections = n;
    // Optional: Add a hard cap for maximum sections if desired, e.g.:
    // if (active_num_sort_sections > 96) active_num_sort_sections = 96;

    // Fallback to serial sort for small N or if chunks would be too small
    int estimated_avg_chunk_size = (n > 0 && active_num_sort_sections > 0) ? (n / active_num_sort_sections) : n;
    if (n < PARALLEL_SORT_MIN_CHUNK_SIZE_THRESHOLD || \
        active_num_sort_sections <= 1 || \
        estimated_avg_chunk_size < PARALLEL_SORT_MIN_CHUNK_SIZE_THRESHOLD) {
        // If PARALLEL_SORT_MIN_CHUNK_SIZE_THRESHOLD is set carefully (e.g. >= 2 * PARALLEL_SORT_MIN_CORRECTNESS_OVERLAP),
        // this also helps ensure chunks are large enough for meaningful overlap.
        quadsort_wrapper(columns, n); // Call serial quadsort wrapper
        return;
    }

    // Determine chunk boundaries
    int base_chunk_size = n / active_num_sort_sections;
    int remainder = n % active_num_sort_sections;
    int *startIdx = (int *)malloc(active_num_sort_sections * sizeof(int));
    int *endIdx = (int *)malloc(active_num_sort_sections * sizeof(int));
    if (!startIdx || !endIdx) { /* Handle error */ CLEAN_EXIT(1); }

    int offset = 0;
    for (int c = 0; c < active_num_sort_sections; c++)
    {
        int size_c = base_chunk_size + (c < remainder ? 1 : 0);
        startIdx[c] = offset;
        endIdx[c] = offset + size_c - 1;
        offset += size_c;
    }

    // Sort each chunk in parallel using quadsort
#pragma omp parallel for schedule(dynamic)
    for (int c = 0; c < active_num_sort_sections; c++)
    {
        quadsort(&columns[startIdx[c]], endIdx[c] - startIdx[c] + 1, sizeof(double *), compare_particles);
    }

    // Calculate minChunkSize based on actual chunk distribution using active_num_sort_sections
    int minChunkSize = n;
    if (active_num_sort_sections > 0 && n > 0 && endIdx && startIdx) { // Check endIdx/startIdx validity
        minChunkSize = (endIdx[0] - startIdx[0] + 1);
        for (int c = 1; c < active_num_sort_sections; c++) {
            int csize = endIdx[c] - startIdx[c] + 1;
            if (csize < minChunkSize) minChunkSize = csize;
        }
    }
    if (minChunkSize <= 0 && n > 0) minChunkSize = 1; // Safety for valid n

    int overlapSize;
    if (n <= 1 || active_num_sort_sections <= 1 || minChunkSize <= 0) {
        overlapSize = 0;
    } else {
        int proportional_overlap = minChunkSize / PARALLEL_SORT_OVERLAP_DIVISOR;
        if (proportional_overlap == 0 && minChunkSize > 0) {
            proportional_overlap = 1;
        }

        // Ensure overlap is at least the minimum required for correctness,
        // but only if that minimum isn't itself making the overlap too large for the chunk.
        if (PARALLEL_SORT_MIN_CORRECTNESS_OVERLAP > 0 && proportional_overlap < PARALLEL_SORT_MIN_CORRECTNESS_OVERLAP) {
            overlapSize = PARALLEL_SORT_MIN_CORRECTNESS_OVERLAP;
        } else {
            overlapSize = proportional_overlap;
        }

        // Cap the overlap: It should not be an excessive fraction of the smallest chunk.
        // This also handles cases where MIN_CORRECTNESS_OVERLAP might be too large for a small chunk.
        int max_permissible_relative_overlap = minChunkSize / 2; // Example: Cap at 50% of chunk
        if (max_permissible_relative_overlap < 1 && minChunkSize > 0) max_permissible_relative_overlap = 1; // Ensure cap is at least 1 if chunk exists

        if (overlapSize > max_permissible_relative_overlap && minChunkSize > 1) {
            overlapSize = max_permissible_relative_overlap;
        }

        // If, after all logic, overlap is 0 but we have multiple sections and data, ensure minimal overlap.
        if (overlapSize == 0 && minChunkSize > 0 && active_num_sort_sections > 1) {
             overlapSize = 1;
        }
    }
    if (overlapSize < 0) overlapSize = 0; // Final safety check
    // Additional absolute cap based on total N, mostly for sanity with very few sections.
    if (n > 1 && overlapSize > n / 2) overlapSize = n / 2;


    // Optional debug print (controlled by -DDEBUG_SORT_PARAMS compile flag)
    #ifdef DEBUG_SORT_PARAMS
    #ifdef _OPENMP
    if (omp_get_thread_num() == 0) // Print only from one thread
    #endif
    {
        printf("[NSPHERE_QS_PARALLEL_DEBUG] N=%d, Sections=%d, MinChunkSz=%d, OverlapSize=%d (Using DIV:%d, MIN_CORRECT:%d)\n",
               n, active_num_sort_sections, minChunkSize, overlapSize,
               PARALLEL_SORT_OVERLAP_DIVISOR, PARALLEL_SORT_MIN_CORRECTNESS_OVERLAP);
        fflush(stdout);
    }
    #endif

    // Merge/fix the seams in parallel using quadsort
    int nSeams = active_num_sort_sections - 1;
    if (overlapSize > 0 && nSeams > 0) {
#pragma omp parallel for schedule(dynamic)
        for (int s = 0; s < nSeams; s++)
        {
            int c_left = s;
            int c_right = s + 1;

            // Robust boundary calculations for seam sorting
            int seam_sort_start = endIdx[c_left] - overlapSize + 1;
            if (seam_sort_start < startIdx[c_left]) seam_sort_start = startIdx[c_left];

            int seam_sort_end = startIdx[c_right] + overlapSize - 1;
            if (seam_sort_end > endIdx[c_right]) seam_sort_end = endIdx[c_right];

            // Sort the combined overlap region
            if (seam_sort_start <= seam_sort_end)
            {
                int seam_len = seam_sort_end - seam_sort_start + 1;
                if (seam_len > 1) {  // Only sort if there's more than one element
                    quadsort(&columns[seam_sort_start], seam_len, sizeof(double *), compare_particles);
                }
            }
        }
    }

    free(startIdx);
    free(endIdx);
}

/**
 * Validates sorting results by comparing with standard qsort.
 *
 * Creates a copy of the input array, sorts it with standard qsort,
 * then compares the results element by element to verify that the
 * sorting algorithm produced the expected ordering.
 *
 * @param columns 2D array of particle data that has been sorted
 * @param n Number of elements in the array
 * @param label Name of the sorting algorithm for diagnostic output
 */
void verify_sort_results(double **columns, int n, const char *label)
{
    // Create a temporary array of pointers to the columns
    double **tempCopy = (double **)malloc(n * sizeof(double *));
    if (!tempCopy) { fprintf(stderr, "Malloc failed in verify_sort_results\n"); return; }
    for (int i = 0; i < n; i++)
    {
        tempCopy[i] = columns[i];
    }

    // Sort the temporary pointer array using standard qsort
    stdlib_qsort_wrapper(tempCopy, n);

    // Compare the original sorted array with the qsort-ed copy
    long mismatches = 0;
    for (int i = 0; i < n; i++)
    {
        // Compare based on the actual data pointed to
        if (compare_particles(&columns[i], &tempCopy[i]) != 0)
        {
            mismatches++;
        }
    }

    if (mismatches == 0)
    {
        fprintf(stderr, "[DEBUG] SortAlg='%s': Verified => Results match standard qsort.\n", label);
    }
    else
    {
        fprintf(stderr, "[DEBUG] SortAlg='%s': *** MISMATCH *** => %ld rows differ from qsort.\n",
                label, mismatches);
    }

    free(tempCopy);
}

/**
 * Main function for sorting particle data using a specified algorithm.
 *
 * Transposes the data format from particles[component][particle] to
 * columns[particle][component], applies the specified sorting algorithm,
 * and then transposes back to the original format.
 *
 * @param particles 2D array of particle data to be sorted [component][particle]
 * @param npts Number of particles to sort
 * @param sortAlg String identifier of the sorting algorithm to use ("quadsort", "quadsort_parallel", or default to "insertion_parallel")
 */
/**
 * @brief Sorts particle data using the specified sorting algorithm.
 * @details Performs a three-phase particle sorting operation:
 *   1. Memory allocation and data transposition to column-major format
 *   2. Application of the selected sorting algorithm
 *   3. Reverse transposition of sorted data and memory cleanup
 *
 * The function uses per-call local buffer allocation for transposing data.
 *
 * @param particles 2D array of particle data to be sorted [component][particle]
 * @param npts Number of particles to sort
 * @param sortAlg Sorting algorithm to use ("quadsort", "quadsort_parallel",
 *               "insertion", or "insertion_parallel")
 */
void sort_particles_with_alg(double **particles, int npts, const char *sortAlg)
{

    /**
     * Phase 1: Memory allocation and data transposition
     * Prepares the column-major data format required for efficient sorting.
     */

    double **columns_to_sort_on; // Will point to the buffer used for sorting

    /**
     * Persistent buffer allocation strategy.
     * Uses a global buffer to reduce allocation overhead across multiple sort operations.
     */

    // Allocate or reallocate only if needed
    if (g_sort_columns_buffer == NULL || g_sort_columns_buffer_npts != npts) {
        // Free existing buffer if size has changed
        if (g_sort_columns_buffer != NULL) {
            for (int i = 0; i < g_sort_columns_buffer_npts; i++) {
                if (g_sort_columns_buffer[i]) free(g_sort_columns_buffer[i]);
            }
            free(g_sort_columns_buffer);
        }

        // Allocate new buffer with the required size
        g_sort_columns_buffer = (double **)malloc(npts * sizeof(double *));
        if (!g_sort_columns_buffer) {
            fprintf(stderr, "ERROR: Malloc failed for g_sort_columns_buffer in sort_particles_with_alg\n");
            CLEAN_EXIT(1);
        }

        // Allocate sub-arrays for each particle's components
        for (int i = 0; i < npts; i++) {
            g_sort_columns_buffer[i] = (double *)malloc(5 * sizeof(double));
            if (!g_sort_columns_buffer[i]) {
                fprintf(stderr, "ERROR: Malloc failed for g_sort_columns_buffer[%d] in sort_particles_with_alg\n", i);
                // Clean up partial allocation
                for(int k=0; k<i; ++k) free(g_sort_columns_buffer[k]);
                free(g_sort_columns_buffer);
                g_sort_columns_buffer = NULL;
                CLEAN_EXIT(1);
            }
        }
        g_sort_columns_buffer_npts = npts;
    }

    // Transpose data from row-major (particles) to column-major (g_sort_columns_buffer)
    #pragma omp parallel for
    for (int i = 0; i < npts; i++) {
        for (int j = 0; j < 5; j++) {
            g_sort_columns_buffer[i][j] = particles[j][i];
        }
    }
    columns_to_sort_on = g_sort_columns_buffer;


    /**
     * Phase 2: Apply the selected sorting algorithm
     * Uses one of several sorting algorithms based on input parameter or default.
     * Available algorithms include quadsort (sequential), quadsort_parallel,
     * insertion sort (sequential), and insertion_parallel (default).
     */

    // Select and apply the sorting algorithm
    const char *method = (sortAlg ? sortAlg : "insertion_parallel"); // Default if NULL

    if (strcmp(method, "quadsort") == 0) {
        quadsort_wrapper(columns_to_sort_on, npts);
    }
    else if (strcmp(method, "quadsort_parallel") == 0) {
        quadsort_parallel_sort(columns_to_sort_on, npts);
    }
    else if (strcmp(method, "insertion") == 0) {
        insertion_sort(columns_to_sort_on, npts);
    }
    else { // Default to parallel insertion sort
        insertion_parallel_sort(columns_to_sort_on, npts);
    }


    /**
     * Phase 3: Data transposition and memory cleanup
     * Restores sorted data to original format and performs appropriate cleanup.
     */

    // Transpose sorted data back to the original format
    #pragma omp parallel for
    for (int i = 0; i < npts; i++) {
        for (int j = 0; j < 5; j++) {
            particles[j][i] = columns_to_sort_on[i][j];
        }
    }
    // Note: The persistent buffer g_sort_columns_buffer is NOT freed here.
    // It will be reused for subsequent sort operations and freed at program exit.


}

/**
 * @brief Convenience wrapper function for sorting particles with the default algorithm.
 * @details Calls sort_particles_with_alg using the default sorting algorithm specified in g_defaultSortAlg.
 *
 * Parameters
 * ----------
 * particles : double**
 *     2D array of particle data to be sorted [component][particle]
 * npts : int
 *     Number of particles to sort
 */
void sort_particles(double **particles, int npts)
{
    sort_particles_with_alg(particles, npts, g_defaultSortAlg);
}

// =========================================================================
// SPLINE DATA SORTING UTILITIES
// =========================================================================
//
// Utility functions and structures for sorting spline data arrays.
// Provides mechanisms to sort arrays used for GSL spline creation (like radius `r`
// and potential `Psi`) while maintaining the correct correspondence between
// paired values after sorting based on one of the arrays (typically radius).

/**
 * @brief Sorts radius and potential arrays in tandem, maintaining their correspondence.
 * @details This function takes an array of radial coordinates (`rrA_spline`) and an
 *          array of corresponding potential values (`psiAarr_spline`). It sorts
 *          `rrA_spline` in ascending order and applies the identical swaps to
 *          `psiAarr_spline`, ensuring that `psiAarr_spline[i]` still corresponds to
 *          `rrA_spline[i]` after sorting. This is crucial for creating GSL splines
 *          where the x-array must be strictly monotonic and the y-array must maintain
 *          its pairing with the x-values. The arrays are assumed to have `npts + 1` elements,
 *          indexed from 0 to `npts`.
 *
 * @param rrA_spline    [in,out] Array of radial coordinates to be sorted. Modified in-place.
 * @param psiAarr_spline [in,out] Array of corresponding Psi values. Modified in-place in tandem with `rrA_spline`.
 * @param npts          The number of points, typically meaning arrays are of size `npts + 1`.
 */
void sort_rr_psi_arrays(double *rrA_spline, double *psiAarr_spline, int npts)
{
    // Allocate temporary array of pairs
    struct RrPsiPair *pairs = (struct RrPsiPair *)malloc((npts + 1) * sizeof(struct RrPsiPair));
    if (!pairs)
    {
        perror("malloc failed in sort_rr_psi_arrays");
        CLEAN_EXIT(EXIT_FAILURE);
    }

    // Populate the pairs array
    for (int i = 0; i <= npts; i++)
    {
        pairs[i].rr = rrA_spline[i];
        pairs[i].psi = psiAarr_spline[i];
    }

    // Sort the pairs based on the radius value
    qsort(pairs, npts + 1, sizeof(struct RrPsiPair), compare_by_rr);

    // Copy the sorted data back into the original arrays
    for (int i = 0; i <= npts; i++)
    {
        rrA_spline[i] = pairs[i].rr;
        psiAarr_spline[i] = pairs[i].psi;
    }

    free(pairs);
}
