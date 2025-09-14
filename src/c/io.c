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
#include "logging.h"
#include "cli.h"
#include "exit.h"
#include "io.h"
#include "density.h"
#include "density_nfw.h"
#include "utils.h"
#include <string.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <math.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/statvfs.h>
#endif

/**
 * @brief Gets the available disk space for the filesystem containing the given path.
 * @details This function uses platform-specific APIs to determine the free space
 *          available to the current user on the filesystem where `path` resides.
 *          On Windows, it uses `GetDiskFreeSpaceEx`. On POSIX-compliant systems
 *          (Linux, macOS), it uses `statvfs`.
 *
 * @param path [in] A path to a file or directory on the filesystem to check.
 *                For Windows, this can be a root directory like "C:\\".
 *                For POSIX, any path within the target filesystem, e.g., "data/".
 * @return long long Available disk space in bytes. Returns -1 on error or if the
 *                   functionality is not implemented for the current platform.
 */
long long get_available_disk_space(const char *path) {
    #ifdef _WIN32
        ULARGE_INTEGER freeBytesAvailable;
        if (GetDiskFreeSpaceEx(path, &freeBytesAvailable, NULL, NULL))
            return (long long)freeBytesAvailable.QuadPart;
    #else
        struct statvfs stat;
        if (statvfs(path, &stat) == 0)
            return (long long)stat.f_bavail * (long long)stat.f_frsize;
    #endif
    return -1;
}

/**
 * @brief Applies the global suffix to a filename.
 * @details For .dat files, inserts suffix before the extension.
 *          For other files, appends suffix to the end of filename.
 *
 * Parameters
 * ----------
 * base_filename : const char*
 *     Original filename.
 * with_suffix : int
 *     Flag indicating whether to apply the suffix (1=yes, 0=no).
 * buffer : char*
 *     Output buffer for the resulting filename.
 * bufsize : size_t
 *     Size of the output buffer.
 *
 * Returns
 * -------
 * None
 */
void get_full_filename(const char *base_filename, int with_suffix, char *buffer, size_t bufsize) {
    if (!with_suffix || g_file_suffix[0] == '\0') {
        // No suffix to apply
        strncpy(buffer, base_filename, bufsize - 1);
        buffer[bufsize - 1] = '\0';
        return;
    }

    const char *ext = strrchr(base_filename, '.');
    if (ext && strcmp(ext, ".dat") == 0) {
        // For .dat files: insert suffix before extension
        size_t basename_len = ext - base_filename;
        if (basename_len + strlen(g_file_suffix) + strlen(ext) + 1 > bufsize) {
            // Buffer too small
            strncpy(buffer, base_filename, bufsize - 1);
            buffer[bufsize - 1] = '\0';
            return;
        }

        strncpy(buffer, base_filename, basename_len);
        buffer[basename_len] = '\0';
        strcat(buffer, g_file_suffix);
        strcat(buffer, ext);
    } else {
        // For non-dat files: append suffix to filename
        if (strlen(base_filename) + strlen(g_file_suffix) + 1 > bufsize) {
            // Buffer too small
            strncpy(buffer, base_filename, bufsize - 1);
            buffer[bufsize - 1] = '\0';
            return;
        }

        strcpy(buffer, base_filename);
        strcat(buffer, g_file_suffix);
    }
}

/**
 * @brief Formats a byte count into a human-readable string with appropriate units.
 *
 * Parameters
 * ----------
 * size_in_bytes : long
 *     The size in bytes to format.
 * buffer : char*
 *     Output buffer for the formatted string.
 * buffer_size : size_t
 *     Size of the output buffer.
 *
 * Returns
 * -------
 * None
 */
void format_file_size(long size_in_bytes, char *buffer, size_t buffer_size) {
    const char *units[] = {"B", "KB", "MB", "GB", "TB"};
    int unit_index = 0;
    double size = (double)size_in_bytes;

    // Find appropriate unit
    while (size >= 1024.0 && unit_index < 4) {
        size /= 1024.0;
        unit_index++;
    }

    // Format with appropriate precision based on size
    if (unit_index == 0) // Bytes: no decimal places needed
        snprintf(buffer, buffer_size, "%ld %s", (long)size, units[unit_index]);
    else if (size >= 10) // Larger sizes: one decimal place
        snprintf(buffer, buffer_size, "%.1f %s", size, units[unit_index]);
    else // Small sizes: two decimal places
        snprintf(buffer, buffer_size, "%.2f %s", size, units[unit_index]);
}

// =========================================================================
// BINARY FILE I/O UTILITIES
// =========================================================================

/**
 * @brief Writes binary data to a file using a printf-like format string.
 * @details Parses a format string containing simplified specifiers (`%d`, `%f`, `%g`, `%e`).
 *          Writes corresponding arguments from the variadic list (`...`) as binary data.
 *          Integer types (`%d`) are written as `int`.
 *          Floating-point types (`%f`, `%g`, `%e`) are read as `double` from args but
 *          written as `float` to the file for storage efficiency.
 *
 * Parameters
 * ----------
 * fp : FILE*
 *     File pointer to write to (must be opened in binary mode).
 * format : const char*
 *     Format string with specifiers (`%d`, `%f`, `%g`, `%e`). Other characters are ignored.
 * ... :
 *     Variable arguments matching the format specifiers.
 *
 * Returns
 * -------
 * int
 *     The number of items successfully written according to the format string.
 *
 * @note Skips optional width/precision specifiers in the format string.
 * @see fscanf_bin
 */
int fprintf_bin(FILE *fp, const char *format, ...)
{
    va_list args;
    va_start(args, format);

    int count_items = 0;
    const char *p = format;

    while (*p != '\0') {
        if (*p == '%') {
            p++; // Move past '%'

            // Skip format modifiers until we find a type specifier
            while (*p && !strchr("dfge", *p) && !(*p == 'l'))
                p++;

            if (*p == 'd') {
                // Process integer format
                int val = va_arg(args, int);
                fwrite(&val, sizeof(val), 1, fp);
                count_items++;
            }
            else if (*p == 'f' || *p == 'g' || *p == 'e') {
                // Process float format (stored as 4-byte float)
                double tmp = va_arg(args, double);
                float val = (float)tmp;
                fwrite(&val, sizeof(float), 1, fp);
                count_items++;
            }
        }
        if (*p)
            p++;
    }

    va_end(args);
    return count_items;
}

/**
 * @brief Reads binary data from a file using a scanf-like format string.
 * @details Parses a format string containing simplified specifiers (`%d`, `%f`, `%g`, `%e`).
 *          Reads corresponding binary data from the file and stores it in the
 *          pointer arguments provided in the variadic list (`...`).
 *          Integer types (`%d`) are read as `int`.
 *          Floating-point types (`%f`, `%g`, `%e`) are read as `float` from the file
 *          but stored into `double*` arguments provided by the caller.
 *
 * Parameters
 * ----------
 * fp : FILE*
 *     File pointer to read from (must be opened in binary mode).
 * format : const char*
 *     Format string with specifiers (`%d`, `%f`, `%g`, `%e`). Other characters are ignored.
 * ... :
 *     Variable pointer arguments matching the format specifiers (e.g., `int*`, `double*`).
 *
 * Returns
 * -------
 * int
 *     The number of items successfully read and assigned according to the format string.
 *     Stops reading on the first failure or EOF.
 *
 * @note Skips optional width/precision specifiers in the format string.
 * @see fprintf_bin
 */
int fscanf_bin(FILE *fp, const char *format, ...)
{
    va_list args;
    va_start(args, format);

    int count_items = 0;
    const char *p = format;

    while (*p != '\0') {
        if (*p == '%') {
            p++; // Move past '%'

            // Skip format modifiers until we reach a type specifier
            while (*p && !strchr("dfge", *p) && !(*p == 'l'))
                p++;

            if (*p == 'd') {
                // Process integer format
                int *iptr = va_arg(args, int *);
                size_t nread = fread(iptr, sizeof(int), 1, fp);
                if (nread == 1)
                    count_items++;
                else
                    return count_items;
            } else if (*p == 'f' || *p == 'g' || *p == 'e') {
                // Read 4-byte float from file
                float fval;
                size_t nread = fread(&fval, sizeof(float), 1, fp);
                if (nread == 1) {
                    // Store in caller's double pointer with type conversion
                    double *dptr = va_arg(args, double *);
                    *dptr = (double)fval;
                    count_items++;
                } else
                    return count_items; // Stop on read failure
            }
        }
        if (*p)
            p++;
    }

    va_end(args);
    return count_items;
}

// =========================================================================
// RESTART AND RECOVERY MANAGEMENT
// =========================================================================

/**
 * @brief Finds the index of the last successfully processed and written snapshot in restart mode.
 * @details This function is called when `--restart` is active. It checks for the existence
 *          and basic integrity of snapshot data files (e.g., `Rank_Mass_Rad_VRad_unsorted_t%05d.dat`
 *          and `Rank_Mass_Rad_VRad_sorted_t%05d.dat`) corresponding to the timesteps listed
 *          in `snapshot_steps`.
 *          Integrity is checked by comparing file sizes against those of the first snapshot's
 *          files (snapshot_steps[0]), allowing for a small percentage tolerance (typically +/- 5%).
 *          The function aims to determine from which snapshot index `s` (in `snapshot_steps`)
 *          the post-processing (e.g., Rank file generation) should resume.
 *          It assumes `g_file_suffix` is correctly set to identify the relevant run's files.
 *
 * @param snapshot_steps [in] Array of integer timestep numbers for which snapshot data
 *                           was expected to be written (often these are the `dtwrite` intervals).
 * @param noutsnaps      [in] The total number of snapshot indices in the `snapshot_steps` array.
 *
 * @return int The index `s` into `snapshot_steps` corresponding to the *first snapshot that
 *             needs to be processed* (i.e., last valid snapshot index + 1).
 *             - Returns -1 if no valid snapshot files are found (implying processing should start from index 0).
 *             - Returns -2 if all expected snapshot files for all unique snapshot numbers exist and
 *               appear valid (implying no further snapshot processing is needed for this phase).
 *             Prints status messages to stdout and logs warnings/errors.
 */
int find_last_processed_snapshot(int *snapshot_steps, int noutsnaps) {
    // Validate input parameters before proceeding
    if (snapshot_steps == NULL || noutsnaps <= 0) {
        printf("ERROR: Invalid snapshot_steps or noutsnaps in find_last_processed_snapshot\n");
        return -1; // Start from beginning.
    }

    int last_valid_snap = -1;
    int last_valid_index = -1;
    long reference_unsorted_size = 0;
    long reference_sorted_size = 0;

    printf("Restart mode: Checking for existing data products with suffix '%s'...\n\n", g_file_suffix);

    // Initialize to track if we've checked all snapshots.
    int all_snapshots_checked = 1; // Assume all checked until we find a problem.
    int checked_count = 0;         // Count how many files we've actually checked.
    int unique_snapshots = 0;      // Count how many unique snapshot numbers we have.

    // Calculate how many unique snapshot numbers there are (might be less than noutsnaps).
    int total_writes = 0; // Determine max snapshot number.
    for (int i = 0; i < noutsnaps; i++)
        if (snapshot_steps[i] > total_writes)
            total_writes = snapshot_steps[i];
    total_writes++; // Convert from 0-based index to count.

    int *seen = (int *)calloc(total_writes, sizeof(int)); // Use calloc to initialize to 0.
    if (seen) {
        for (int i = 0; i < noutsnaps; i++) {
            int snap = snapshot_steps[i];
            if (snap >= 0 && snap < total_writes && !seen[snap]) {
                seen[snap] = 1;
                unique_snapshots++;
            }
        }
        free(seen);
    }
    else
        unique_snapshots = noutsnaps; // Fallback if memory allocation fails.

    int total_expected = unique_snapshots * 2; // Total files expected (unsorted + sorted for each unique snapshot).
    printf("Detected %d unique snapshot numbers out of %d indices. Expecting %d files total.\n\n", unique_snapshots, noutsnaps, total_expected);

    // First, see if we even have the reference snapshot file.
    if (noutsnaps > 0) {
        int snap = snapshot_steps[0];
        char fname_unsorted[256];
        char fname_sorted[256];

        char base_filename_1[256];
        snprintf(base_filename_1, sizeof(base_filename_1), "data/Rank_Mass_Rad_VRad_unsorted_t%05d.dat", snap);
        get_full_filename(base_filename_1, 1, fname_unsorted, sizeof(fname_unsorted));
        char base_filename_2[256];
        snprintf(base_filename_2, sizeof(base_filename_2), "data/Rank_Mass_Rad_VRad_sorted_t%05d.dat", snap);
        get_full_filename(base_filename_2, 1, fname_sorted, sizeof(fname_sorted));

        FILE *fun = fopen(fname_unsorted, "rb");
        FILE *fsort = fopen(fname_sorted, "rb");

        if (fun && fsort) {
            // Get file sizes.
            fseek(fun, 0, SEEK_END);
            reference_unsorted_size = ftell(fun);
            fseek(fsort, 0, SEEK_END);
            reference_sorted_size = ftell(fsort);

            // If both files have non-zero size, use as reference.
            if (reference_unsorted_size > 0 && reference_sorted_size > 0) {
                printf("Found reference file sizes from snapshot %d (unsorted: %ld bytes, sorted: %ld bytes)\n\n", snap, reference_unsorted_size,
                       reference_sorted_size);
                last_valid_snap = snap;
                last_valid_index = 0;
                checked_count += 2; // Count these two files.
            }

            fclose(fun);
            fclose(fsort);
        } else {
            // Only close if non-NULL.
            if (fun)
                fclose(fun);
            if (fsort)
                fclose(fsort);

            // If the very first files cannot be opened, processing should not continue.
            printf("Could not find initial snapshot files for index 0. Will start from the beginning.\n");
            return -1;
        }
    } else// No snapshots to check.
        return -1;

    // If reference sizes cannot be found, proper validation is not possible.
    if (reference_unsorted_size == 0 || reference_sorted_size == 0) {
        printf("Could not find valid reference file sizes. Will start from the beginning.\n\n");
        return -1;
    }

    printf("Checking all %d snapshots for completeness...\n\n", noutsnaps);

    // Now check ALL snapshots (except index 0 which we already checked).
    for (int i = 1; i < noutsnaps; i++) {
        int snap = snapshot_steps[i];

        // If this snapshot index maps to the same snapshot number as a previous index,.
        // We might be seeing duplicated snapshot numbers in the calculation.
        if (snap == snapshot_steps[0])
            printf("Warning: Duplicate snapshot number %d (index 0 and %d)\n", snap, i);

        char fname_unsorted[256];
        char fname_sorted[256];

        char base_filename_3[256];
        snprintf(base_filename_3, sizeof(base_filename_3), "data/Rank_Mass_Rad_VRad_unsorted_t%05d.dat", snap);
        get_full_filename(base_filename_3, 1, fname_unsorted, sizeof(fname_unsorted));

        char base_filename_4[256];
        snprintf(base_filename_4, sizeof(base_filename_4), "data/Rank_Mass_Rad_VRad_sorted_t%05d.dat", snap);
        get_full_filename(base_filename_4, 1, fname_sorted, sizeof(fname_sorted));

        FILE *fun = fopen(fname_unsorted, "rb");
        FILE *fsort = fopen(fname_sorted, "rb");

        if (fun && fsort) {
            // Get file sizes.
            fseek(fun, 0, SEEK_END);
            long unsorted_size = ftell(fun);
            fseek(fsort, 0, SEEK_END);
            long sorted_size = ftell(fsort);

            // Compare to reference sizes.
            // Allow for some small variation (±5%).
            double unsorted_ratio = (double)unsorted_size / reference_unsorted_size;
            double sorted_ratio = (double)sorted_size / reference_sorted_size;

            if (unsorted_size > 0 && sorted_size > 0 && in_range(unsorted_ratio,0.95,1.05) && in_range(sorted_ratio,0.95,1.05)) {
                last_valid_snap = snap;
                last_valid_index = i;
                log_message("INFO", "Verified valid files for snapshot %d (unsorted: %ld bytes, sorted: %ld bytes)", snap, unsorted_size, sorted_size);
                checked_count += 2; // Count both files as checked.
            } else {
                log_message("WARNING", "Found invalid files for snapshot %d (unsorted: %ld bytes, sorted: %ld bytes) - expected ~%ld and ~%ld bytes",
                            snap, unsorted_size, sorted_size, reference_unsorted_size, reference_sorted_size);
                fclose(fun);
                fclose(fsort);
                all_snapshots_checked = 0; // Not all snapshots are valid.
                printf("Invalid snapshot found. Will restart from this point.\n\n");
                break; // Stop at first invalid snapshot.
            }

            fclose(fun);
            fclose(fsort);
        } else {
            // File doesn't exist at all for this snapshot.
            log_message("WARNING", "Missing files for snapshot %d", snap);
            // Only close if non-NULL.
            if (fun)
                fclose(fun);
            if (fsort)
                fclose(fsort);
            all_snapshots_checked = 0; // Not all snapshots are valid.
            printf("Missing snapshot found. Will restart from this point.\n");
            break; // Stop at first missing snapshot.
        }
    }

    printf("End of file check: last_valid_snap=%d, last_valid_index=%d, noutsnaps=%d, all_snapshots_checked=%d\n\n", last_valid_snap, last_valid_index,
           noutsnaps, all_snapshots_checked);
    printf("Files checked: %d out of %d expected (unique snapshots: %d)\n\n", checked_count, total_expected, unique_snapshots);

    if (last_valid_snap == -1) {
        printf("No valid data products found. Will start from the beginning.\n\n");
        return -1;
    }

    // Check if we've verified more files than expected - this can happen if we have duplicate snapshot numbers.
    if (checked_count > total_expected)
        printf("WARNING: Checked more files (%d) than expected (%d) - likely due to duplicate snapshot numbers.\n", checked_count, total_expected);

    // Determine if we need to proceed with Rank file generation.
    if (unique_snapshots == 1 && all_snapshots_checked) {
        // Special case: Only one unique snapshot number (usually 0), and it's valid.
        printf("WARNING: Only one unique snapshot number found (%d). There's likely an issue with the calculation.\n", snapshot_steps[0]);
        printf("Only 1 Rank file (snapshot %d) exists. Starting from the beginning to create all files.\n", snapshot_steps[0]);
        return -1; // Start from beginning.
    }
    // Only say "all files exist" if:
    // 1. We have more than one unique snapshot, and.
    // 2. We've checked all expected files and found them valid.
    else if (unique_snapshots > 1 && checked_count >= total_expected && all_snapshots_checked) {
        printf("All %d data product files (for %d unique snapshots) already exist and are valid. Nothing to do.\n\n", checked_count, unique_snapshots);
        return -2; // Special code for "all done".
    } else {
        printf("Will restart processing from snapshot index %d (after snapshot %d)\n\n", last_valid_index + 1, last_valid_snap);
        return last_valid_index + 1; // Return the index of the NEXT snapshot to process.
    }
}

// =========================================================================
// INITIAL CONDITION FILE I/O FUNCTIONS
// =========================================================================

/**
 * @brief Writes particle initial conditions to a binary file.
 * @details Stores the complete initial particle state (radius, velocity, angular momentum,
 *          original index, orientation parameter mu) and the particle count (`npts`)
 *          to a binary file for later retrieval via `read_initial_conditions`.
 *          Opens the file in write binary mode ("wb"). First writes the integer `npts`,
 *          then writes the 5 double-precision values for each particle sequentially.
 *
 * Parameters
 * ----------
 * particles : double**
 *     2D array containing particle properties [component][particle_index].
 *     Expected components: 0=rad, 1=vel, 2=angmom, 3=orig_idx, 4=mu.
 * npts : int
 *     Number of particles to write.
 * filename : const char*
 *     Path to the output binary file.
 *
 * Returns
 * -------
 * None
 *
 * @note The binary file format is: `npts` (int32), followed by `npts` records,
 *       each consisting of 5 `double` values (radius, velocity, ang. mom., orig. index, mu).
 * @warning Prints an error to stderr if the file cannot be opened.
 *
 * @see read_initial_conditions
 */
void write_initial_conditions(double **particles, int npts, const char *filename) {
    printf("Saving initial conditions to %s...\n", filename);
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        fprintf(stderr, "Error: cannot open '%s' for writing initial conditions.\n", filename);
        return;
    }

   // Write the number of particles first
    fwrite(&npts, sizeof(int), 1, fp);

    for (int i = 0; i < npts; i++) {
        double r_val = particles[0][i];
        double v_val = particles[1][i];
        double ell_val = particles[2][i];
        double idx_val = particles[3][i];
        double mu_val = particles[4][i];

        fwrite(&r_val, sizeof(double), 1, fp);
        fwrite(&v_val, sizeof(double), 1, fp);
        fwrite(&ell_val, sizeof(double), 1, fp);
        fwrite(&idx_val, sizeof(double), 1, fp);
        fwrite(&mu_val, sizeof(double), 1, fp);
    }

    fclose(fp);
    log_message("INFO", "Wrote initial conditions (%d particles) to '%s'", npts, filename);
}

/**
 * @brief Reads particle initial conditions from a binary file.
 * @details Loads the complete particle state (radius, velocity, angular momentum,
 *          original index, orientation parameter mu) from a binary file previously
 *          created by `write_initial_conditions`. Verifies that the number of
 *          particles read from the file matches the expected count `npts`.
 *          Opens the file in read binary mode ("rb").
 *
 * Parameters
 * ----------
 * particles : double**
 *     2D array to store the loaded particle properties [component][particle_index].
 *     Must be pre-allocated with dimensions [5][npts].
 * npts : int
 *     Expected number of particles to read.
 * filename : const char*
 *     Path to the input binary file.
 *
 * Returns
 * -------
 * None (populates the `particles` array).
 *
 * @note See `write_initial_conditions` for file format details.
 * @warning Prints an error to stderr if the file cannot be opened, if the particle
 *          count doesn't match `npts` (returns early), or if a read error occurs.
 *
 * @see write_initial_conditions
 */
void read_initial_conditions(double **particles, int npts, const char *filename) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        fprintf(stderr, "Error: cannot open '%s' for reading initial conditions.\n", filename);
        return;
    }

    int file_npts;
    if (fread(&file_npts, sizeof(int), 1, fp) != 1) {
        fprintf(stderr, "Error: failed to read npts from '%s'.\n", filename);
        fclose(fp);
        return;
    }
    if (file_npts != npts) {
        fprintf(stderr, "Warning: file npts=%d doesn't match current npts=%d.\n", file_npts, npts);
        fclose(fp);
        return;
    }

    for (int i = 0; i < npts; i++) {
        double r_val, v_val, ell_val, idx_val, mu_val;
        if (fread(&r_val, sizeof(double), 1, fp) != 1 || fread(&v_val, sizeof(double), 1, fp) != 1 || fread(&ell_val, sizeof(double), 1, fp) != 1 || fread(&idx_val, sizeof(double), 1, fp) != 1 || fread(&mu_val, sizeof(double), 1, fp) != 1) {
            fprintf(stderr, "Error: partial read at i=%d in '%s'.\n", i, filename);
            fclose(fp);
            return;
        }

        particles[0][i] = r_val;
        particles[1][i] = v_val;
        particles[2][i] = ell_val;
        particles[3][i] = idx_val;
        particles[4][i] = mu_val;
    }

    fclose(fp);
    log_message("INFO", "Read initial conditions (%d particles) from '%s'", npts, filename);
}

// =========================================================================
// PARTICLE I/O AND OPERATIONS
// =========================================================================

/**
 * @brief Appends a block of full particle data to the specified output file.
 * @details This function writes a chunk of particle data, corresponding to `block_size`
 *          timesteps, to the given binary file. The data for each particle (rank, radius,
 *          radial velocity, angular momentum) is written sequentially for each timestep
 *          within the block (step-major order). This means all particle data for step `s`
 *          is contiguous, followed by all data for step `s+1`, etc., within the block.
 *          The file is opened in append binary mode ("ab").
 *          This is used for creating the `all_particle_data.dat` file which stores
 *          the complete evolution history when `g_doAllParticleData` is enabled.
 *
 * @param filename   [in] Path to the output binary file (e.g., "data/all_particle_data<suffix>.dat").
 * @param npts       [in] Number of particles.
 * @param block_size [in] Number of timesteps of data contained in the provided `_block` arrays.
 * @param L_block    [in] Pointer to the block of angular momentum data (float array).
 *                        Assumed to be `[step_in_block * npts + particle_orig_id]`.
 * @param Rank_block [in] Pointer to the block of particle rank data (int array).
 *                        Assumed to be `[step_in_block * npts + particle_orig_id]`.
 * @param R_block    [in] Pointer to the block of radial position data (float array).
 *                        Assumed to be `[step_in_block * npts + particle_orig_id]`.
 * @param Vrad_block [in] Pointer to the block of radial velocity data (float array).
 *                        Assumed to be `[step_in_block * npts + particle_orig_id]`.
 * @note Exits via `CLEAN_EXIT(1)` if the file cannot be opened for appending.
 */
void append_all_particle_data_chunk_to_file(const char *filename, int npts, int block_size, float *L_block, int *Rank_block, float *R_block,
                                            float *Vrad_block) {
    FILE *f = fopen(filename, "ab");
    if (!f)
        raise_error("Error: cannot open %s for appending all_particle_data\n", filename);

    // Write particle data in step-major order
    for (int step = 0; step < block_size; step++)
        for (int i = 0; i < npts; i++) {
            int rankval = Rank_block[step * npts + i];
            float rval = R_block[step * npts + i];
            float vval = Vrad_block[step * npts + i];
            float lval = L_block[step * npts + i];

            // Write particle data in fixed order
            fwrite(&rankval, sizeof(int), 1, f);
            fwrite(&rval, sizeof(float), 1, f);
            fwrite(&vval, sizeof(float), 1, f);
            fwrite(&lval, sizeof(float), 1, f);
        }
    fclose(f);
}

/**
 * @brief Retrieves particle data for a specific snapshot from the `all_particle_data.dat` binary file.
 * @details This function reads the data (rank, radius, radial velocity, angular momentum)
 *          for all `npts` particles corresponding to a single snapshot number (`snap`)
 *          from the specified binary file. The file is expected to be in step-major order,
 *          where each record per particle consists of an int (rank) and three floats (R, Vrad, L).
 *          It calculates the correct file offset to seek to the desired snapshot.
 *          To ensure thread safety when called in parallel (e.g., during post-processing
 *          of snapshots), file I/O (seeking and reading) is performed within an
 *          OpenMP critical section named `file_access`. Temporary local buffers are used
 *          for reading, and data is then copied to the caller-provided output arrays.
 *
 * @param filename    [in] Path to the binary data file (e.g., "data/all_particle_data<suffix>.dat").
 * @param snap        [in] The snapshot number (0-indexed, corresponding to write events) to retrieve.
 * @param npts        [in] Number of particles per snapshot.
 * @param block_size  [in] The number of snapshots that were written per block to the file by
 *                       `append_all_particle_data_chunk_to_file`. Used for calculating seek offset.
 * @param L_out       [out] Pointer to an array (size `npts`) to store the retrieved angular momentum values.
 * @param Rank_out    [out] Pointer to an array (size `npts`) to store the retrieved particle ranks.
 * @param R_out       [out] Pointer to an array (size `npts`) to store the retrieved radial positions.
 * @param Vrad_out    [out] Pointer to an array (size `npts`) to store the retrieved radial velocities.
 * @note Exits via `CLEAN_EXIT(1)` on memory allocation failure, file open failure, fseek failure, or unexpected EOF.
 */
void retrieve_all_particle_snapshot(const char *filename, int snap, int npts, int block_size, float *L_out, int *Rank_out, float *R_out, float *Vrad_out) {
    // Allocate local (thread-private) arrays
    float *tmpL = (float *)malloc(npts * sizeof(float));
    int *tmpRank = (int *)malloc(npts * sizeof(int));
    float *tmpR = (float *)malloc(npts * sizeof(float));
    float *tmpV = (float *)malloc(npts * sizeof(float));

    if (!tmpL || !tmpRank || !tmpR || !tmpV)
        raise_error("Error: out of memory in retrieve_all_particle_snapshot!\n");
    // Status messages are handled in the ordered section of the parallel loop.

    // Read from file in a critical section
    #pragma omp critical(file_access)
    {
        FILE *f = fopen(filename, "rb");
        if (!f)
            raise_error("Error: cannot open %s for reading\n", filename);

        // Compute offset in file.
        int block_number = snap / block_size;
        int index_in_block = snap % block_size;

        long long step_data_size = (long long)npts * 16; // 16 bytes per record.
        long long block_data_size = (long long)block_size * step_data_size;
        long long offset = block_data_size * block_number + step_data_size * index_in_block;

        if (fseek(f, offset, SEEK_SET) != 0) {
            fclose(f);
            raise_error("Error: fseek failed for snap=%d\n", snap);
            CLEAN_EXIT(1);
        }

        // Read npts records into local buffers.
        for (int i = 0; i < npts; i++) {
            int rankval;
            float rval = 0.0, vval = 0.0, lval = 0.0;  // Just to silence the compiler warnings

            if (fread(&rankval, sizeof(int), 1, f) != 1 || fread(&rval, sizeof(float), 1, f) != 1 || fread(&vval, sizeof(float), 1, f) != 1 || fread(&lval, sizeof(float), 1, f) != 1) {
                fclose(f);
                raise_error("Error: unexpected EOF while reading snap=%d (particle %d)\n", snap, i);
            }
            tmpRank[i] = rankval;
            tmpR[i] = rval;
            tmpV[i] = vval;
            tmpL[i] = lval;
        }
        fclose(f);
    } // End critical(file_access).

    // Copy from local buffers to output arrays
    memcpy(Rank_out, tmpRank, npts * sizeof(int));
    memcpy(R_out, tmpR, npts * sizeof(float));
    memcpy(Vrad_out, tmpV, npts * sizeof(float));
    memcpy(L_out, tmpL, npts * sizeof(float));

    // Free the temporary local arrays
    free(tmpL);
    free(tmpRank);
    free(tmpR);
    free(tmpV);
}

/**
 * @brief Build the filename tag string based on options for display purposes
 */
void compile_filename_tag() {
    if (custom_tag[0] != '\0')
        strcat(filename_tag, custom_tag);

    if (include_method_in_suffix) {
        if (filename_tag[0] != '\0')
            strcat(filename_tag, "_");
        strcat(filename_tag, method_filename);
    }
}

/**
 * @brief Create the 'init' directory if it doesn't exist.
 */
void mkdir_init(){
    /** @brief Create the output 'data' directory if it doesn't exist. */
    struct stat st = {0};
    if (stat("data", &st) == -1)
        mkdir("data", 0755); // POSIX standard, works on most systems including MinGW/Cygwin

    struct stat st_init = {0};
    if (stat("init", &st_init) == -1) {
        #if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
            if (mkdir("init") != 0) // Decide if this is fatal - perhaps not if only writing
                 perror("Error creating init directory");
            else
                 log_message("INFO", "Created init/ directory.");
        #else
            if (mkdir("init", 0755) != 0) // Decide if this is fatal
                 perror("Error creating init directory"); // POSIX standard
            else
                 log_message("INFO", "Created init/ directory.");
        #endif
    }
}

/**
 * @brief Write current run parameters to `data/lastparams<suffix>.dat` and create a standard link `data/lastparams.dat`.
 */
void write_to_lastparams() {
    char file_tag[512] = "";
    if (custom_tag[0] != '\0') {
        strcat(file_tag, custom_tag);
        if (include_method_in_suffix) {
            strcat(file_tag, "_");
            strcat(file_tag, method_filename);
        }
    }
    else if (include_method_in_suffix)
        strcat(file_tag, method_filename);

    /** @note Create filename with suffix for the specific run parameters file. */
    char filename[512]; // Holds suffixed filename, e.g., data/lastparams_run1_100k_10k_5.dat
    get_full_filename("data/lastparams.dat", 1, filename, sizeof(filename));
    printf("Saving parameters: %s\n", filename);

    FILE *fp_params = fopen(filename, "w"); // Text mode for regular fprintf
    if (!fp_params)
        raise_error("Error: cannot open %s\n", filename);

    fprintf(fp_params, "%d %d %d %s\n", npts, Ntimes, tfinal_factor, file_tag);
    fclose(fp_params);

    /** @note Create a standard-named link `data/lastparams.dat` pointing to the
             suffixed version for compatibility with scripts. Uses copy on Windows. */
    char linkname[512] = "data/lastparams.dat"; // Standard name

    /* Platform detection using standard predefined macros */
    #if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
        /** @note Windows: Copy file content as symlinks can be unreliable or require special privileges. */
        // Windows or Windows-like environment: create a direct file copy.
            // Ensure compatibility with various Windows environments.

        // Use lower-level file operations instead of system commands for better compatibility.
        FILE *source, *dest;
        source = fopen(filename, "rb");
        if (!source)
            printf("Warning: Failed to open source file %s for copying\n", filename);
        else {
            dest = fopen(linkname, "wb");
            if (!dest) {
                printf("Warning: Failed to create destination file %s\n", linkname);
                fclose(source);
            } else {
                // Copy file content.
                char buffer[4096];
                size_t bytes_read;

                while ((bytes_read = fread(buffer, 1, sizeof(buffer), source)) > 0)
                    fwrite(buffer, 1, bytes_read, dest);

                fclose(dest);
                fclose(source);

                // Extract the basename for display purposes
                const char *basename = strrchr(filename, '/');
                basename = basename ? basename + 1 : filename; // Skip the '/' or use full name if no '/'

                printf("Created link: %s -> %s\n\n", basename, linkname);
            }
        }
    #else
        /** @note Unix: Create symbolic link from basename(filename) to 'linkname', fallback to copy. */
        // Unix-like systems (Linux, macOS, etc.) and fallback for other platforms: use symbolic links.
        char command[1024];

        // First, remove any existing link or file.
        snprintf(command, sizeof(command), "rm -f \"%s\" 2>/dev/null", linkname);
        system(command);

        // Then create the symbolic link - use the basename of the file, not the full path
        // Extract the basename from filename
        const char *basename = strrchr(filename, '/');
        basename = basename ? basename + 1 : filename; // Skip the '/' or use full name if no '/'

        snprintf(command, sizeof(command), "ln -s \"%s\" \"%s\"", basename, linkname);
        if (system(command) != 0) {
            // If symbolic link fails, fall back to copying the file.
            snprintf(command, sizeof(command), "cp \"%s\" \"%s\"", filename, linkname);
            if (system(command) != 0)
                printf("Warning: Failed to create link or copy %s to %s\n", filename, linkname);
            else
                printf("Created link: %s -> %s\n\n", filename, linkname);
        } else
            printf("Created link: %s -> %s\n\n", filename, linkname);
    #endif
}

/**
 * @brief Write trajectory data (radius, energy, angular momentum) for selected low-L particles.
 * @details Outputs the time evolution of radius, relative energy, and angular momentum
 *          for `nlowest` particles selected based on their initial angular momentum
 *          (either lowest absolute L or closest to a reference L, per `use_closest_to_Lcompare`).
 *          The specific particles are stored in the `chosen` array (by their original IDs).
 *          Written to `lowest_l_trajectories.dat`.
 */
void write_low_l_particles(double dt, int nlowest, double **lowestL_r, double **lowestL_E, double **lowestL_L){
    // Write trajectories for selected lowest-L particles if simulation was run
    if (skip_file_writes)
        return;
    char full_filename[256];
    get_full_filename("data/lowest_l_trajectories.dat", 1, full_filename, sizeof(full_filename));
    FILE *fp_lowest = fopen(full_filename, "wb"); // Binary mode for fprintf_bin
    if (!fp_lowest)
        raise_error("Error: cannot open data/lowest_l_trajectories.dat\n");

    // Write data for Ntimes steps:
    for (int step = 0; step < Ntimes; step++) {
        double tval = step * dt;
        fprintf_bin(fp_lowest, "%f", tval);
        for (int p = 0; p < nlowest; p++) {
            double rr = lowestL_r[p][step];
            double Ecur = lowestL_E[p][step];
            double lcur = lowestL_L[p][step];
            fprintf_bin(fp_lowest, " %f %f %f", rr, Ecur, lcur);
        }
        fprintf_bin(fp_lowest, "\n");
    }
    fclose(fp_lowest);
}

/**
 * @brief OUTPUT FILE INITIALIZATION block (`all_particle_data.dat`).
 * @details Creates (or overwrites) an empty binary file `all_particle_data<suffix>.dat`
 *          if saving all particle data (`g_doAllParticleData` is true) AND the simulation
 *          is *not* being skipped (`skip_simulation` is false). This file will be appended to
 *          incrementally during the simulation timestepping loop via block writes.
 * @see append_all_particle_data_chunk_to_file
 * @see apd_filename
 * @see g_doAllParticleData
 * @see skip_simulation
 */
void initialize_output_file(int noutsnaps) {
    if (!g_doAllParticleData || skip_simulation)
        return;

    // Create/truncate the output file in binary write mode.
    FILE *fapd = fopen(apd_filename, "wb");
    if (!fapd)
        raise_error("Error: cannot create all_particle_data output file %s\n", apd_filename);
    fclose(fapd); // Close immediately, file is now ready for appending.
    printf("Initialized empty file for all particle data: %s\n", apd_filename);

    // Calculate and display expected file size
    long long expected_size = (long long)total_writes() * (long long)npts * 16LL; // 16 bytes per particle record

    printf("All particle data file requires: ");
    log_disk_space(expected_size);

    // Calculate and display expected snapshot file sizes
    // Each snapshot has 2 files: unsorted (28 bytes/particle) and sorted (32 bytes/particle)
    long long snapshot_size = (long long)npts * (28LL + 32LL); // Total per snapshot pair
    long long total_snapshot_size = snapshot_size * (long long)noutsnaps;

    printf("%d time snapshot files will require: ", noutsnaps);
    log_disk_space(total_snapshot_size);

    // Calculate and display total disk space
    long long total_disk_space = expected_size + total_snapshot_size;

    printf("Total disk space required: ");
    log_disk_space(total_disk_space);

    // Check available disk space
    long long available_space = get_available_disk_space("data/");
    if (available_space > 0) {
        printf("Available disk space: ");
        log_disk_space(available_space);

        // Check if we're within 5% of total available or insufficient
        double usage_after = (double)(available_space - total_disk_space) / (double)available_space;

        if (available_space < total_disk_space)
            raise_insufficient_memory(total_disk_space, available_space); // Insufficient space
        else if (usage_after < 0.05) { // Within 5% of capacity after simulation
            warn_low_memory(total_disk_space, available_space, usage_after);
            prompt_quit("Continue", "Aborting simulation.\n");
        }
    } else {
        fprintf(stderr, "Warning: Could not determine available disk space.\n");
        prompt_quit("Continue without disk space check", "Aborting simulation.\n");
    }

    printf("\n");

    /** @brief Display initial simulation progress. */
    printf("0%% complete, timestep 0/%d, time=0.0000 Myr, elapsed=0.00 s\n", Ntimes);
}

/** @brief Initialize the global file suffix string `g_file_suffix` based on cmd line args. */
void fill_suffix_tags() {
    g_file_suffix[0] = '\0';

    /** @brief Add custom tag to suffix if provided via `--tag`. */
    if (custom_tag[0] != '\0')
        snprintf(g_file_suffix, sizeof(g_file_suffix), "_%s", custom_tag);

    /** @brief Add method/parameter tag to suffix. */
    char temp[256];
    if (include_method_in_suffix)
        snprintf(temp, sizeof(temp), "_%s_%d_%d_%d", method_filename, npts, Ntimes, tfinal_factor);
    else
        snprintf(temp, sizeof(temp), "_%d_%d_%d", npts, Ntimes, tfinal_factor);

    /** @brief Append the parameter tag to the global suffix. */
    strcat(g_file_suffix, temp); // Append temp to g_file_suffix
}

int validate_snapshot_memory_allocation(double *r_grid, double *log_r_grid, double *mass_grid, double *density_grid, double *density_sorted,
                                        double *R_decimated, double *Mass_decimated, double *R_filtered, double *Mass_filtered, int r_violations,
                                        int snap) {
    if (!r_grid || !log_r_grid || !mass_grid || !density_grid || !density_sorted) {
        log_message("ERROR", "Thread %d: Failed to allocate grid arrays for snapshot %d", omp_get_thread_num(), snap);
        // Free previously allocated resources
        free(R_decimated); free(Mass_decimated);
        if (r_violations > 0) {
            free(R_filtered);
            free(Mass_filtered);
        }
        free(r_grid);
        free(log_r_grid);
        free(mass_grid);
        free(density_grid);
        free(density_sorted);
        density_sorted = NULL;
        return 1;
    }
    return 0;
}

/**
 * @def Write theoretical potential profile
 */
void write_potential_profile(FILE *fp) {
    get_full_filename("data/Psiprofile.dat", 1, full_filename, sizeof(full_filename));
    fp = fopen(full_filename, "wb");
    if (fp) {
        for (double r_plot = 0.0; r_plot < radius[num_points - 1]; r_plot += (radius[num_points - 1] / 900.0))
            if (r_plot >= radius[0])
                 fprintf_bin(fp, "%f %f\n", r_plot, evaluatespline(splinePsi, Psiinterp, r_plot));
         if (num_points > 0)
             fprintf_bin(fp, "%f %f\n", radius[num_points-1], evaluatespline(splinePsi, Psiinterp, radius[num_points-1]));
        fclose(fp);
    } else
        log_message("ERROR", "Failed to open %s for final Psi profile", full_filename);
}

/**
 * @def Write theoretical mass profile
 */
void write_mass_profile(FILE *fp) {
    get_full_filename("data/massprofile.dat", 1, full_filename, sizeof(full_filename));
    fp = fopen(full_filename, "wb");
    if (fp) {
        for (double r_plot = 0.0; r_plot < radius[num_points - 1]; r_plot += (radius[num_points - 1] / 900.0))
            if (r_plot >= radius[0])
                fprintf_bin(fp, "%f %f\n", r_plot, gsl_spline_eval(splinemass, r_plot, enclosedmass));
        if (num_points > 0)
             fprintf_bin(fp, "%f %f\n", radius[num_points-1], gsl_spline_eval(splinemass, radius[num_points-1], enclosedmass));
        fclose(fp);
    } else
        log_message("ERROR", "Failed to open %s for final mass profile", full_filename);
}

/**
 * @def Write theoretical density profile
 */
void write_density_profile(FILE *fp) { //TODO - clean this up to make it generic and not nfw specific
    get_full_filename("data/density_profile.dat", 1, full_filename, sizeof(full_filename));
    fp = fopen(full_filename, "wb");
    if (fp) {
        if (g_use_nfw_profile) {
            double nt_nfw_scaler_final = g_nfw_profile_halo_mass / (4.0 * PI * normalization);
            for (int i = 0; i < num_points; i++) {
                double rr = radius[i];
                double rs_k = rr / g_nfw_profile_rc;
                double term_s_k = rs_k + 0.01;
                if (term_s_k <= 1e-9) term_s_k = 1e-9;
                double term_n_k = (1.0 + rs_k) * (1.0 + rs_k);
                double term_c_base_k = rs_k / g_nfw_profile_falloff_factor;
                double term_c_k = 1.0 + pow(term_c_base_k, 10.0);
                double rho_shape_k = (term_s_k < 1e-9 || term_n_k < 1e-9 || term_c_k < 1e-9) ? 0.0 : (1.0 / (term_s_k * term_n_k * term_c_k));
                if (rr < 1e-6 && term_s_k < 1e-3 && rho_shape_k == 0.0)
                   rho_shape_k = 1.0 / (term_s_k * term_n_k * term_c_k);
                double rho_r_k = nt_nfw_scaler_final * rho_shape_k;
                fprintf_bin(fp, "%f %f\n", rr, rho_r_k);
            }
        } else
            for (int i = 0; i < num_points; i++) {
                double rr = radius[i];
                double rho_r = g_cored_profile_halo_mass / normalization * (1.0 / cube(1.0 + sqr(rr / g_cored_profile_rc)));
                fprintf_bin(fp, "%f %f\n", rr, rho_r);
        }
        fclose(fp);
    } else
        log_message("ERROR", "Failed to open %s for final density profile", full_filename);
}

/**
 * @def Write theoretical dPsi/dr profile
 */
void write_dPsidr_profile(FILE *fp) {
    get_full_filename("data/dpsi_dr.dat", 1, full_filename, sizeof(full_filename));
    fp = fopen(full_filename, "wb");
    if (fp) {
        for (int i = 0; i < num_points; i++) {
            double rr = radius[i];
            if (rr > 0.0) {
                double Menc = gsl_spline_eval(splinemass, rr, enclosedmass);
                double dpsidr = -(G_CONST * Menc) / sqr(rr);
                fprintf_bin(fp, "%f %f\n", rr, dpsidr);
            }
        }
        fclose(fp);
    } else
        log_message("ERROR", "Failed to open %s for final dpsi/dr profile", full_filename);
}

/**
 * @def Write theoretical drho/dPsi profile
 */
void write_drhodPsi_profile(FILE *fp) {
    get_full_filename("data/drho_dpsi.dat", 1, full_filename, sizeof(full_filename));
    fp = fopen(full_filename, "wb");
    if (fp) {
        for (int i = 1; i < num_points - 1; i++) {
            double Psi_val = evaluatespline(splinePsi, Psiinterp, radius[i]);
            double drhodPsi = calculate_density_drhodPsi(i);
            if (isnan(drhodPsi))
                continue;
            fprintf_bin(fp, "%f %f\n", Psi_val, drhodPsi);
        }
        fclose(fp);
    } else
        log_message("ERROR", "Failed to open %s for final drho/dpsi profile", full_filename);
}

/**
 * @def Write theoretical f(E) profile
 */
void write_f_of_E_profile(FILE *fp) {
    get_full_filename("data/f_of_E.dat", 1, full_filename, sizeof(full_filename));
    fp = fopen(full_filename, "wb");
    if (fp) {
        for (int i = 0; i <= num_points; i++) {
            double E = Evalues[i];
            double deriv = 0.0;
            if (i > 0 && i < num_points + 1) {
                if (i > 0 && i < num_points)
                    deriv = (innerintegrandvalues[i + 1] - innerintegrandvalues[i - 1]) / (Evalues[i + 1] - Evalues[i - 1]);
                else if (i == 0)
                    deriv = (innerintegrandvalues[i + 1] - innerintegrandvalues[i]) / (Evalues[i + 1] - Evalues[i]);
                else if (i == num_points)
                    deriv = (innerintegrandvalues[i] - innerintegrandvalues[i - 1]) / (Evalues[i] - Evalues[i - 1]);
            }
            double fE = fabs(deriv) / (sqrt(8.0) * sqr(PI));
            if (E == 0.0 || !isfinite(fE))
                fE = 0.0;
            fprintf_bin(fp, "%f %f\n", E, fE);
        }
        fclose(fp);
    } else
        log_message("ERROR", "Failed to open %s for final f(E) profile", full_filename);
}

/**
 * @def Write distribution function at a fixed radius if simulation was run
 */
void write_distribution_after_simulation(FILE *fp) {
    if (skip_file_writes)
        return;
    get_full_filename("data/df_fixed_radius.dat", 1, full_filename, sizeof(full_filename));
    fp = fopen(full_filename, "wb");
    if (fp) {
        double r_F = 2.0 * density_rc();
        double Psi_rf = evaluatespline(splinePsi, Psiinterp, r_F);
        Psi_rf *= VEL_CONV_SQ;
        double Psimin_test = VEL_CONV_SQ * Psimin; // Convert to (km/s)² for velocity calculation

        int vsteps = 10000;
        int reduce_vsteps = 300;
        for (int vv = 0; vv <= vsteps - reduce_vsteps; vv++) {
            double sqrt_arg_v = Psi_rf - Psimin_test;
            if (sqrt_arg_v < 0)
                sqrt_arg_v = 0;
            double vtest = (double)vv * (sqrt(2.0 * sqrt_arg_v) / (vsteps));
            double Etest = Psi_rf - 0.5 * vtest * vtest;
            Etest = Etest / VEL_CONV_SQ; // Convert back to code units for bounds check
            double fEval = 0.0;
            if (Etest >= Psimin && Etest <= Psimax) { // Bounds check in code units
                double derivative;
                int status = gsl_interp_eval_deriv_e(g_main_fofEinterp, Evalues, innerintegrandvalues, Etest, g_main_fofEacc, &derivative);
                if (status == GSL_SUCCESS)
                    fEval = derivative / (sqrt(8.0) * sqr(PI)) * sqr(vtest) * sqr(r_F);
            }
            if (!isfinite(fEval))
                fEval = 0.0;
            fprintf_bin(fp, "%f %f\n", vtest, fEval);
        }
        fclose(fp);
    } else
        log_message("ERROR", "Failed to open %s for final df_fixed_radius", full_filename);
}

/**
 * @brief Write final theoretical profile characteristics to .dat files.
 * @details This block outputs several files (massprofile, Psiprofile, density_profile,
 *          dpsi_dr, drho_dpsi, f_of_E, df_fixed_radius) using the splines
 *          (e.g., splinemass, splinePsi, g_main_fofEinterp) and parameters
 *          (e.g., num_points, radius, normalization, g_nfw_profile_rc, etc.)
 *          that were established during the main density initial condition generation phase.
 *          Analytical formulas for the density and its derivatives are used where appropriate.
 */
void write_full_density_output(FILE *fp) {
    if (skip_file_writes)
        return;
    log_message("INFO", "Writing density theoretical profiles to final .dat files...");
    write_mass_profile(fp);
    write_potential_profile(fp);
    write_density_profile(fp);
    write_dPsidr_profile(fp);
    write_drhodPsi_profile(fp);
    write_f_of_E_profile(fp);
    write_distribution_after_simulation(fp);
}
