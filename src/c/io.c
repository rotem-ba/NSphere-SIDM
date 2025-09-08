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
 #include "exit.h"
 #include <string.h>
 #include <stdarg.h>
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
         if (GetDiskFreeSpaceEx(path, &freeBytesAvailable, NULL, NULL)) {
             return (long long)freeBytesAvailable.QuadPart;
         }
     #else
         struct statvfs stat;
         if (statvfs(path, &stat) == 0) {
             return (long long)stat.f_bavail * (long long)stat.f_frsize;
         }
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
 void get_suffixed_filename(const char *base_filename, int with_suffix, char *buffer, size_t bufsize)
 {
     if (!with_suffix || g_file_suffix[0] == '\0')
     {
         // No suffix to apply
         strncpy(buffer, base_filename, bufsize - 1);
         buffer[bufsize - 1] = '\0';
         return;
     }

     const char *ext = strrchr(base_filename, '.');
     if (ext && strcmp(ext, ".dat") == 0)
     {
         // For .dat files: insert suffix before extension
         size_t basename_len = ext - base_filename;
         if (basename_len + strlen(g_file_suffix) + strlen(ext) + 1 > bufsize)
         {
             // Buffer too small
             strncpy(buffer, base_filename, bufsize - 1);
             buffer[bufsize - 1] = '\0';
             return;
         }

         strncpy(buffer, base_filename, basename_len);
         buffer[basename_len] = '\0';
         strcat(buffer, g_file_suffix);
         strcat(buffer, ext);
     }
     else
     {
         // For non-dat files: append suffix to filename
         if (strlen(base_filename) + strlen(g_file_suffix) + 1 > bufsize)
         {
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
 void format_file_size(long size_in_bytes, char *buffer, size_t buffer_size)
 {
     const char *units[] = {"B", "KB", "MB", "GB", "TB"};
     int unit_index = 0;
     double size = (double)size_in_bytes;

     // Find appropriate unit
     while (size >= 1024.0 && unit_index < 4)
     {
         size /= 1024.0;
         unit_index++;
     }

     // Format with appropriate precision based on size
     if (unit_index == 0)
     {
         // Bytes: no decimal places needed
         snprintf(buffer, buffer_size, "%ld %s", (long)size, units[unit_index]);
     }
     else if (size >= 10)
     {
         // Larger sizes: one decimal place
         snprintf(buffer, buffer_size, "%.1f %s", size, units[unit_index]);
     }
     else
     {
         // Small sizes: two decimal places
         snprintf(buffer, buffer_size, "%.2f %s", size, units[unit_index]);
     }
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

     while (*p != '\0')
     {
         if (*p == '%')
         {
             p++; // Move past '%'

             // Skip format modifiers until we find a type specifier
             while (*p && !strchr("dfge", *p) && !(*p == 'l'))
             {
                 p++;
             }

             if (*p == 'd')
             {
                 // Process integer format
                 int val = va_arg(args, int);
                 fwrite(&val, sizeof(val), 1, fp);
                 count_items++;
             }
             else if (*p == 'f' || *p == 'g' || *p == 'e')
             {
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

     while (*p != '\0')
     {
         if (*p == '%')
         {
             p++; // Move past '%'

             // Skip format modifiers until we reach a type specifier
             while (*p && !strchr("dfge", *p) && !(*p == 'l'))
             {
                 p++;
             }

             if (*p == 'd')
             {
                 // Process integer format
                 int *iptr = va_arg(args, int *);
                 size_t nread = fread(iptr, sizeof(int), 1, fp);
                 if (nread == 1)
                     count_items++;
                 else
                     return count_items;
             }
             else if (*p == 'f' || *p == 'g' || *p == 'e')
             {
                 // Read 4-byte float from file
                 float fval;
                 size_t nread = fread(&fval, sizeof(float), 1, fp);
                 if (nread == 1)
                 {
                     // Store in caller's double pointer with type conversion
                     double *dptr = va_arg(args, double *);
                     *dptr = (double)fval;
                     count_items++;
                 }
                 else
                 {
                     return count_items; // Stop on read failure
                 }
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
 int find_last_processed_snapshot(int *snapshot_steps, int noutsnaps)
 {
     // Validate input parameters before proceeding
     if (snapshot_steps == NULL || noutsnaps <= 0)
     {
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
     {
         if (snapshot_steps[i] > total_writes)
         {
             total_writes = snapshot_steps[i];
         }
     }
     total_writes++; // Convert from 0-based index to count.

     {
         int *seen = (int *)calloc(total_writes, sizeof(int)); // Use calloc to initialize to 0.
         if (seen)
         {
             for (int i = 0; i < noutsnaps; i++)
             {
                 int snap = snapshot_steps[i];
                 if (snap >= 0 && snap < total_writes && !seen[snap])
                 {
                     seen[snap] = 1;
                     unique_snapshots++;
                 }
             }
             free(seen);
         }
         else
         {
             unique_snapshots = noutsnaps; // Fallback if memory allocation fails.
         }
     }

     int total_expected = unique_snapshots * 2; // Total files expected (unsorted + sorted for each unique snapshot).
     printf("Detected %d unique snapshot numbers out of %d indices. Expecting %d files total.\n\n",
            unique_snapshots, noutsnaps, total_expected);

     // First, see if we even have the reference snapshot file.
     if (noutsnaps > 0)
     {
         int snap = snapshot_steps[0];
         char fname_unsorted[256];
         char fname_sorted[256];

         char base_filename_1[256];
         snprintf(base_filename_1, sizeof(base_filename_1), "data/Rank_Mass_Rad_VRad_unsorted_t%05d.dat", snap);
         get_suffixed_filename(base_filename_1, 1, fname_unsorted, sizeof(fname_unsorted));
         char base_filename_2[256];
         snprintf(base_filename_2, sizeof(base_filename_2), "data/Rank_Mass_Rad_VRad_sorted_t%05d.dat", snap);
         get_suffixed_filename(base_filename_2, 1, fname_sorted, sizeof(fname_sorted));

         FILE *fun = fopen(fname_unsorted, "rb");
         FILE *fsort = fopen(fname_sorted, "rb");

         if (fun && fsort)
         {
             // Get file sizes.
             fseek(fun, 0, SEEK_END);
             reference_unsorted_size = ftell(fun);
             fseek(fsort, 0, SEEK_END);
             reference_sorted_size = ftell(fsort);

             // If both files have non-zero size, use as reference.
             if (reference_unsorted_size > 0 && reference_sorted_size > 0)
             {
                 printf("Found reference file sizes from snapshot %d (unsorted: %ld bytes, sorted: %ld bytes)\n\n",
                        snap, reference_unsorted_size, reference_sorted_size);
                 last_valid_snap = snap;
                 last_valid_index = 0;
                 checked_count += 2; // Count these two files.
             }

             fclose(fun);
             fclose(fsort);
         }
         else
         {
             // Only close if non-NULL.
             if (fun)
                 fclose(fun);
             if (fsort)
                 fclose(fsort);

             // If the very first files cannot be opened, processing should not continue.
             printf("Could not find initial snapshot files for index 0. Will start from the beginning.\n");
             return -1;
         }
     }
     else
     {
         // No snapshots to check.
         return -1;
     }

     // If reference sizes cannot be found, proper validation is not possible.
     if (reference_unsorted_size == 0 || reference_sorted_size == 0)
     {
         printf("Could not find valid reference file sizes. Will start from the beginning.\n\n");
         return -1;
     }

     printf("Checking all %d snapshots for completeness...\n\n", noutsnaps);

     // Now check ALL snapshots (except index 0 which we already checked).
     for (int i = 1; i < noutsnaps; i++)
     {
         int snap = snapshot_steps[i];

         // If this snapshot index maps to the same snapshot number as a previous index,.
         // We might be seeing duplicated snapshot numbers in the calculation.
         if (snap == snapshot_steps[0])
         {
             printf("Warning: Duplicate snapshot number %d (index 0 and %d)\n", snap, i);
         }

         char fname_unsorted[256];
         char fname_sorted[256];

         char base_filename_3[256];
         snprintf(base_filename_3, sizeof(base_filename_3), "data/Rank_Mass_Rad_VRad_unsorted_t%05d.dat", snap);
         get_suffixed_filename(base_filename_3, 1, fname_unsorted, sizeof(fname_unsorted));

         char base_filename_4[256];
         snprintf(base_filename_4, sizeof(base_filename_4), "data/Rank_Mass_Rad_VRad_sorted_t%05d.dat", snap);
         get_suffixed_filename(base_filename_4, 1, fname_sorted, sizeof(fname_sorted));

         FILE *fun = fopen(fname_unsorted, "rb");
         FILE *fsort = fopen(fname_sorted, "rb");

         if (fun && fsort)
         {
             // Get file sizes.
             fseek(fun, 0, SEEK_END);
             long unsorted_size = ftell(fun);
             fseek(fsort, 0, SEEK_END);
             long sorted_size = ftell(fsort);

             // Compare to reference sizes.
             // Allow for some small variation (±5%).
             double unsorted_ratio = (double)unsorted_size / reference_unsorted_size;
             double sorted_ratio = (double)sorted_size / reference_sorted_size;

             if (unsorted_size > 0 && sorted_size > 0 &&
                 unsorted_ratio >= 0.95 && unsorted_ratio <= 1.05 &&
                 sorted_ratio >= 0.95 && sorted_ratio <= 1.05)
             {
                 last_valid_snap = snap;
                 last_valid_index = i;
                 log_message("INFO", "Verified valid files for snapshot %d (unsorted: %ld bytes, sorted: %ld bytes)",
                             snap, unsorted_size, sorted_size);
                 checked_count += 2; // Count both files as checked.
             }
             else
             {
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
         }
         else
         {
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

     printf("End of file check: last_valid_snap=%d, last_valid_index=%d, noutsnaps=%d, all_snapshots_checked=%d\n\n",
            last_valid_snap, last_valid_index, noutsnaps, all_snapshots_checked);
     printf("Files checked: %d out of %d expected (unique snapshots: %d)\n\n",
            checked_count, total_expected, unique_snapshots);

     if (last_valid_snap == -1)
     {
         printf("No valid data products found. Will start from the beginning.\n\n");
         return -1;
     }

     // Check if we've verified more files than expected - this can happen if we have duplicate snapshot numbers.
     if (checked_count > total_expected)
     {
         printf("WARNING: Checked more files (%d) than expected (%d) - likely due to duplicate snapshot numbers.\n",
                checked_count, total_expected);
     }

     // Determine if we need to proceed with Rank file generation.
     if (unique_snapshots == 1 && all_snapshots_checked)
     {
         // Special case: Only one unique snapshot number (usually 0), and it's valid.
         printf("WARNING: Only one unique snapshot number found (%d). There's likely an issue with the calculation.\n",
                snapshot_steps[0]);
         printf("Only 1 Rank file (snapshot %d) exists. Starting from the beginning to create all files.\n",
                snapshot_steps[0]);
         return -1; // Start from beginning.
     }
     // Only say "all files exist" if:
     // 1. We have more than one unique snapshot, and.
     // 2. We've checked all expected files and found them valid.
     else if (unique_snapshots > 1 && checked_count >= total_expected && all_snapshots_checked)
     {
         printf("All %d data product files (for %d unique snapshots) already exist and are valid. Nothing to do.\n\n",
                checked_count, unique_snapshots);
         return -2; // Special code for "all done".
     }
     else
     {
         printf("Will restart processing from snapshot index %d (after snapshot %d)\n\n",
                last_valid_index + 1, last_valid_snap);
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
 void write_initial_conditions(double **particles, int npts, const char *filename)
 {
     FILE *fp = fopen(filename, "wb");
     if (!fp)
     {
         fprintf(stderr, "Error: cannot open '%s' for writing initial conditions.\n", filename);
         return;
     }

     // Write the number of particles first
     fwrite(&npts, sizeof(int), 1, fp);

     for (int i = 0; i < npts; i++)
     {
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
 void read_initial_conditions(double **particles, int npts, const char *filename)
 {
     FILE *fp = fopen(filename, "rb");
     if (!fp)
     {
         fprintf(stderr, "Error: cannot open '%s' for reading initial conditions.\n", filename);
         return;
     }

     int file_npts;
     if (fread(&file_npts, sizeof(int), 1, fp) != 1)
     {
         fprintf(stderr, "Error: failed to read npts from '%s'.\n", filename);
         fclose(fp);
         return;
     }
     if (file_npts != npts)
     {
         fprintf(stderr, "Warning: file npts=%d doesn't match current npts=%d.\n", file_npts, npts);
         fclose(fp);
         return;
     }

     for (int i = 0; i < npts; i++)
     {
         double r_val, v_val, ell_val, idx_val, mu_val;
         if (fread(&r_val, sizeof(double), 1, fp) != 1 || fread(&v_val, sizeof(double), 1, fp) != 1 || fread(&ell_val, sizeof(double), 1, fp) != 1 || fread(&idx_val, sizeof(double), 1, fp) != 1 || fread(&mu_val, sizeof(double), 1, fp) != 1)
         {
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
                                             float *Vrad_block)
 {
     FILE *f = fopen(filename, "ab");
     if (!f)
     {
         printf("Error: cannot open %s for appending all_particle_data\n", filename);
         CLEAN_EXIT(1);
     }

     // Write particle data in step-major order
     for (int step = 0; step < block_size; step++)
     {
         for (int i = 0; i < npts; i++)
         {
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
 void retrieve_all_particle_snapshot(const char *filename, int snap, int npts, int block_size, float *L_out, int *Rank_out, float *R_out, float *Vrad_out)
 {
     // Allocate local (thread-private) arrays
     float *tmpL = (float *)malloc(npts * sizeof(float));
     int *tmpRank = (int *)malloc(npts * sizeof(int));
     float *tmpR = (float *)malloc(npts * sizeof(float));
     float *tmpV = (float *)malloc(npts * sizeof(float));

     if (!tmpL || !tmpRank || !tmpR || !tmpV)
     {
         fprintf(stderr, "Error: out of memory in retrieve_all_particle_snapshot!\n");
         CLEAN_EXIT(1);
     }

     // Status messages are handled in the ordered section of the parallel loop.

 // Read from file in a critical section
 #pragma omp critical(file_access)
     {
         FILE *f = fopen(filename, "rb");
         if (!f)
         {
             fprintf(stderr, "Error: cannot open %s for reading\n", filename);
             CLEAN_EXIT(1);
         }

         // Compute offset in file.
         int block_number = snap / block_size;
         int index_in_block = snap % block_size;

         long long step_data_size = (long long)npts * 16; // 16 bytes per record.
         long long block_data_size = (long long)block_size * step_data_size;
         long long offset = block_data_size * block_number + step_data_size * index_in_block;

         if (fseek(f, offset, SEEK_SET) != 0)
         {
             fprintf(stderr, "Error: fseek failed for snap=%d\n", snap);
             fclose(f);
             CLEAN_EXIT(1);
         }

         // Read npts records into local buffers.
         for (int i = 0; i < npts; i++)
         {
             int rankval;
             float rval, vval, lval;

             if (fread(&rankval, sizeof(int), 1, f) != 1 ||
                 fread(&rval, sizeof(float), 1, f) != 1 ||
                 fread(&vval, sizeof(float), 1, f) != 1 ||
                 fread(&lval, sizeof(float), 1, f) != 1)
             {
                 fprintf(stderr, "Error: unexpected EOF while reading snap=%d (particle %d)\n", snap, i);
                 fclose(f);
                 CLEAN_EXIT(1);
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
