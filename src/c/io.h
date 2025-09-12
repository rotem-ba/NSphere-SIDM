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

#ifndef IO_H
#define IO_H

#include <string.h>
#include <stdio.h>
#include "globals.h"

// =========================================================================
// Windows‑compatibility shims
// =========================================================================
// Provide POSIX‑style helpers for MinGW/Clang:
//   • mkdir(path,mode)   → _mkdir(path)
//   • drand48 / srand48  → wrappers around ANSI rand
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
    #include <direct.h>
    #include <stdlib.h>

    /** Accept 1‑ or 2‑argument forms of mkdir on Windows. */
    #define mkdir(path, ...) _mkdir(path)

    #define NSPHERE_WINDOWS_SHIMS_DONE 1
#endif

long long get_available_disk_space(const char *path);
void get_full_filename(const char *base_filename, int with_suffix, char *buffer, size_t bufsize);
void format_file_size(long size_in_bytes, char *buffer, size_t buffer_size);
int fprintf_bin(FILE *fp, const char *format, ...);
int fscanf_bin(FILE *fp, const char *format, ...);

int find_last_processed_snapshot(int *snapshot_steps, int noutsnaps);
void write_initial_conditions(double **particles, int npts, const char *filename);
void read_initial_conditions(double **particles, int npts, const char *filename);

void append_all_particle_data_chunk_to_file(const char *filename, int npts, int block_size, float *L_block, int *Rank_block, float *R_block,
                                            float *Vrad_block);
void retrieve_all_particle_snapshot(const char *filename, int snap, int npts, int block_size, float *L_out, int *Rank_out, float *R_out, float *Vrad_out);

void compile_filename_tag();
void mkdir_init();
void write_to_lastparams();

void write_low_l_particles(double dt, int nlowest, double **lowestL_r, double **lowestL_E, double **lowestL_L);

/** @brief Calculate total number of write events and steps between major snapshots. */
inline int total_writes() {
    return ((Ntimes - 1) / dtwrite) + 1; // Total potential write points
}

void initialize_output_file(int noutsnaps);
void fill_suffix_tags();
int validate_snapshot_memory_allocation(double *r_grid, double *log_r_grid, double *mass_grid, double *density_grid, double *density_sorted,
                                        double *R_decimated, double *Mass_decimated, double *R_filtered, double *Mass_filtered, int r_violations, int snap);

void write_potential_profile(FILE *fp);
void write_mass_profile(FILE *fp);
void write_density_profile(FILE *fp);
void write_dPsidr_profile(FILE *fp);
void write_drhodPsi_profile(FILE *fp);
void write_f_of_E_profile(FILE *fp);
void write_distribution_after_simulation(FILE *fp);
void write_full_density_output(FILE *fp);

#endif // IO_H
