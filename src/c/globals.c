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

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

/**
 * @brief Global simulation feature flags.
 * @details These control various optional behaviors and optimizations
 *          in the simulation. Each flag can be set via command line arguments.
 */
int g_doDebug = 1;
int g_doDynPsi = 1;
int g_doDynRank = 1;
int g_doAllParticleData = 1;
int g_doRestart = 0;
int skip_file_writes = 0;
int g_enable_logging = 0;
int g_enable_sidm_scattering = 0;
int g_sidm_execution_mode = 1;
long long g_total_sidm_scatters = 0;
int g_sidm_kappa_provided = 0;

// Global variable definitions
double g_sidm_kappa = 50.0; ///< In cm^2/g
int *g_particle_scatter_state = NULL;

// Seed Management Globals
unsigned long int g_master_seed = 0;
unsigned long int g_initial_cond_seed = 0;
unsigned long int g_sidm_seed = 0;
int g_master_seed_provided = 0;
int g_initial_cond_seed_provided = 0;
int g_sidm_seed_provided = 0;
int g_attempt_load_seeds = 0;

const char* g_initial_cond_seed_filename_base = "data/last_initial_seed";
const char* g_sidm_seed_filename_base = "data/last_sidm_seed";
char g_file_suffix[256] = ""; ///< Global file suffix string.

// Profile parameter macros used by profile variables
#define RC 100.0                            ///< In kpc
#define RC_NFW_DEFAULT 1.18                 ///< In kpc
#define HALO_MASS_NFW 1.15e9                ///< In Msun
#define HALO_MASS 1.0e12                    ///< In Msun
#define CUTOFF_FACTOR_NFW_DEFAULT 85.0      ///< In Rc
#define CUTOFF_FACTOR_CORED_DEFAULT 85.0    ///< In Rc
#define FALLOFF_FACTOR_NFW_DEFAULT 19.0
#define NUM_MINI_SUBSTEPS_BOOTSTRAP 20

// Generalized Profile Parameters (set by new command line flags)
double g_scale_radius_param = RC;
double g_halo_mass_param = HALO_MASS;
double g_cutoff_factor_param = CUTOFF_FACTOR_CORED_DEFAULT;
char   g_profile_type_str[16] = "nfw";
int g_scale_radius_param_provided = 0;
int g_halo_mass_param_provided = 0;
int g_cutoff_factor_param_provided = 0;
double g_falloff_factor_param = FALLOFF_FACTOR_NFW_DEFAULT;
int    g_falloff_factor_param_provided = 0;
int g_profile_type_str_provided = 0;

// Profile Selection and NFW-Specific Parameters
int g_use_nfw_profile = 0;
double g_nfw_profile_rc = RC_NFW_DEFAULT;
double g_nfw_profile_halo_mass = HALO_MASS_NFW;
double g_nfw_profile_rmax_norm_factor = CUTOFF_FACTOR_NFW_DEFAULT;
double g_nfw_profile_falloff_factor = FALLOFF_FACTOR_NFW_DEFAULT;

// Cored Plummer-like Profile Specific Parameters (populated from generalized flags)
double g_cored_profile_rc = RC;
double g_cored_profile_halo_mass = HALO_MASS;
double g_cored_profile_rmax_factor = CUTOFF_FACTOR_CORED_DEFAULT;

/**
 * @brief Gravitational force calculation control flag.
 * @details Used for testing and debugging orbital dynamics:
 *          - 0 = Normal gravitational force calculation (default).
 *          - 1 = Zero gravity (particles move in straight lines).
 *
 * @note Setting this to 1 is useful for validating the integration scheme
 *       independent of gravitational physics.
 */
int use_identity_gravity = 0;

/**
 * @brief Convolution method selection for density smoothing.
 * @details 0 = FFT-based convolution (default, faster for large datasets).
 *          1 = Direct spatial convolution (more accurate but slower).
 */
int debug_direct_convolution = 0;

double g_active_halo_mass = HALO_MASS;

/**
 * @brief Angular momentum selection configuration for particle filtering.
 * @details Mode 0: Select particles with the 5 lowest L values.
 *          Mode 1: Select particles with L values closest to Lcompare.
 */
int use_closest_to_Lcompare = 1;
double Lcompare = 0.05;

/* ========================================================================= */

gsl_rng *g_rng = NULL;
gsl_rng **g_rng_per_thread = NULL;
int g_max_omp_threads_for_rng = 1;

// =========================================================================
// PERSISTENT SORT BUFFER CONFIGURATION
// =========================================================================
//
// The simulation exclusively uses a persistent global buffer (`g_sort_columns_buffer`)
// for particle data transposition during sorting operations. This strategy
// minimizes memory allocation/deallocation overhead.

/** Global persistent buffer for particle data transposition during sorting. */
double **g_sort_columns_buffer = NULL;
/** Number of particles the persistent buffer was allocated for; updated if npts changes. */
int g_sort_columns_buffer_npts = 0;

/** Global NFW mass spline for force calculations. */
gsl_spline *g_nfw_splinemass_for_force = NULL;
/** Global NFW mass spline accelerator for force calculations. */
gsl_interp_accel *g_nfw_enclosedmass_accel_for_force = NULL;

// =========================================================================
// SORT ALGORITHM CONFIGURATION
// =========================================================================
//
///< Default sorting algorithm identifier string. Set based on command-line options.
const char *g_defaultSortAlg = "quadsort_parallel";

// Constants controlling the behavior of parallel sorting algorithms.
// These parameters tune the parallel sorting operations used when
// OpenMP is available, affecting the partitioning of data across threads
// and the overlap required for correct merging of sorted sections.

/**
 * Default number of sections when OpenMP is unavailable or reports few threads.
 * Provides a baseline level of partitioning even in limited thread environments.
 */
const int PARALLEL_SORT_DEFAULT_SECTIONS = 8;

/**
 * Number of sort sections per OpenMP thread for workload distribution.
 * Multiplier used to determine total section count from available threads.
 */
const int PARALLEL_SORT_SECTIONS_PER_THREAD = 2;

/**
 * Divisor for calculating proportional overlap between sort sections.
 * Overlap is calculated as chunk_size / OVERLAP_DIVISOR to scale with data size.
 */
const int PARALLEL_SORT_OVERLAP_DIVISOR = 8;

/**
 * Minimum required overlap between adjacent sort sections (in elements).
 * This ensures sufficient overlap for correct merging of sorted sections,
 * even with sparse data distributions and when calculated proportional
 * overlap would be too small.
 */
const int PARALLEL_SORT_MIN_CORRECTNESS_OVERLAP = 32;

/**
 * Minimum average chunk size threshold for parallel sort operation.
 * If the calculated size of each chunk (total elements / number of sections)
 * falls below this threshold, the algorithm reverts to serial sorting to
 * avoid overhead from managing very small parallel tasks.
 */
const int PARALLEL_SORT_MIN_CHUNK_SIZE_THRESHOLD = 128;

// =========================================================================
// FILE I/O AND DATA MANAGEMENT SUBSYSTEM
// =========================================================================
//
// Functions for saving, loading, and managing simulation data including:
// - Initial condition generation and I/O
// - Snapshot file management
// - Binary file format utilities

int doReadInit = 0;
int doWriteInit = 0;
const char *readInitFilename = NULL;
const char *writeInitFilename = NULL;
