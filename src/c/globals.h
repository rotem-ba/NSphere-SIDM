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

#ifndef GLOBALS_H
#define GLOBALS_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

#ifdef _OPENMP
#include <omp.h>
#else
// OpenMP function stubs when compiled without OpenMP
// Use __attribute__((unused)) to prevent unused function warnings
#include <sys/time.h>  /* For gettimeofday */
static int __attribute__((unused)) omp_get_max_threads(void) { return 1; }
static int __attribute__((unused)) omp_get_num_procs(void) { return 1; }
static int __attribute__((unused)) omp_get_thread_num(void) { return 0; }
static void __attribute__((unused)) omp_set_max_active_levels(int x) { (void)x; }
static void __attribute__((unused)) omp_set_num_threads(int x) { (void)x; }
static double __attribute__((unused)) omp_get_wtime(void) {
    // Use higher precision time function for non-OpenMP builds
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}
#endif

// =========================================================================
// PHYSICAL CONSTANTS AND ASTROPHYSICAL PARAMETERS
// =========================================================================
//
// Core constants and unit conversion factors for astrophysical calculations
#define PI 3.14159265358979323846 ///< Mathematical constant Pi.
#define G_CONST 4.3e-6           ///< Newton's gravitational constant in kpc (km/sec)^2/Msun.
#define kmsec_to_kpcmyr 1.02271e-3 ///< Conversion factor: km/s to kpc/Myr.
#define VEL_CONV_SQ (kmsec_to_kpcmyr * kmsec_to_kpcmyr) ///< Velocity conversion squared (kpc/Myr)^2 per (km/s)^2.

/**
 * @brief Global simulation feature flags.
 * @details These control various optional behaviors and optimizations
 *          in the simulation. Each flag can be set via command line arguments.
 */
extern int g_doDebug;                   ///< Enable detailed debug output (default: on).
extern int g_doDynPsi;                  ///< Enable dynamic potential recalculation (default: on).
extern int g_doDynRank;                 ///< Enable dynamic rank calculation per step (default: on).
extern int g_doAllParticleData;         ///< Save complete particle evolution history (default: on).
extern int g_doRestart;                 ///< Enable simulation restart from checkpoint.
extern int skip_file_writes;            ///< Skip file writes during simulation restart.
extern int g_enable_logging;            ///< Enable logging to file (controlled by `--log` flag).
extern int g_enable_sidm_scattering;    ///< Enable SIDM scattering physics (0=no, 1=yes). Default is OFF.
extern int g_sidm_execution_mode;       ///< SIDM execution mode: 0 for serial, 1 for parallel (default).
extern long long g_total_sidm_scatters; ///< Global counter for total SIDM scatters.
extern int g_sidm_kappa_provided;       ///< Flag: 1 if `--sidm-kappa` was given by the user.

// Global variable definitions
extern double g_sidm_kappa;                     ///< SIDM opacity kappa (cm^2/g), default 50.0.
extern int *g_particle_scatter_state;           ///< Tracks recent scatter history for Adams-Bashforth integrator state reset.

// Seed Management Globals
extern unsigned long int g_master_seed;         ///< Master seed for the simulation, if provided.
extern unsigned long int g_initial_cond_seed;   ///< Seed used for generating initial conditions.
extern unsigned long int g_sidm_seed;           ///< Seed used for SIDM calculations.
extern int g_master_seed_provided;              ///< Flag: 1 if `--master-seed` was given by the user.
extern int g_initial_cond_seed_provided;        ///< Flag: 1 if `--init-cond-seed` was given by the user.
extern int g_sidm_seed_provided;                ///< Flag: 1 if `--sidm-seed` was given by the user.
extern int g_attempt_load_seeds;                ///< Flag: 1 if we should try to load seeds from files if not provided.

extern const char* g_initial_cond_seed_filename_base;     ///< Base name for IC seed file.
extern const char* g_sidm_seed_filename_base;             ///< Base name for SIDM seed file.
extern char g_file_suffix[256]; ///< Global file suffix string.

// Profile parameter macros used by profile variables
#define RC 100.0                            ///< Core radius in kpc.
#define RC_NFW_DEFAULT 1.18                 ///< Default NFW profile scale radius (kpc).
#define HALO_MASS_NFW 1.15e9                ///< Default NFW profile halo mass in solar masses (Msun).
#define HALO_MASS 1.0e12                    ///< Default general halo mass (used for Cored profile by default) in Msun.
#define CUTOFF_FACTOR_NFW_DEFAULT 85.0      ///< Default rmax factor for NFW profile (rmax = factor * rc).
#define CUTOFF_FACTOR_CORED_DEFAULT 85.0    ///< Default rmax factor for Cored profile (rmax = factor * rc).
#define FALLOFF_FACTOR_NFW_DEFAULT 19.0     ///< Default falloff factor for NFW profile (C_cutoff_factor).
#define NUM_MINI_SUBSTEPS_BOOTSTRAP 20      ///< Number of mini Euler steps per bootstrap full step.

// Profile Selection and NFW-Specific Parameters
extern int g_use_nfw_profile;               ///< Flag to use NFW-like profile for ICs: 0 = Cored (default), 1 = NFW.

extern double g_nfw_profile_rc;                 ///< NFW-specific scale radius (kpc); set by `--scale-radius` if NFW active, else defaults to RC_NFW_DEFAULT.
extern double g_nfw_profile_halo_mass;          ///< NFW-specific halo mass (Msun); set by `--halo-mass` if NFW active, else defaults to HALO_MASS_NFW.
extern double g_nfw_profile_rmax_norm_factor;   ///< NFW-specific r_max factor for IC norm/grid; set by `--cutoff-factor` if NFW active, else defaults to CUTOFF_FACTOR_NFW_DEFAULT.
extern double g_nfw_profile_falloff_factor;     ///< NFW-specific falloff transition C factor; set by `--falloff-factor` if NFW active, else defaults to FALLOFF_FACTOR_NFW_DEFAULT.


// Generalized Profile Parameters (set by new command line flags)
extern double g_scale_radius_param;         ///< Generalized scale radius (kpc), defaults to RC macro.
extern double g_halo_mass_param;            ///< Generalized halo mass (Msun), defaults to HALO_MASS macro.
extern double g_cutoff_factor_param;        ///< Generalized rmax factor, defaults to Cored's default.
extern char   g_profile_type_str[16];       ///< Profile type string ("nfw" or "cored"), default "nfw".

extern int g_scale_radius_param_provided;   ///< Flag: 1 if `--scale-radius` was given by the user.
extern int g_halo_mass_param_provided;      ///< Flag: 1 if `--halo-mass` was given by the user.
extern int g_cutoff_factor_param_provided;  ///< Flag: 1 if `--cutoff-factor` was given by the user.
extern double g_falloff_factor_param;       ///< Generalized falloff factor, defaults to NFW's default.
extern int g_falloff_factor_param_provided; ///< Flag: 1 if `--falloff-factor` was given by the user.
extern int g_profile_type_str_provided;     ///< Flag: 1 if `--profile` was given by the user.

// Cored Plummer-like Profile Specific Parameters (populated from generalized flags)
extern double g_cored_profile_rc;           ///< Cored-profile-specific scale radius (kpc); set by `--scale-radius` if Cored active, else defaults to RC macro.
extern double g_cored_profile_halo_mass;    ///< Cored-profile-specific halo mass (Msun); set by `--halo-mass` if Cored active, else defaults to HALO_MASS macro.
extern double g_cored_profile_rmax_factor;  ///< Cored-profile-specific r_max factor for IC norm/grid; set by `--cutoff-factor` if Cored active, else defaults to CUTOFF_FACTOR_CORED_DEFAULT.

/**
 * @brief Conditional compilation macros for feature flags.
 * @details These macros provide a cleaner syntax for conditional code blocks
 *          that depend on the global feature flags.
 *
 * @def IF_DEBUG
 * @brief Macro for code blocks executed only if `g_doDebug` is true.
 * @def IF_DYNPSI
 * @brief Macro for code blocks executed only if `g_doDynPsi` is true.
 * @def IF_DYNRANK
 * @brief Macro for code blocks executed only if `g_doDynRank` is true.
 * @def IF_ALL_PART
 * @brief Macro for code blocks executed only if `g_doAllParticleData` is true.
 */
#define IF_DEBUG if (g_doDebug)
#define IF_DYNPSI if (g_doDynPsi)
#define IF_DYNRANK if (g_doDynRank)
#define IF_ALL_PART if (g_doAllParticleData)

/**
 * @brief Gravitational force calculation control flag.
 * @details Used for testing and debugging orbital dynamics:
 *          - 0 = Normal gravitational force calculation (default).
 *          - 1 = Zero gravity (particles move in straight lines).
 *
 * @note Setting this to 1 is useful for validating the integration scheme
 *       independent of gravitational physics.
 */
extern int use_identity_gravity;

/**
 * @brief Convolution method selection for density smoothing.
 * @details 0 = FFT-based convolution (default, faster for large datasets).
 *          1 = Direct spatial convolution (more accurate but slower).
 */
extern int debug_direct_convolution;

extern double g_active_halo_mass; ///< Active halo mass for N-body force calculations.
extern double normalization; ///< Mass normalization factor for energy calculations. Calculated based on the integral of the density profile.

/**
 * @brief Angular momentum selection configuration for particle filtering.
 * @details Mode 0: Select particles with the 5 lowest L values.
 *          Mode 1: Select particles with L values closest to Lcompare.
 */
extern int use_closest_to_Lcompare; ///< Mode selector (0 or 1).
extern double Lcompare;          ///< Reference L value for closest-match mode (Mode 1).

/* ========================================================================= */

extern gsl_rng *g_rng; ///< GSL Random Number Generator state.
extern gsl_rng **g_rng_per_thread; ///< Array of GSL RNG states, one per OpenMP thread.
extern int g_max_omp_threads_for_rng; ///< Number of threads for which RNGs are allocated.

// =========================================================================
// SORTING ALGORITHM CONFIGURATION
// =========================================================================
///< Default sorting algorithm identifier string. Set based on command-line options.
extern const char *g_defaultSortAlg;

// =========================================================================
// PERSISTENT SORT BUFFER CONFIGURATION
// =========================================================================
//
// The simulation exclusively uses a persistent global buffer (`g_sort_columns_buffer`)
// for particle data transposition during sorting operations. This strategy
// minimizes memory allocation/deallocation overhead.

/** Global persistent buffer for particle data transposition during sorting. */
extern double **g_sort_columns_buffer;
/** Number of particles the persistent buffer was allocated for; updated if npts changes. */
extern int g_sort_columns_buffer_npts;

/** Global NFW mass spline for force calculations. */
extern gsl_spline *g_nfw_splinemass_for_force;
/** Global NFW mass spline accelerator for force calculations. */
extern gsl_interp_accel *g_nfw_enclosedmass_accel_for_force;

// =========================================================================
// PARALLEL SORT ALGORITHM CONFIGURATION
// =========================================================================
// Constants controlling the behavior of parallel sorting algorithms.
// These parameters tune the parallel sorting operations used when
// OpenMP is available, affecting the partitioning of data across threads
// and the overlap required for correct merging of sorted sections.

/**
 * Default number of sections when OpenMP is unavailable or reports few threads.
 * Provides a baseline level of partitioning even in limited thread environments.
 */
extern const int PARALLEL_SORT_DEFAULT_SECTIONS;

/**
 * Number of sort sections per OpenMP thread for workload distribution.
 * Multiplier used to determine total section count from available threads.
 */
extern const int PARALLEL_SORT_SECTIONS_PER_THREAD;

/**
 * Divisor for calculating proportional overlap between sort sections.
 * Overlap is calculated as chunk_size / OVERLAP_DIVISOR to scale with data size.
 */
extern const int PARALLEL_SORT_OVERLAP_DIVISOR;

/**
 * Minimum required overlap between adjacent sort sections (in elements).
 * This ensures sufficient overlap for correct merging of sorted sections,
 * even with sparse data distributions and when calculated proportional
 * overlap would be too small.
 */
extern const int PARALLEL_SORT_MIN_CORRECTNESS_OVERLAP;

/**
 * Minimum average chunk size threshold for parallel sort operation.
 * If the calculated size of each chunk (total elements / number of sections)
 * falls below this threshold, the algorithm reverts to serial sorting to
 * avoid overhead from managing very small parallel tasks.
 */
extern const int PARALLEL_SORT_MIN_CHUNK_SIZE_THRESHOLD;

// =========================================================================
// FILE I/O AND DATA MANAGEMENT SUBSYSTEM
// =========================================================================
//
// Functions for saving, loading, and managing simulation data including:
// - Initial condition generation and I/O
// - Snapshot file management
// - Binary file format utilities

extern int doReadInit;                  ///< Flag indicating whether to read initial conditions from file (1=yes, 0=no).
extern int doWriteInit;                 ///< Flag indicating whether to write initial conditions to file (1=yes, 0=no).
extern const char *readInitFilename;    ///< Filename to read initial conditions from (if doReadInit=1).
extern const char *writeInitFilename;   ///< Filename to write initial conditions to (if doWriteInit=1).

// =========================================================================
// RUNTIME PARAMETERS
// =========================================================================

extern int npts;                      ///< Number of particles in the simulation (default 100,000)
extern int Ntimes;                    ///< Number of time steps the simulation will run over (default 10,000)
extern int tfinal_factor;             ///< Number of dynamical times scales the simulation will run over (default 5)
extern int nout;                      ///< Number of output snapshots (default 100)
extern int dtwrite;                   ///< Number of steps between writes (default 100)
extern double tidal_fraction;         ///< Tidal stripping fraction (default 0.0 = off)
extern int method_select;             ///< Integration method used (default 1 = Adaptive Leapfrog with Adaptive Levi-Civita)
extern char *method_name;             ///< Name of the integration method used
extern char method_filename[32];      ///< Filename suffix for the integration method used
extern int display_sort;              ///< Sorting algorithm used (default 1 = Parallel Quadsort)
extern int include_method_in_suffix;  ///< Boolean, include method string in output filenames (default 0 = don't show)
extern char custom_tag[256];          ///< Tag added to output files (default: no custom tag)
extern char filename_tag[512];        ///< Filename total tag
extern char full_filename[256];       ///< Full filename (with suffixes)
extern char apd_filename[256];        ///< Filename for data/all_particle_data<suffix>.dat

extern int skip_simulation;           ///< Flag to determine if simulation phase can be skipped


#endif // GLOBALS_H
