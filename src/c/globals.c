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

/**
 * @brief Global simulation feature flags.
 * @details These control various optional behaviors and optimizations
 *          in the simulation. Each flag can be set via command line arguments.
 */
int g_doDebug = 1;           ///< Enable detailed debug output (default: on).
int g_doDynPsi = 1;          ///< Enable dynamic potential recalculation (default: on).
int g_doDynRank = 1;         ///< Enable dynamic rank calculation per step (default: on).
int g_doAllParticleData = 1; ///< Save complete particle evolution history (default: on).
int g_doRestart = 0;         ///< Enable simulation restart from checkpoint.
int skip_file_writes = 0;    ///< Skip file writes during simulation restart.
int g_enable_logging = 0;    ///< Enable logging to file (controlled by `--log` flag).
extern int g_enable_sidm_scattering;    ///< Enable SIDM scattering physics (0=no, 1=yes). Default is OFF.
extern int g_sidm_execution_mode;       ///< SIDM execution mode: 0 for serial, 1 for parallel (default).
extern long long g_total_sidm_scatters; ///< Global counter for total SIDM scatters.

// Seed Management Globals
unsigned long int g_master_seed = 0;         ///< Master seed for the simulation, if provided.
unsigned long int g_initial_cond_seed = 0;   ///< Seed used for generating initial conditions.
unsigned long int g_sidm_seed = 0;           ///< Seed used for SIDM calculations.
int g_master_seed_provided = 0;              ///< Flag: 1 if `--master-seed` was given by the user.
int g_initial_cond_seed_provided = 0;        ///< Flag: 1 if `--init-cond-seed` was given by the user.

int g_attempt_load_seeds = 0;                ///< Flag: 1 if we should try to load seeds from files if not provided.

const char* g_initial_cond_seed_filename_base = "data/last_initial_seed"; ///< Base name for IC seed file.

// Profile parameter macros used by profile variables
#define RC 100.0                  ///< Core radius in kpc.
#define RC_NFW_DEFAULT 1.18       ///< Default NFW profile scale radius (kpc).
#define HALO_MASS_NFW 1.15e9      ///< Default NFW profile halo mass in solar masses (Msun).
#define HALO_MASS 1.0e12          ///< Default general halo mass (used for Cored profile by default) in Msun.
#define CUTOFF_FACTOR_NFW_DEFAULT 85.0  ///< Default rmax factor for NFW profile (rmax = factor * rc).
#define CUTOFF_FACTOR_CORED_DEFAULT 85.0 ///< Default rmax factor for Cored profile (rmax = factor * rc).
#define FALLOFF_FACTOR_NFW_DEFAULT 19.0 ///< Default falloff factor for NFW profile (C_cutoff_factor).
#define NUM_MINI_SUBSTEPS_BOOTSTRAP 20 ///< Number of mini Euler steps per bootstrap full step.

// Profile Selection and NFW-Specific Parameters
int g_use_nfw_profile = 0; ///< Flag to use NFW-like profile for ICs: 0 = Cored (default), 1 = NFW.

double g_nfw_profile_rc = RC_NFW_DEFAULT;    ///< NFW-specific scale radius (kpc); set by `--scale-radius` if NFW active, else defaults to RC_NFW_DEFAULT.
double g_nfw_profile_halo_mass = HALO_MASS_NFW; ///< NFW-specific halo mass (Msun); set by `--halo-mass` if NFW active, else defaults to HALO_MASS_NFW.
double g_nfw_profile_rmax_norm_factor = CUTOFF_FACTOR_NFW_DEFAULT; ///< NFW-specific r_max factor for IC norm/grid; set by `--cutoff-factor` if NFW active, else defaults to CUTOFF_FACTOR_NFW_DEFAULT.
double g_nfw_profile_falloff_factor = FALLOFF_FACTOR_NFW_DEFAULT; ///< NFW-specific falloff transition C factor; set by `--falloff-factor` if NFW active, else defaults to FALLOFF_FACTOR_NFW_DEFAULT.


// Generalized Profile Parameters (set by new command line flags)
double g_scale_radius_param = RC;        ///< Generalized scale radius (kpc), defaults to RC macro.
double g_halo_mass_param = HALO_MASS;    ///< Generalized halo mass (Msun), defaults to HALO_MASS macro.
double g_cutoff_factor_param = CUTOFF_FACTOR_CORED_DEFAULT; ///< Generalized rmax factor, defaults to Cored's default.
char   g_profile_type_str[16] = "nfw";   ///< Profile type string ("nfw" or "cored"), default "nfw".

int g_scale_radius_param_provided = 0;   ///< Flag: 1 if `--scale-radius` was given by the user.
int g_halo_mass_param_provided = 0;      ///< Flag: 1 if `--halo-mass` was given by the user.
int g_cutoff_factor_param_provided = 0;  ///< Flag: 1 if `--cutoff-factor` was given by the user.
double g_falloff_factor_param = FALLOFF_FACTOR_NFW_DEFAULT; ///< Generalized falloff factor, defaults to NFW's default.
int    g_falloff_factor_param_provided = 0;                  ///< Flag: 1 if `--falloff-factor` was given by the user.
int g_profile_type_str_provided = 0;     ///< Flag: 1 if `--profile` was given by the user.

// Cored Plummer-like Profile Specific Parameters (populated from generalized flags)
double g_cored_profile_rc = RC;              ///< Cored-profile-specific scale radius (kpc); set by `--scale-radius` if Cored active, else defaults to RC macro.
double g_cored_profile_halo_mass = HALO_MASS; ///< Cored-profile-specific halo mass (Msun); set by `--halo-mass` if Cored active, else defaults to HALO_MASS macro.
double g_cored_profile_rmax_factor = CUTOFF_FACTOR_CORED_DEFAULT;   ///< Cored-profile-specific r_max factor for IC norm/grid; set by `--cutoff-factor` if Cored active, else defaults to CUTOFF_FACTOR_CORED_DEFAULT.

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
 * @brief Convolution method selection for density smoothing.
 * @details 0 = FFT-based convolution (default, faster for large datasets).
 *          1 = Direct spatial convolution (more accurate but slower).
 */
int debug_direct_convolution = 0;
