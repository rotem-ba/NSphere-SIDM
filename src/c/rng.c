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
#include "io.h"

/**
 * @brief Initialize random number generator for particle generation.
 * @details Sets up the GSL Random Number Generator environment and allocates
 *          the global GSL RNG state used throughout the simulation.
 *          Needed for the Sample Generator if not reading ICs from file.
 */
void init_rng() {
    gsl_rng_env_setup();                           // Setup GSL RNG environment
    const gsl_rng_type *T_rng = gsl_rng_default;  // Use default RNG type
    g_rng = gsl_rng_alloc(T_rng);                  // Allocate global GSL RNG state
    if (g_rng == NULL)                             // Check allocation success
        raise_error("Error allocating GSL RNG.\n");

    // Seed the global g_rng (used for ICs and Serial SIDM)
    gsl_rng_set(g_rng, g_initial_cond_seed);       // Use the determined IC seed for g_rng
    log_message("INFO", "Global g_rng (intended primarily for IC generation) seeded with %lu", g_initial_cond_seed);
}

/**
 * @brief Initialize random number generator per thread..
 */
void init_rng_per_thread() {
    #ifdef _OPENMP
        g_max_omp_threads_for_rng = omp_get_max_threads();
        if (g_max_omp_threads_for_rng <= 0) g_max_omp_threads_for_rng = 1; // Safety
    #else
        g_max_omp_threads_for_rng = 1;
    #endif

    g_rng_per_thread = (gsl_rng **)malloc(g_max_omp_threads_for_rng * sizeof(gsl_rng *));
    if (g_rng_per_thread == NULL)
        raise_error("Error: Failed to allocate memory for per-thread RNG array.\n");

    const gsl_rng_type *T_rng_thread = gsl_rng_default;

    for (int i_rng = 0; i_rng < g_max_omp_threads_for_rng; ++i_rng) {
        g_rng_per_thread[i_rng] = gsl_rng_alloc(T_rng_thread);
        if (g_rng_per_thread[i_rng] == NULL) {
            for (int k_rng = 0; k_rng < i_rng; ++k_rng) gsl_rng_free(g_rng_per_thread[k_rng]);
            free(g_rng_per_thread);
            raise_error("Error: Failed to allocate GSL RNG for thread %d.\n", i_rng);
        }
        gsl_rng_set(g_rng_per_thread[i_rng], g_sidm_seed + (unsigned long int)i_rng);
    }
    log_message("INFO", "Initialized %d per-thread GSL RNGs (for SIDM) using base SIDM seed %lu", g_max_omp_threads_for_rng, g_sidm_seed);
}

/**
 * @brief Determine SIDM Seed
 */
void set_sidm_seed(unsigned long current_time_pid_seed, char *seed_filepath, FILE * fp_seed) {
    if (!g_sidm_seed_provided) {
        if (g_master_seed_provided)
            g_sidm_seed = g_master_seed + 2; // Deterministic offset, different from IC seed
        else if (g_attempt_load_seeds) {
            get_full_filename(g_sidm_seed_filename_base, 1, seed_filepath, sizeof(seed_filepath));
            fp_seed = fopen(seed_filepath, "r");
            if (fp_seed) {
                if (fscanf(fp_seed, "%lu", &g_sidm_seed) == 1)
                    log_message("INFO", "Loaded SIDM seed %lu from %s", g_sidm_seed, seed_filepath);
                else {
                    g_sidm_seed = current_time_pid_seed + 200; // Fallback
                    log_message("WARNING", "Failed to read SIDM seed from %s, generating new: %lu", seed_filepath, g_sidm_seed);
                }
                fclose(fp_seed);
            } else {
                g_sidm_seed = current_time_pid_seed + 200; // File not found, generate
                log_message("INFO", "No SIDM seed file found, generating new: %lu", g_sidm_seed);
            }
        } else {
            g_sidm_seed = current_time_pid_seed + 200; // Default generation
            log_message("INFO", "Generating new SIDM seed: %lu", g_sidm_seed);
        }
    } else
        log_message("INFO", "Using user-provided SIDM seed: %lu", g_sidm_seed);
}

/**
 * @brief Determine generic Seed
 */
void set_seed(unsigned long *seed, char *seed_type, int seed_provided, const char *seed_filename_base, unsigned long current_time_pid_seed) {
    if (seed_provided)
        log_message("INFO", "Using user-provided %s seed: %lu", seed_type, *seed);
    else if (g_master_seed_provided)
        *seed = g_master_seed + 2; // Deterministic offset, different from IC seed
    else if (!g_attempt_load_seeds) {
        *seed = current_time_pid_seed + 200; // Default generation
        log_message("INFO", "Generating new %s seed: %lu", seed_type, *seed);
    } else {
        char seed_filepath[512];
        FILE *fp_seed;
        get_full_filename(seed_filename_base, 1, seed_filepath, sizeof(seed_filepath));
        fp_seed = fopen(seed_filepath, "r");
        if (!fp_seed) {
            *seed = current_time_pid_seed + 200; // File not found, generate
            log_message("INFO", "No %s seed file found, generating new: %lu", seed_type, *seed);
        } else {
            if (fscanf(fp_seed, "%lu", seed) == 1)
                log_message("INFO", "Loaded %s seed %lu from %s", seed_type, seed, seed_filepath);
            else {
                *seed = current_time_pid_seed + 200; // Fallback
                log_message("WARNING", "Failed to read %s seed from %s, generating new: %lu", seed_type, seed_filepath, *seed);
            }
            fclose(fp_seed);
        }
    }
}
