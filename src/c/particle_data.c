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

#include <string.h>
#include <stdlib.h>
#include "logging.h"
#include "garbage_collection.h"

// =========================================================================
// GLOBAL PARTICLE DATA ARRAYS
// =========================================================================

/** @brief Global particle data arrays for snapshot processing. */
int *Rank_partdata_snap = NULL;
float *R_partdata_snap = NULL;
float *Vrad_partdata_snap = NULL;
float *L_partdata_snap = NULL;

/** @brief Global arrays for particle data processing (block storage). */
float *L_block = NULL;
int *Rank_block = NULL;
float *R_block = NULL;
float *Vrad_block = NULL;

/** @brief Variables for tracking low angular momentum particles. */
int nlowest = 5;
int *chosen = NULL;
double **lowestL_r = NULL;
double **lowestL_E = NULL;
double **lowestL_L = NULL;

double **particles = NULL;           ///< Main particle data array
double **baryons = NULL;             ///< Main baryons data array

/**
 * @brief Allocate main particle data array.
 */
void allocate_particles_memory(double **particles, int npts_initial) {
    particles = (double **)malloc(5 * sizeof(double *));
    if (particles == NULL)
        raise_error("ERROR: Memory allocation failed for particle array pointer\n");
    for (int i = 0; i < 5; i++) {
        particles[i] = (double *)malloc(npts_initial * sizeof(double));
        if (particles[i] == NULL)
            raise_error("ERROR: Memory allocation failed for particles[%d]\n", i);
    }
}


/**
 * @brief Allocate dark matter and baryon particle data array.
 */
void allocate_all_particles_memory(int npts_dark_matter, int npts_baryons) {
    allocate_particles_memory(particles, npts_dark_matter);
    allocate_particles_memory(baryons, npts_baryons);
}

/**
 * @brief Allocate new smaller arrays (`final_particles`) for the `npts` particles to keep.
 */
void trim_particles(double **particles, int npts) {
    double **final_particles = (double **)malloc(5 * sizeof(double *));
    if (final_particles == NULL)
        raise_error("Memory allocation failed for final_particles\n");

    /** @brief Copy innermost `npts` particles to final arrays and replace `particles` pointers. */
    for (int i = 0; i < 5; i++){ // Loop over components
        final_particles[i] = (double *)malloc(npts * sizeof(double));
        if (final_particles[i] == NULL)
            raise_error("Memory allocation failed for final_particles[%d]\n", i);
        /** @note Copy only the first `npts` elements (innermost after sort). */
        memcpy(final_particles[i], particles[i], npts * sizeof(double));

        /** @note Free original oversized array and update `particles[i]` pointer. */
        free(particles[i]);                // Free the original oversized array
        particles[i] = final_particles[i]; // particles[i] now points to the smaller array
    }
    free(final_particles); // Free the temporary ** structure, not the data arrays
}

void free_particles_memory(double **particles) {
    free_double_array(particles, 5);
}

void free_all_particles_memory() {
    free_particles_memory(particles);
    free_particles_memory(baryons);
}
