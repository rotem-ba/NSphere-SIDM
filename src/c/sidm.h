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

#ifndef SIDM_H
#define SIDM_H

#include "globals.h"
#include "utils.h"
#include <gsl/gsl_rng.h>
#include <math.h>

/**
 * @brief Structure for buffering SIDM scatter events in parallel execution.
 * @details Used to store scatter event information before applying updates
 *          to the main particle array, ensuring thread safety in parallel SIDM.
 */
struct ScatterEvent{
    int i;              ///< Index of the first particle in the scattering pair.
    int m_offset;       ///< Offset of the scattering partner relative to particle i.
    threevector Vifinal;///< Final 3D velocity vector of particle i.
    threevector Vmfinal;///< Final 3D velocity vector of the partner particle.
};
typedef struct ScatterEvent ScatterEvent;

/**
 * @brief Global simulation feature flags.
 * @details These control various optional behaviors and optimizations
 *          in the simulation. Each flag can be set via command line arguments.
 */

// Global variables (extern declarations - defined in sidm.c)
extern double g_sidm_kappa;                    ///< SIDM opacity kappa (cm^2/g).
extern int *g_particle_scatter_state;          ///< Particle scatter state tracking.

// SIDM physics functions
double sigmatotal(double vrel, int npts);

// SIDM scattering functions
void handle_sidm_step();

void perform_sidm_scattering_serial(gsl_rng *rng, long long *Nscatter_total_step);

void perform_sidm_scattering_parallel(gsl_rng **rng_per_thread_list, int num_threads_for_rng, long long *Nscatter_total_step);

int compare_scatter_events(const void *a, const void *b);

#endif // SIDM_H
