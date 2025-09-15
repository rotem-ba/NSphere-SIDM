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
#include "particle_data.h"
#include <math.h>

/**
 * @brief TIDAL STRIPPING IMPLEMENTATION block.
 * @details If `tidal_fraction` > 0 and not in restart mode, simulates tidal stripping by removing the outermost
 *          fraction of particles based on radius. It first sorts the `npts_initial`
 *          particles by radius, then keeps only the innermost `npts` particles.
 *          It reallocates the `particles` array to the final size `npts` and remaps
 *          the original indices stored in `particles[3]` to ranks [0, npts-1] using
 *          `reassign_orig_ids_with_rank`.
 * @see reassign_orig_ids_with_rank
 * @see sort_particles_with_alg
 */
void tidal_strip(int npts_initial) {
    if (g_doRestart)
        return;
    if (tidal_fraction == 0.0) {
        sort_particles_with_alg(particles, npts, "quadsort");
    }
    else {
        printf("Tidal stripping: sorting and retaining inner %.1f%% of particles...\n", (1.0 - tidal_fraction) * 100.0);
        sort_particles_with_alg(particles, npts_initial, "quadsort"); // Sorts by particles[0]
        trim_particles(npts);
        printf("Tidal stripping complete: %d particles retained.\n\n", npts);
    }
    /** @note Remap original IDs (now in `particles[3]` for the kept particles) to ranks [0, npts-1]. */
    reassign_orig_ids_with_rank(particles[3], npts);
}

/**
 * @brief Oversample initial conditions based on tidal fraction.
 * @details Calculate the number of initial particles (`npts_initial`) needed
 *          before tidal stripping to ensure `npts` particles remain afterwards.
 *          If `tidal_fraction` is 0, `npts_initial` equals `npts`.
 */
int get_npts_initial_with_tidal_fraction() {
    return (tidal_fraction > 0.0) ? ceil(npts / (1.0 - tidal_fraction)) : npts;
}
