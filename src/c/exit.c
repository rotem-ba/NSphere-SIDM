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

#include <stdlib.h>
#include "particle_data.h"

/**
 * @brief Frees all global arrays used for particle data processing.
 * @details Frees L_block, Rank_block, R_block, Vrad_block, and chosen.
 *          Does NOT free lowestL_* arrays as they are handled elsewhere.
 */
void cleanup_all_particle_data(void){
    free(L_block);
    free(Rank_block);
    free(R_block);
    free(Vrad_block);

    // Free low angular momentum tracking arrays
    free(chosen);

    // lowestL arrays are freed in the main function after they're used
}

/**
 * @brief Frees local arrays used for snapshot processing.
 * @details Frees Rank_partdata_snap, R_partdata_snap, Vrad_partdata_snap, L_partdata_snap.
 */
void free_local_snap_arrays(void) {
    free(Rank_partdata_snap);
    free(R_partdata_snap);
    free(Vrad_partdata_snap);
    free(L_partdata_snap);
}
