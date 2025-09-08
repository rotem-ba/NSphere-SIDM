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

 #ifndef EXIT_H
 #define EXIT_H

/**
 * @brief Frees all global arrays used for particle data processing.
 * @details Frees L_block, Rank_block, R_block, Vrad_block, and chosen.
 *          Does NOT free lowestL_* arrays as they are handled elsewhere.
 */
void cleanup_all_particle_data(void);

/**
 * @brief Frees local arrays used for snapshot processing.
 * @details Frees Rank_partdata_snap, R_partdata_snap, Vrad_partdata_snap, L_partdata_snap.
 */
void free_local_snap_arrays(void);

/**
 * @def CLEAN_EXIT(code)
 * @brief Thread-safe exit macro that properly cleans up allocated resources.
 *
 * @details This macro ensures proper resource cleanup before program termination:
 *          - Uses OpenMP critical section to ensure only one thread performs cleanup.
 *          - Calls cleanup_all_particle_data() and free_local_snap_arrays()
 *            to free dynamically allocated memory.
 *          - Exits with the specified error code.
 *
 * @param code The exit code for the program.
 *
 * @par Example
 * @code
 * if (error_condition) {
 *     CLEAN_EXIT(1);
 * }
 * @endcode
 */
#define CLEAN_EXIT(code)                 \
    do                                   \
    {                                    \
        _Pragma("omp critical")          \
        {                                \
            cleanup_all_particle_data(); \
            free_local_snap_arrays();    \
            exit(code);                  \
        }                                \
    } while (0)

#endif // EXIT_H
