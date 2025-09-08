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

#ifndef PARTICLE_DATA_H
#define PARTICLE_DATA_H

// =========================================================================
// GLOBAL PARTICLE DATA ARRAYS
// =========================================================================

/** @brief Global particle data arrays for snapshot processing. */
extern int *Rank_partdata_snap;   ///< Particle rank (sorted position) data for a snapshot.
extern float *R_partdata_snap;    ///< Radial position data for a snapshot.
extern float *Vrad_partdata_snap; ///< Radial velocity data for a snapshot.
extern float *L_partdata_snap;    ///< Angular momentum data for a snapshot.

/** @brief Global arrays for particle data processing (block storage). */
extern float *L_block;    ///< Angular momentum block.
extern int *Rank_block;   ///< Particle rank (sorted position) block.
extern float *R_block;    ///< Radial position block.
extern float *Vrad_block; ///< Radial velocity block.

/** @brief Variables for tracking low angular momentum particles. */
extern int nlowest;           ///< Number of lowest angular momentum particles to track.
extern int *chosen;        ///< Array of indices (original IDs) for selected low-L particles.
extern double **lowestL_r; ///< Radial positions of tracked low-L particles over time [particle][time_step].
extern double **lowestL_E; ///< Energy values of tracked low-L particles over time [particle][time_step].
extern double **lowestL_L; ///< Angular momenta of tracked low-L particles over time [particle][time_step].

#endif // PARTICLE_DATA_H
