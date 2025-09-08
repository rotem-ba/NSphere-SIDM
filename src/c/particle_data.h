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


/**
 * @brief Structure to track angular momentum with particle index and direction.
 * @details Used in sorting and selection of particles by angular momentum,
 *          especially when finding particles closest to a reference L.
 */
struct LAndIndex
{
    double L; ///< Angular momentum value (or squared difference from Lcompare).
    int idx;  ///< Original particle index (before sorting by L).
    int sign; ///< Direction indicator (+1 or -1) or sign of (L - Lcompare).
};
typedef struct LAndIndex LAndIndex;

/**
 * @brief Compact data structure for particle properties used in sorting and analysis.
 * @details Contains essential physical properties (radial position, velocity,
 *          angular momentum) and tracking metadata (rank, original index) for each
 *          particle in the simulation. Used extensively for sorting, file I/O, and
 *          data analysis operations.
 *
 * @note Uses compact `float` types for physical quantities to reduce memory usage
 *       when processing large particle counts.
 * @note The `rank` field is assigned during radial sorting, while `original_index`
 *       preserves the initial array position for tracking particles across snapshots.
 */
struct PartData
{
   int rank;           // Particle rank (sorted position)
   float rad;          // Radial position
   float vrad;         // Radial velocity
   float angmom;       // Angular momentum
   int original_index; // Original position in array before sorting
};
typedef struct PartData PartData;

struct RrPsiPair
{
   double rr;  ///< Radius value or x-axis value for sorting
   double psi; ///< Corresponding potential value or y-axis value
};
typedef struct RrPsiPair RrPsiPair;


#endif // PARTICLE_DATA_H
