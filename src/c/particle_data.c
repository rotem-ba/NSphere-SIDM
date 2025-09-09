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
