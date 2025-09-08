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

 #ifndef PARTICLE_ARRAY_OPS_H
 #define PARTICLE_ARRAY_OPS_H

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

int check_strict_monotonicity(const double *arr, int n, const char *name);

int compare_partdata_by_rad(const void *a, const void *b);
int compare_particles(const void *a, const void *b);
int compare_by_rr(const void *a, const void *b);

void sort_by_rad(struct PartData *array, int npts);
int double_cmp(const void *a, const void *b);
void insertion_sort(double **columns, int n);
void stdlib_qsort_wrapper(double **columns, int n);
void quadsort_wrapper(double **columns, int n);
void insertion_parallel_sort(double **columns, int n);
void quadsort_parallel_sort(double **columns, int n);

void verify_sort_results(double **columns, int n, const char *label);

void sort_particles_with_alg(double **particles, int npts, const char *sortAlg);
void sort_particles(double **particles, int npts);

void sort_rr_psi_arrays(double *rrA_spline, double *psiAarr_spline, int npts);

#endif // PARTICLE_ARRAY_OPS_H
