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

 #include <math.h>

 /**
  * @brief Adjusts the total number of timesteps to align with desired output snapshot intervals.
  * @details This function calculates an adjusted number of total simulation timesteps, \f$N'_{times}\f$,
  *          such that it is greater than or equal to the initially requested `Ntimes_initial` (\f$N\f$)
  *          and satisfies the constraint: \f$(N'_{times} - 1)\f$ must be an integer multiple of
  *          \f$(M - 1) \times p\f$. Here, \f$M\f$ is `nout` (number of desired output snapshot points,
  *          which means \f$M-1\f$ intervals) and \f$p\f$ is `dtwrite` (the low-level write interval
  *          in terms of simulation timesteps).
  *          This alignment ensures that exactly `nout` snapshots can be produced at intervals
  *          that are multiples of `dtwrite` and that also evenly span the total adjusted simulation duration.
  *
  * @param Ntimes_initial [in] Initially requested total number of simulation timesteps (\f$N\f$).
  * @param nout           [in] Number of desired output snapshot points (\f$M\f$). Must be >= 2 for adjustment to apply.
  * @param dtwrite        [in] The interval (in timesteps) at which low-level data is potentially written (\f$p\f$). Must be >= 1.
  * @return int The adjusted total number of timesteps (\f$N'_{times}\f$). Returns `Ntimes_initial`
  *             if `nout < 2` or `dtwrite < 1` or other edge cases where the constraint cannot be met.
  */
int adjust_ntimesteps(int Ntimes_initial, int nout, int dtwrite)
{
     // Notation: M = nout, p = dtwrite, N = Ntimes_initial.
     // Find the smallest N' >= N such that (N' - 1) is a multiple of (M - 1) * p.

     int M = nout;
     int p = dtwrite;
     int N = Ntimes_initial;

     // Edge cases:
     if (M < 2)
     {
         // If only 0 or 1 snapshot requested, no interval constraint applies.
         return N;
     }
     if (p < 1)
     {
         // Invalid write interval.
         return N;
     }

     // The total number of intervals between M snapshots is (M - 1).
     // The total number of steps spanning these intervals must be a multiple of p.
     // Therefore, the total number of steps (N' - 1) must be a multiple of (M - 1) * p.
     // Find the smallest integer k >= 1 such that (M - 1) * k * p >= (N - 1).
     double required_steps = (double)(N - 1);
     double steps_per_output_cycle = (M - 1) * (double)p;

     // Handle case where denominator is zero (e.g., M=1 or p=0, caught above but added safety)
     if (steps_per_output_cycle <= 0) {
         return N; // Cannot satisfy constraint
     }

     double ratio = required_steps / steps_per_output_cycle;
     int k = (int)ceil(ratio);
     if (k < 1)
     {
         k = 1; // Ensure at least one full output cycle.
     }

     int Nprime_minus_1 = (M - 1) * k * p;
     int Nprime = Nprime_minus_1 + 1;

     return Nprime;
 }
