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
#include "globals.h"
#include "utils.h"

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
int adjust_ntimesteps(int Ntimes_initial, int nout, int dtwrite) {
    // Find the smallest N' >= N such that (Ntimes_initial' - 1) is a multiple of (nout - 1) * dtwrite.

    // Edge cases:
    if (nout < 2) // If only 0 or 1 snapshot requested, no interval constraint applies.
        return Ntimes_initial;
    if (dtwrite < 1) // Invalid write interval.
        return Ntimes_initial;

    // The total number of intervals between nout snapshots is (nout - 1).
    // The total number of steps spanning these intervals must be a multiple of dtwrite.
    // Therefore, the total number of steps (Ntimes_initial' - 1) must be a multiple of (nout - 1) * dtwrite.
    // Find the smallest integer k >= 1 such that (nout - 1) * k * dtwrite >= (Ntimes_initial - 1).
    double required_steps = (double)(Ntimes_initial - 1);
    double steps_per_output_cycle = (nout - 1) * (double)dtwrite;

    // Handle case where denominator is zero (e.g., nout=1 or dtwrite=0, caught above but added safety)
    if (steps_per_output_cycle <= 0)
        return Ntimes_initial; // Cannot satisfy constraint

    double ratio = required_steps / steps_per_output_cycle;
    int k = (int)ceil(ratio);
    if (k < 1)
        k = 1; // Ensure at least one full output cycle.

    int Nprime_minus_1 = (nout - 1) * k * dtwrite;
    int Nprime = Nprime_minus_1 + 1;

    return Nprime;
}

/**
 * @brief Adjust Ntimes using adjust_ntimesteps to align with output schedule.
 * @details Ensures (Ntimes - 1) is a multiple of (noutsnaps - 1) * dtwrite.
 * @see adjust_ntimesteps
 */
void adjust_Ntimes() {
    int oldN = Ntimes;
    Ntimes = adjust_ntimesteps(Ntimes, nout + 1, dtwrite); // Use noutsnaps here
    if (Ntimes != oldN)
        printf("Adjusted Number of Time Steps to %d to satisfy parameter constraints.\n", Ntimes);
    /** @brief Calculate total number of write events and steps between major snapshots. */
    ext_Ntimes = Ntimes + dtwrite; // Allocate trajectory arrays slightly larger
}

/**
 * @brief SIMULATION TIMESTEP CALCULATION block.
 * @details Calculates the characteristic dynamical time (`tdyn`) based on core radius (`RC`)
 *          and total mass (`HALO_MASS`). Uses this to determine the total simulation
 *          duration (`totaltime = tfinal_factor * tdyn`) and the individual timestep
 *          size (`dt = totaltime / Ntimes`) used in the integration loop.
 * @see tdyn
 * @see totaltime
 * @see dt
 */
void initialize_simulation_time() {
     double characteristic_radius_for_tdyn = g_use_nfw_profile ? g_nfw_profile_rc : g_cored_profile_rc;
     tdyn = 1.0 / sqrt((VEL_CONV_SQ * G_CONST) * g_active_halo_mass / cube(characteristic_radius_for_tdyn));
     totaltime = (double)tfinal_factor * tdyn; ///< Total simulation time (Myr)
     dt = totaltime / ((double)Ntimes);        ///< Individual timestep size (Myr)
     printf("Dynamical time tdyn = %.4f Myr\n", tdyn);
     printf("Total simulation time = %.4f Myr (%.1f tdyn)\n", totaltime, (double)tfinal_factor);
     printf("Timestep dt = %.6f Myr\n\n", dt);

     start_time = omp_get_wtime(); ///< Wall-clock start time for timing
     /** @brief Setup progress reporting steps (array `print_steps` holding step numbers for 0%, 5%, ..., 100%). */
     for (int k = 0; k <= 20; k++)
         print_steps[k] = (int)floor(k * 0.05 * Ntimes); // Calculate steps for progress output

}
