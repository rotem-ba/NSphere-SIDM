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
#include "logging.h"
#include "io.h"
#include "debug.h"
#include "density.h"
#include <math.h>

// =========================================================================
// Energy Debugging and Validation Subsystem
// =========================================================================


/** @brief Arrays for tracking energy components through simulation for debugging. */
double dbg_approxE[DEBUG_MAX_STEPS]; ///< Theoretical model energy (per unit mass).
double dbg_dynE[DEBUG_MAX_STEPS];    ///< Actual dynamical energy (per unit mass).
double dbg_kinE[DEBUG_MAX_STEPS];    ///< Kinetic energy component (per unit mass).
double dbg_potE[DEBUG_MAX_STEPS];    ///< Potential energy component (per unit mass).
double dbg_time[DEBUG_MAX_STEPS];    ///< Simulation time at each snapshot (Myr).
double dbg_radius[DEBUG_MAX_STEPS];  ///< Particle radius at each snapshot (kpc).
int dbg_count = 0;                   ///< Number of debug snapshots recorded.

/**
 * @brief Records the theoretical model energy for a debug snapshot.
 *
 * Parameters
 * ----------
 * snapIndex : int
 *     Index of the snapshot (0 to DEBUG_MAX_STEPS - 1).
 * E_value : double
 *     Theoretical energy value (per unit mass).
 * time_val : double
 *     Simulation time (Myr).
 *
 * Returns
 * -------
 * None
 *
 * @note Updates `dbg_count` if `snapIndex` extends the recorded range.
 * @see store_debug_dynE_components
 * @see finalize_debug_energy_output
 */
void store_debug_approxE(int snapIndex, double E_value, double time_val) {
    if (snapIndex < 0 || snapIndex >= DEBUG_MAX_STEPS)
        return;

    dbg_approxE[snapIndex] = E_value;
    dbg_time[snapIndex] = time_val;

    // Update the total count if needed
    if (snapIndex >= dbg_count)
        dbg_count = snapIndex + 1;
}

/**
 * @brief Records the actual dynamical energy and its components for a debug snapshot.
 *
 * Parameters
 * ----------
 * snapIndex : int
 *     Index of the snapshot (0 to DEBUG_MAX_STEPS - 1).
 * totalE : double
 *     Total energy (KE + PE) per unit mass.
 * kinE : double
 *     Kinetic energy component (per unit mass).
 * potE : double
 *     Potential energy component (per unit mass).
 * time_val : double
 *     Simulation time (Myr).
 * radius_val : double
 *     Particle radius (kpc).
 *
 * Returns
 * -------
 * None
 *
 * @note Updates `dbg_count` if `snapIndex` extends the recorded range.
 * @see store_debug_approxE
 * @see finalize_debug_energy_output
 */
void store_debug_dynE_components(int snapIndex, double totalE, double kinE, double potE, double time_val, double radius_val)
{
    if (snapIndex < 0 || snapIndex >= DEBUG_MAX_STEPS)
        return;

    dbg_dynE[snapIndex] = totalE;
    dbg_kinE[snapIndex] = kinE;
    dbg_potE[snapIndex] = potE;
    dbg_time[snapIndex] = time_val;
    dbg_radius[snapIndex] = radius_val;

    // Update the total count if needed
    if (snapIndex >= dbg_count)
        dbg_count = snapIndex + 1;
}

/**
 * @brief Writes the collected debug energy comparison data to a file.
 * @details This function is called at the end of the simulation if debug mode (`g_doDebug`)
 *          is enabled. It iterates through the stored debug snapshots collected during
 *          the simulation and post-processing, writing the following data for the
 *          tracked particle (DEBUG_PARTICLE_ID) to `data/debug_energy_compare.dat`
 *          (with the appropriate file suffix applied):
 *          - Snapshot index
 *          - Simulation time (Myr)
 *          - Particle radius (kpc)
 *          - Approximate (theoretical) energy (per unit mass)
 *          - Dynamic (calculated) energy (per unit mass)
 *          - Difference between dynamic and approximate energy
 *          - Kinetic energy component (per unit mass)
 *          - Potential energy component (per unit mass)
 *
 * @note The file writing operation is skipped if the simulation is in restart mode
 *       and file writes are disabled (`skip_file_writes` is true).
 *
 * @see store_debug_approxE
 * @see store_debug_dynE_components
 */
void finalize_debug_energy_output(void)
{
    // Skip file writes if in restart mode and file writes are skipped
    if (!skip_file_writes)
    {
        // Create filename with suffix
        char filename[256];
        get_full_filename("data/debug_energy_compare.dat", 1, filename, sizeof(filename));

        FILE *fp = fopen(filename, "wb"); // Binary mode for fprintf_bin
        if (!fp) {
            // Use log_message for consistency if logging is enabled
            log_message("ERROR", "Cannot open %s for writing debug energy data.", filename);
            // Fallback to printf if logging might be off or failed
            printf("Error: cannot open %s\n", filename);
            return;
        }

        /* Output 8 columns, including KE & PE. */
        fprintf_bin(fp, "# snapIdx  time(Myr)  radius(kpc)   approxE   dynE   (dyn-approx)    KE       PE\n");

        for (int i = 0; i < dbg_count; i++) {
            double eA = dbg_approxE[i];
            double eD = dbg_dynE[i];
            double ke = dbg_kinE[i];
            double pe = dbg_potE[i];
            double tVal = dbg_time[i];
            double rVal = dbg_radius[i];

            double dff = eD - eA;
            fprintf_bin(fp, "%4d  %.6f  %.6f  %.8g  %.8g  %.8g  %.8g  %.8g\n", i, tVal, rVal, eA, eD, dff, ke, pe);
        }
        fclose(fp);

        // Only log the debug message, never print to console
        log_message("DEBUG", "Wrote %s with %d lines.", filename, dbg_count);
    } else
        log_message("INFO", "Skipped writing debug energy file due to restart/skip_file_writes flag.");
}

/**
 * @brief Log debug Psi spline values.
 */
void debug_log_Psivalues() {
    if (!g_doDebug)
        return;
    log_message("DEBUG", "Psi(r) and r(Psi) spline data summary (num_points=%d):", num_points);
    if (num_points > 1) {
        log_message("DEBUG", "  Psivalues[0] (Psimax candidate) = %.6e", Psivalues[0]);
        log_message("DEBUG", "  Psivalues[num_points-1] (Psimin candidate) = %.6e", Psivalues[num_points-1]);
        if (Psivalues[0] <= Psivalues[num_points-1])
            log_message("WARNING", "Psi(r) may not be monotonic decreasing (Psivalues[0]=%.3e <= Psivalues[end]=%.3e)!", Psivalues[0], Psivalues[num_points-1]);
    }
    log_message("DEBUG", "End of Psi(r) data summary.");
}

/**
 * @brief Log debug NFW M spline values.
 */
void debug_log_nfw_M(double current_profile_halo_mass, const double nfw_params[4]) {
    if (!g_doDebug)
        return;
    log_message("DEBUG", "M(r) spline data summary (num_points=%d):", num_points);
    log_message("DEBUG", "  Target M_total for sampling = %.3e Msun", current_profile_halo_mass);
    log_message("DEBUG", "  nt_nfw used for M(r) calcs = %.3e", nfw_params[2]);
    log_message("DEBUG", "  rmax for M(r) array = %.3e kpc", rmax);
    if (num_points > 0) {
        log_message("DEBUG", "  Final Mass at rmax (radius[num_points-1]=%.3e kpc): %.3e Msun", radius[num_points-1], mass[num_points-1]);
        if (fabs(mass[num_points-1] - current_profile_halo_mass) / current_profile_halo_mass > 0.1)
            log_message("WARNING", "Mass at rmax (%.3e) differs significantly from target halo mass (%.3e)!", mass[num_points-1], current_profile_halo_mass);
    }
    log_message("DEBUG", "End of M(r) data summary.");
}
