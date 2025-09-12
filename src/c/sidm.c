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
#include "particle_data.h"
#include "sidm.h"
#include "logging.h"
#include "utils.h"
#include "density.h"
#include <math.h>
#include <gsl/gsl_rng.h>

/**
 * @brief Calculates the total SIDM cross-section for a given relative velocity.
 * @details Computes the self-interaction cross-section using the opacity parameter
 *          kappa and the particle mass derived from the total halo mass and number
 *          of particles.
 *
 * @param vrel [in] Relative velocity between particles (unused in current implementation).
 * @param npts [in] Number of particles in the simulation.
 * @return double Cross-section in kpc^2.
 */
double sigmatotal(double vrel __attribute__((unused)), int npts) {
    double kappa = g_sidm_kappa; // Self-interaction opacity parameter (cm²/g)
    // double active_profile_rc = density_rc(); //Unused at the moment
    // Ensure npts is positive to prevent division by zero or negative particle mass
    if (npts <= 0)
        return 0.0;
    double particle_mass_Msun = g_halo_mass_param / ((double)npts);
    if (particle_mass_Msun <= 0)
        return 0.0;
    return 2.089e-10 * kappa * particle_mass_Msun; // Cross-section (kpc²)
}

/**
 * @brief Handles the SIDM scattering phase for a single timestep.
 * @details Checks if SIDM is enabled and if not in a bootstrap phase that should skip SIDM.
 *          If proceeding, it resets particle scatter flags, selects serial or parallel execution
 *          based on `g_sidm_execution_mode`, calls the appropriate core scattering function
 *          (`perform_sidm_scattering_serial` or `perform_sidm_scattering_parallel`),
 *          updates the global total scatter count, and logs debug information if scatters occurred
 *          and debugging is enabled. The core scattering functions are responsible for updating
 *          the `g_particle_scatter_state` flags for particles that underwent scattering.
 *
 * @note This function modifies the `particles` array in-place.
 * @note It uses global variables: `g_enable_sidm_scattering`, `g_sidm_execution_mode`,
 *       `g_rng_per_thread`, `g_max_omp_threads_for_rng`, `g_rng`, `g_total_sidm_scatters`,
 *       `g_active_halo_mass`, `g_doDebug`, and `g_particle_scatter_state`.
 */
void handle_sidm_step() {
    if (!g_enable_sidm_scattering || !bootstrap_phase_done)
        return; // Skip SIDM if disabled or in a bootstrap phase that should skip SIDM

    long long Nscatters_in_this_step = 0;

    if (g_sidm_execution_mode == 1) { // Parallel
        #ifdef _OPENMP
            if (g_rng_per_thread != NULL && g_max_omp_threads_for_rng > 0)
                perform_sidm_scattering_parallel(g_rng_per_thread, g_max_omp_threads_for_rng, &Nscatters_in_this_step);
            else {
                log_message("ERROR", "SIDM Parallel mode selected but per-thread RNGs not available. Skipping SIDM for step.");
                Nscatters_in_this_step = 0;
            }
        #else
            // Serial fallback if OpenMP not compiled but parallel mode selected
            log_message("WARNING", "SIDM Parallel mode selected but OpenMP not enabled. Running SIDM serially.");
            gsl_rng *rng_for_serial_fallback = (g_rng_per_thread != NULL && g_rng_per_thread[0] != NULL) ? g_rng_per_thread[0] : g_rng;
            if (rng_for_serial_fallback != NULL)
                perform_sidm_scattering_serial(rng_for_serial_fallback, &Nscatters_in_this_step);
            else {
                log_message("ERROR", "SIDM Serial fallback: No suitable RNG available. Skipping SIDM for step.");
                Nscatters_in_this_step = 0;
            }
        #endif
    } else { // Serial SIDM execution
        gsl_rng *rng_for_serial = (g_rng_per_thread != NULL && g_rng_per_thread[0] != NULL) ? g_rng_per_thread[0] : g_rng;
        if (rng_for_serial != NULL)
            perform_sidm_scattering_serial(rng_for_serial, &Nscatters_in_this_step);
        else {
            log_message("ERROR", "SIDM Serial mode: No suitable RNG available. Skipping SIDM for step.");
            Nscatters_in_this_step = 0;
        }
    }

    g_total_sidm_scatters += Nscatters_in_this_step;

    if (Nscatters_in_this_step > 0 && g_doDebug)
        log_message("DEBUG", "Method %d Step: %lld SIDM scatters this step, %lld total", method_select, Nscatters_in_this_step,
                    g_total_sidm_scatters);
}

/**
 * @brief Executes self-interacting dark matter (SIDM) scattering for one simulation timestep (Serial version).
 * @details Implements a serial SIDM scattering algorithm. For each particle `i` (the primary scatterer),
 *          it considers up to `nscat` (typically 10) subsequent particles in the array as potential
 *          scattering partners. (The `particles` array is assumed to be sorted by radius, so these are
 *          spatially nearby neighbors). The algorithm proceeds as follows for each primary particle `i`:
 *
 *          1. **Velocity Construction**: Constructs 3D velocity vectors for the primary particle and
 *             its potential scattering partners using radial velocity and angular momentum components.
 *          2. **Interaction Rate Calculation**: Computes the scattering rate (σ × v_rel) for each
 *             potential partner based on their relative velocity and the SIDM cross-section.
 *          3. **Probability Calculation**: Determines the total scattering probability using a shell
 *             volume approximation and the timestep.
 *          4. **Stochastic Selection**: Uses random sampling to decide whether a scatter occurs and,
 *             if so, selects the scattering partner based on weighted probabilities.
 *          5. **Scattering Dynamics**: Implements isotropic scattering in the center-of-mass frame,
 *             transforming the scattered velocities back to the lab frame.
 *          6. **State Update**: Updates particle velocities and angular momenta, and marks scattered
 *             particles for potential integrator state reset.
 *
 * @param particles [in,out] The main particle data array: `particles[component][current_sorted_index]`.
 *                           - Component 0: radial position (kpc)
 *                           - Component 1: radial velocity (km/s)
 *                           - Component 2: angular momentum (kpc × km/s)
 *                           - Component 3: original particle ID
 * @param npts [in] Number of particles in the simulation.
 * @param rng [in] GSL random number generator instance for all stochastic processes.
 * @param Nscatter_total_step [out] Pointer to a long long to accumulate total scattering events this timestep.
 */
void perform_sidm_scattering_serial(gsl_rng *rng, long long *Nscatter_total_step) {
    long long Nscatters_this_call = 0;
    int i;

    // Iterate through each particle as potential scatterer
    for (i = 0; i < npts - 1; i++) {
        int nscat = 10; // Consider 10 nearest neighbors as scattering candidates
        if (npts - 1 - i < nscat)
            nscat = npts - 1 - i; // Limit to available particles
        if (nscat <= 0) continue;

        double partialprobability[nscat + 1]; // Interaction rates for each candidate
        double probability_sum_term = 0.0;   // Total interaction rate sum

        // Construct 3D velocity vector for primary particle
        // Random azimuthal orientation for transverse velocity component
        double phii = 2.0 * PI * gsl_rng_uniform(rng);
        double Viperp = particles[2][i] / particles[0][i]; // v_perp = L/r
        threevector Vi = to_vector(Viperp * cos(phii), Viperp * sin(phii), particles[1][i]);

        // Calculate interaction rates with neighboring particles
        for (int m = 1; m <= nscat; m++) {
            int partner_idx = i + m;

            // Construct 3D velocity for scattering partner
            // Assumes fixed azimuthal alignment for partner particle
            double Vmperp = particles[2][partner_idx] / particles[0][partner_idx];
            threevector Vm = to_vector(Vmperp, 0.0, particles[1][partner_idx]);

            threevector Vrel_vec = vector_diff(Vi, Vm);
            double vrel_val = sqrt(dotproduct(Vrel_vec, Vrel_vec));

            // Calculate interaction rate: σ × v_rel
            partialprobability[m] = sigmatotal(vrel_val, npts) * vrel_val;
            probability_sum_term += partialprobability[m];
        }

        // Determine the outer radius of the shell containing these nscat neighbors
        // Option 1 (Consistent): Shell defined by the nscat-th summed neighbor
        // Option 2 (Current): Shell defined by the (nscat+1)-th particle
        int use_nscat_plus_1_for_shell = 1; // Set to 1 for current method, 0 for alternative
        int outer_shell_particle_idx_for_vol;

        if (use_nscat_plus_1_for_shell)
            outer_shell_particle_idx_for_vol = i + nscat + 1;
        else
            outer_shell_particle_idx_for_vol = i + nscat;

        // Ensure the chosen outer index is within bounds
        if (outer_shell_particle_idx_for_vol >= npts) {
            // If out of bounds, try to use the last available particle as the boundary
            if (i + nscat < npts)
                outer_shell_particle_idx_for_vol = i + nscat;
            else // No valid shell can be formed
                probability_sum_term = 0.0; // Force no scatter, skip probability calculation
        }

        double radius_diff = 0.0;
        // Calculate radius_diff only if there's a chance to scatter and a valid shell
        if (probability_sum_term > 1e-30 && (outer_shell_particle_idx_for_vol > i))
            radius_diff = particles[0][outer_shell_particle_idx_for_vol] - particles[0][i];
        else
            probability_sum_term = 0.0; // Ensure no scatter if shell is invalid

        double probability = 0.0;
        // Calculate scattering probability using shell volume approximation
        if (radius_diff > 1e-15 && particles[0][i] > 1e-15 && probability_sum_term > 1e-30)
            probability = probability_sum_term * (0.5) * dt / (4.0 * PI * sqr(particles[0][i]) * radius_diff);

        // Stochastic scattering determination
        if (gsl_rng_uniform(rng) < probability) {
            Nscatters_this_call++;
            int m_scatter = 1; // Default to first neighbor

            // Weighted selection among multiple neighbors
            if (nscat > 1 && probability_sum_term > 1e-15) {
                double cumulative_prob[nscat + 1];
                cumulative_prob[0] = 0.0;
                for (int k = 1; k <= nscat; k++)
                    cumulative_prob[k] = (k > 1 ? cumulative_prob[k - 1] : 0.0) + partialprobability[k] / probability_sum_term;
                if (nscat > 0) cumulative_prob[nscat] = 1.0;

                double random_select = gsl_rng_uniform(rng);
                m_scatter = 1;
                // Select partner based on cumulative probability distribution
                while (m_scatter < nscat && random_select > cumulative_prob[m_scatter])
                    m_scatter++;
            }

            int actual_partner_idx = i + m_scatter;
            if (actual_partner_idx >= npts) {
                Nscatters_this_call--;
                continue;
            }

            // Reconstruct velocities for selected scattering pair
            double Vmperp_scatter = particles[2][actual_partner_idx] / particles[0][actual_partner_idx];
            threevector Vm_scatter = to_vector(Vmperp_scatter, 0.0, particles[1][actual_partner_idx]);
            threevector Vrel_scatter_vec = vector_diff(Vi, Vm_scatter);
            double vrel_scatter_val = sqrt(dotproduct(Vrel_scatter_vec, Vrel_scatter_vec));

            if (vrel_scatter_val < 1e-15) {
                Nscatters_this_call--;
                continue;
            }

            // Generate isotropic scattering angles in center-of-mass frame
            double costheta = 2.0 * gsl_rng_uniform(rng) - 1.0;
            double sintheta = sqrt(fmax(0.0, 1.0 - costheta * costheta));
            double phif_scatter = 2.0 * PI * gsl_rng_uniform(rng);
            double cf = cos(phif_scatter);
            double sf = sin(phif_scatter);

            // Construct orthonormal coordinate system for scattering transformation
            threevector nhat0, nhat1, nhat2, nhatref;
            nhat0 = vector_scalar(Vrel_scatter_vec, 1/vrel_scatter_val);

            nhatref = (fabs(nhat0.z) < 0.999) ? to_vector(0.0, 0.0, 1.0) : to_vector(1.0, 0.0, 0.0);

            nhat1 = crossproduct(nhat0, nhatref);
            double normnhat1 = sqrt(dotproduct(nhat1, nhat1));
            if (normnhat1 < 1e-15) {
                // Fallback for parallel vectors
                nhatref = (fabs(nhat0.z) < 0.999) ? to_vector(1.0, 0.0, 0.0) : to_vector(0.0, 1.0, 0.0);
                nhat1 = crossproduct(nhat0, nhatref);
                normnhat1 = sqrt(dotproduct(nhat1, nhat1));
                if (normnhat1 < 1e-15) {
                     Nscatters_this_call--;
                     continue;
                }
            }
            nhat1 = vector_scalar(nhat1, 1/normnhat1);
            nhat2 = crossproduct(nhat0, nhat1);

            // Transform scattered velocities from CM frame to lab frame
            threevector nhat_perp_rotated = vector_sum(vector_scalar(nhat1, cf), vector_scalar(nhat2, sf));
            threevector V_rel_final_half = vector_scalar(vector_sum(vector_scalar(nhat0,costheta), vector_scalar(nhat_perp_rotated,sintheta)),vrel_scatter_val/2.0);
            threevector V_cm = vector_scalar(vector_sum(Vi,Vm_scatter), 0.5);

            threevector Vifinal_vec = vector_sum(V_cm, V_rel_final_half);
            threevector Vmfinal_vec = vector_diff(V_cm, V_rel_final_half);

            // Apply velocity changes to particle data arrays
            particles[1][i] = Vifinal_vec.z;
            particles[1][actual_partner_idx] = Vmfinal_vec.z;

            double Vperp_i_final_mag = sqrt(sqr(Vifinal_vec.x) + sqr(Vifinal_vec.y));
            double Vperp_m_final_mag = sqrt(sqr(Vmfinal_vec.x) + sqr(Vmfinal_vec.y));

            particles[2][i] = particles[0][i] * Vperp_i_final_mag; // Update angular momentum
            particles[2][actual_partner_idx] = particles[0][actual_partner_idx] * Vperp_m_final_mag;

            // Mark both scattered particles in case it is needed elsewhere
            int orig_id1 = (int)particles[3][i];
            int orig_id2 = (int)particles[3][actual_partner_idx];
            if (g_particle_scatter_state != NULL) {
                if (orig_id1 >= 0 && orig_id1 < npts)
                    g_particle_scatter_state[orig_id1] = 1;
                if (orig_id2 >= 0 && orig_id2 < npts)
                    g_particle_scatter_state[orig_id2] = 1;
            }
        }
    }

    *Nscatter_total_step = Nscatters_this_call;
}

/**
 * @brief Performs SIDM scattering calculations for one timestep using OpenMP for parallelism.
 * @details This function implements a two-phase parallel algorithm for SIDM scattering:
 *          Phase 1 (Parallel Particle Evaluation):
 *            - The main particle loop (over `i`) is parallelized using OpenMP.
 *            - Each thread processes its assigned subset of primary particles (`i`).
 *            - For each particle `i`, it considers `nscat` neighbors (by sorted rank) as potential scattering partners.
 *            - Interaction rates and total scattering probability for particle `i` with its neighbors
 *              are calculated using a shell volume approximation for local density.
 *            - A per-thread GSL RNG (`local_rng` from `rng_per_thread_list`) is used for all
 *              stochastic decisions (scatter occurrence, partner selection, scattering angles).
 *            - If a scatter occurs for particle `i` with a chosen partner `i+m_scatter`, the
 *              resulting final 3D velocities for both particles are computed.
 *            - These outcomes (indices `i`, `m_offset`, and final velocities) are stored in a
 *              `ScatterEvent` structure and added to a dynamically resizing global buffer
 *              (`global_scatter_results`) under an OpenMP critical section to ensure thread-safe appending.
 *          Phase 2 (Serial Update from Buffered Results):
 *            - After the parallel loop completes, a single thread sorts the `global_scatter_results`
 *              (by primary particle index, then partner offset) to ensure deterministic application order.
 *            - It then iterates through the sorted scatter events and updates the main `particles`
 *              array (radial velocity `particles[1]` and angular momentum `particles[2]`) with the
 *              final post-scatter velocities.
 *            - Updates the `g_particle_scatter_state` flags for particles involved in scattering.
 *          The `particles` array is assumed to be sorted by radius prior to calling this function.
 *
 * @param particles             [in,out] Main particle data array: `particles[component][current_sorted_index]`.
 *                              Modified in-place with post-scattering velocities/angular momenta.
 * @param npts                  [in] Total number of particles.
 * @param rng_per_thread_list   [in] Array of GSL RNG states, one for each OpenMP thread.
 * @param num_threads_for_rng   [in] The number of allocated RNGs in `rng_per_thread_list` (should match max threads).
 * @param Nscatter_total_step   [out] Pointer to a long long to accumulate the total number of scatter events
 *                              that occurred in this timestep.
 */
void perform_sidm_scattering_parallel(gsl_rng **rng_per_thread_list, int num_threads_for_rng, long long *Nscatter_total_step) {
    long long Nscatters_this_call_atomic = 0; // Accumulated in parallel reduction

    // Buffer for storing scattering event outcomes from all threads
    ScatterEvent *global_scatter_results = NULL;
    size_t global_results_count = 0;
    size_t global_results_capacity = 0;
    // Initial capacity can be a small fraction of npts, e.g., npts/100 or a fixed moderate number
    // Adjust if typical scatter rates are known.
    size_t initial_capacity = (npts > 1000) ? (npts / 100) : 100;
    if (initial_capacity == 0)
        initial_capacity = 10; // Ensure non-zero for very small npts

    global_scatter_results = (ScatterEvent *)malloc(initial_capacity * sizeof(ScatterEvent));
    if (global_scatter_results == NULL) {
        fprintf(stderr, "Error: Failed to allocate initial global_scatter_results buffer.\n");
        // Don't CLEAN_EXIT here, try to proceed without SIDM for this step or log error
        *Nscatter_total_step = 0;
        return;
    }
    global_results_capacity = initial_capacity;

    #pragma omp parallel reduction(+:Nscatters_this_call_atomic)
    {
        gsl_rng *local_rng = NULL; // Initialize to NULL
        int thread_id_for_rng = 0;
        #ifdef _OPENMP
            thread_id_for_rng = omp_get_thread_num();
        #endif

        if (rng_per_thread_list != NULL && thread_id_for_rng < num_threads_for_rng && rng_per_thread_list[thread_id_for_rng] != NULL)
            local_rng = rng_per_thread_list[thread_id_for_rng];
        else {
            // Critical issue: Per-thread RNG not available for an active thread.
            // This should not happen if g_rng_per_thread is sized to omp_get_max_threads()
            // and num_threads_for_rng passed to this function matches that.
            // Proceeding with a shared g_rng would be unsafe and non-reproducible.
            // For now, this thread will not perform scattering.
            #pragma omp critical (rng_error_sidm_parallel)
            {
                fprintf(stderr, "CRITICAL SIDM WARNING: Thread %d has no valid per-thread RNG (num_threads_for_rng=%d). This thread will skip SIDM calculations.\n", thread_id_for_rng, num_threads_for_rng);
                log_message("ERROR", "CRITICAL SIDM: Thread %d missing per-thread RNG.", thread_id_for_rng);
            }
            // To make this thread skip its iterations of the omp for loop:
            // One way is to jump past the loop content for this thread.
            // A cleaner way is to check local_rng before using it inside the loop.
        }


        /**
         * Using schedule(static,1) to ensure deterministic assignment of particles
         * to threads, which is crucial for run-to-run reproducibility of the
         * parallel SIDM simulation when using per-thread RNGs seeded identically
         * across runs (for a fixed number of threads).
         */
        #pragma omp for schedule(static,1)
        for (int i = 0; i < npts - 1; i++) {
            // Check if this thread has a valid RNG before proceeding
            if (local_rng == NULL)
                continue; // This thread skips its assigned SIDM work

            int nscat = 10;
            if (npts - 1 - i < nscat)
                nscat = npts - 1 - i;
            if (nscat <= 0)
                continue;

            double partialprobability[nscat + 1]; // Max nscat=10, stack is fine
            double probability_sum_term = 0.0;

            double phii = 2.0 * PI * gsl_rng_uniform(local_rng);
            double Viperp = particles[2][i] / particles[0][i];
            threevector Vi = to_vector(Viperp * cos(phii), Viperp * sin(phii), particles[1][i]);

            for (int m = 1; m <= nscat; m++) {
                int partner_idx = i + m;
                double Vmperp = particles[2][partner_idx] / particles[0][partner_idx];
                threevector Vm = to_vector(Vmperp, 0.0, particles[1][partner_idx]);
                threevector Vrel_vec = vector_diff(Vi, Vm);
                double vrel_val = sqrt(dotproduct(Vrel_vec, Vrel_vec));
                partialprobability[m] = sigmatotal(vrel_val, npts) * vrel_val;
                probability_sum_term += partialprobability[m];
            }

            int use_nscat_plus_1_for_shell_par = 1; // Consistent with serial for now
            int outer_shell_particle_idx_for_vol_par;
            if (use_nscat_plus_1_for_shell_par)
                outer_shell_particle_idx_for_vol_par = i + nscat + 1;
            else
                outer_shell_particle_idx_for_vol_par = i + nscat;
            if (outer_shell_particle_idx_for_vol_par >= npts) {
                if (i + nscat < npts)
                    outer_shell_particle_idx_for_vol_par = i + nscat;
                else
                    probability_sum_term = 0.0;
            }

            double radius_diff_par = 0.0;
            if (probability_sum_term > 1e-30 && (outer_shell_particle_idx_for_vol_par > i) )
                radius_diff_par = particles[0][outer_shell_particle_idx_for_vol_par] - particles[0][i];
            else
                probability_sum_term = 0.0;

            double probability_par = 0.0;
            if (radius_diff_par > 1e-15 && particles[0][i] > 1e-15 && probability_sum_term > 1e-30)
                probability_par = probability_sum_term * (0.5) * dt / (4.0 * PI * sqr(particles[0][i]) * radius_diff_par);

            if (gsl_rng_uniform(local_rng) < probability_par) {
                Nscatters_this_call_atomic++; // Atomically increment shared counter
                int m_scatter = 1;
                if (nscat > 1 && probability_sum_term > 1e-15) {
                    double cumulative_prob[nscat + 1];
                    cumulative_prob[0] = 0.0;
                    for (int k_cs = 1; k_cs <= nscat; k_cs++)
                        cumulative_prob[k_cs] = (k_cs > 1 ? cumulative_prob[k_cs - 1] : 0.0) + partialprobability[k_cs] / probability_sum_term;
                    if (nscat > 0)
                        cumulative_prob[nscat] = 1.0;
                    double random_select = gsl_rng_uniform(local_rng);
                    while (m_scatter < nscat && random_select > cumulative_prob[m_scatter])
                        m_scatter++;
                }
                int actual_partner_idx = i + m_scatter;
                if (actual_partner_idx >= npts)
                    continue; // Should be rare with nscat logic

                double Vmperp_scatter = particles[2][actual_partner_idx] / particles[0][actual_partner_idx];
                threevector Vm_scatter = to_vector(Vmperp_scatter, 0.0, particles[1][actual_partner_idx]);
                threevector Vrel_scatter_vec = vector_diff(Vi, Vm_scatter);
                double vrel_scatter_val = sqrt(dotproduct(Vrel_scatter_vec, Vrel_scatter_vec));
                if (vrel_scatter_val < 1e-15)
                    continue;

                double costheta = 2.0 * gsl_rng_uniform(local_rng) - 1.0;
                double sintheta = sqrt(fmax(0.0, 1.0 - costheta*costheta));
                double phif_scatter = 2.0 * PI * gsl_rng_uniform(local_rng);
                double cf = cos(phif_scatter); double sf = sin(phif_scatter);
                threevector nhat0, nhat1, nhat2, nhatref; // Orthonormal basis construction (as in serial)
                nhat0 = vector_scalar(Vrel_scatter_vec, 1/vrel_scatter_val);
                nhatref = (fabs(nhat0.z) < 0.999) ? to_vector(0.0,0.0,1.0) : to_vector(1.0,0.0,0.0);
                nhat1 = crossproduct(nhat0,nhatref);
                double normnhat1 = sqrt(dotproduct(nhat1,nhat1));
                if (normnhat1 < 1e-15) {
                    nhatref = (fabs(nhat0.x) < 0.999) ? to_vector(1.0,0.0,0.0) : to_vector(0.0,1.0,0.0);
                    nhat1 = crossproduct(nhat0,nhatref);
                    normnhat1 = sqrt(dotproduct(nhat1,nhat1));
                    if (normnhat1 < 1e-15)
                        continue;
                }
                nhat1 = vector_scalar(nhat1, 1/normnhat1);
                nhat2 = crossproduct(nhat0,nhat1);

                threevector nhat_perp_rotated = vector_sum(vector_scalar(nhat1, cf), vector_scalar(nhat2, sf));
                threevector V_rel_final_half = vector_scalar(vector_sum(vector_scalar(nhat0,costheta), vector_scalar(nhat_perp_rotated,sintheta)),vrel_scatter_val/2.0);
                threevector V_cm = vector_scalar(vector_sum(Vi,Vm_scatter),0.5);

                ScatterEvent current_event;
                current_event.i = i;
                current_event.m_offset = m_scatter; // Store offset, not absolute index
                current_event.Vifinal = vector_sum(V_cm,V_rel_final_half);
                current_event.Vmfinal = vector_diff(V_cm,V_rel_final_half);

                #pragma omp critical (add_scatter_result_sidm)
                {
                    if (global_results_count >= global_results_capacity) {
                        size_t new_capacity = (global_results_capacity == 0) ? initial_capacity : global_results_capacity * 2;
                         // Cap growth to avoid excessive memory if many scatters happen (unlikely but safe)
                        if (new_capacity > (size_t)npts && global_results_capacity < (size_t)npts)
                            new_capacity = (size_t)npts;

                        ScatterEvent *new_results_buffer = (ScatterEvent *)realloc(global_scatter_results, new_capacity * sizeof(ScatterEvent));
                        if (!new_results_buffer)
                            // This is a critical error if realloc fails.
                            // For now, we'll just stop adding results, but ideally, log and potentially terminate.
                             fprintf(stderr, "CRITICAL ERROR: Failed to reallocate global_scatter_results buffer in thread %d.\n", thread_id_for_rng);
                            // To prevent further issues, we could try to signal other threads or exit.
                            // This error means we are likely out of memory.
                        else {
                            global_scatter_results = new_results_buffer;
                            global_results_capacity = new_capacity;
                        }
                    }
                    // Only add if capacity is sufficient (realloc might have failed)
                    if (global_results_count < global_results_capacity)
                         global_scatter_results[global_results_count++] = current_event;
                } // end critical section
            } // end if scatter occurs
        } // end omp for loop over particles i
    } // end parallel region

    // Phase 2: Serial Update - Apply buffered scatter results
    // Sort the collected scatter events to ensure deterministic application order
    if (global_results_count > 1)
        qsort(global_scatter_results, global_results_count, sizeof(ScatterEvent), compare_scatter_events);

    // This part is done by a single thread after the parallel computation.
    for (size_t k = 0; k < global_results_count; k++) {
        int p_i = global_scatter_results[k].i;
        int p_m_offset = global_scatter_results[k].m_offset;
        int p_partner_idx = p_i + p_m_offset;

        // Redundant check, but good for safety, especially if realloc failed silently for some threads
        if (p_i < 0 || p_i >= npts || p_partner_idx < 0 || p_partner_idx >= npts || p_m_offset <= 0)
            // log_message("WARNING", "Skipping invalid scatter event from buffer: i=%d, partner_idx=%d, m_offset=%d", p_i, p_partner_idx, p_m_offset);
            continue;

        threevector Vifinal_upd = global_scatter_results[k].Vifinal;
        threevector Vmfinal_upd = global_scatter_results[k].Vmfinal;

        particles[1][p_i] = Vifinal_upd.z; // Update radial velocity for particle i
        particles[1][p_partner_idx] = Vmfinal_upd.z; // Update radial velocity for partner

        double Vperp_i_final_mag_upd = sqrt(sqr(Vifinal_upd.x) + sqr(Vifinal_upd.y));
        double Vperp_m_final_mag_upd = sqrt(sqr(Vmfinal_upd.x) + sqr(Vmfinal_upd.y));

        particles[2][p_i] = particles[0][p_i] * Vperp_i_final_mag_upd; // Update L for particle i
        particles[2][p_partner_idx] = particles[0][p_partner_idx] * Vperp_m_final_mag_upd; // Update L for partner

        // Mark both scattered particles in case it is needed elsewhere
        int orig_id1 = (int)particles[3][p_i];
        int orig_id2 = (int)particles[3][p_partner_idx];
        if (orig_id1 >= 0 && orig_id1 < npts)
            g_particle_scatter_state[orig_id1] = 1;
        if (orig_id2 >= 0 && orig_id2 < npts)
            g_particle_scatter_state[orig_id2] = 1;
    }

    free(global_scatter_results);
    *Nscatter_total_step = Nscatters_this_call_atomic;
}

/**
 * @brief Comparison function for `qsort` to order `ScatterEvent` structures.
 * @details Sorts an array of `ScatterEvent` structures primarily by the first particle's
 *          original index (`i`) in ascending order. If two events have the same primary
 *          particle index `i`, they are then secondarily sorted by the partner's offset
 *          (`m_offset`) in ascending order. This ensures a deterministic (and efficient
 *          for potential cache effects) order when applying buffered scatter updates to
 *          the main particle array, preventing race conditions or non-deterministic outcomes
 *          if multiple scatters involve the same primary particle.
 *
 * @param a [in] Pointer to the first `ScatterEvent` structure.
 * @param b [in] Pointer to the second `ScatterEvent` structure.
 * @return int - An integer less than, equal to, or greater than zero if the first
 *               argument is considered to be respectively less than, equal to,
 *               or greater than the second.
 */
int compare_scatter_events(const void *a, const void *b) {
    const ScatterEvent *event_a = (const ScatterEvent *)a;
    const ScatterEvent *event_b = (const ScatterEvent *)b;

    if (event_a->i < event_b->i)
        return -1;
    if (event_a->i > event_b->i)
        return 1;
    // If i is the same, sort by m_offset
    if (event_a->m_offset < event_b->m_offset)
        return -1;
    if (event_a->m_offset > event_b->m_offset)
        return 1;
    return 0;
}
