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

#include "gravitation_dynamics.h"
#include "globals.h"
#include "particle_array_ops.h"
#include "density.h"
#include "utils.h"
#include <math.h>

// =========================================================================
// ADAPTIVE FULL LEAPFROG STEP: r(n), v(n) --> r(n+1), v(n+1)
// =========================================================================
//
// Physics-based time integration method with adaptive step refinement:
// - Subdivide the time step h = ΔT in powers-of-2 "micro-steps"
// - Compare a (2N+1)-step "coarse" vs. a (4N+1)-step "fine" integration
// - If within tolerance, return one of {coarse, fine, rich (Richardson extrapolation)}
// - Otherwise, double N and repeat until convergence or max subdivision reached
/**
 * @brief Performs a leapfrog integration step using a fixed number of micro-steps.
 * @details This helper function implements the leapfrog (Kick-Drift-Kick) integration
 *          over a total time interval `h` by dividing it into a sequence of micro-steps.
 *          The number of micro-steps is `subSteps` (e.g., \f$2N+1\f$ for a "coarse" pass or
 *          \f$4N+1\f$ for a "fine" pass in an adaptive scheme, where \f$N\f$ is `N_subdivision_factor`).
 *          The micro-timestep sizes for kicks and drifts are adjusted based on `N_subdivision_factor`
 *          and whether it's a coarse or fine integration sequence.
 *          Sequence: Initial half-kick, \f$((\text{subSteps}-1)/2 - 1)\f$ full Drift-Kick pairs,
 *          a final full Drift, and a final half-Kick.
 *
 * @param i                   [in] Particle index (0 to npts-1), for rank in gravitational force.
 * @param npts                [in] Total number of particles.
 * @param r_in                [in] Input radial position (kpc) at the start of the total interval `h`.
 * @param v_in                [in] Input radial velocity (kpc/Myr) at the start of `h`.
 * @param ell                 [in] Angular momentum per unit mass (kpc^2/Myr).
 * @param N_subdivision_factor [in] Base subdivision factor \f$N\f$ used to determine micro-timestep sizes.
 * @param subSteps            [in] Total number of Kicks/Drifts (e.g., \f$2N+1\f$ or \f$4N+1\f$).
 * @param grav                [in] Gravitational constant G (simulation units).
 * @param r_out               [out] Pointer to store the output radial position (kpc) after time `h`.
 * @param v_out               [out] Pointer to store the output radial velocity (kpc/Myr) after time `h`.
 */
void doMicroLeapfrog(int i, int npts, double r_in, double v_in, double ell, int N, int subSteps, double grav, double *r_out, double *v_out) {
    // Initialize current state with input values
    double r_curr = r_in;
    double v_curr = v_in;

    // Set up timestep sizes based on subdivision level
    double halfKick, midStep;
    if (subSteps == (2 * N + 1)) {
        // Coarse integration (2N+1 substeps)
        halfKick = dt / (2.0 * N);
        midStep = dt / (1.0 * N);
    } else {
        // Fine integration (4N+1 substeps)
        halfKick = dt / (4.0 * N);
        midStep = dt / (2.0 * N);
    }

    // Initial half-kick (velocity update)
    double force = gravitational_force(r_curr, i, npts, grav, g_active_halo_mass);
    double dvdt = force + effective_angular_force(r_curr, ell);
    v_curr += halfKick * dvdt;

    // Middle pattern of drift-kick pairs
    int pairs = (subSteps - 1) / 2; // Total number of drift-kick pairs
    // Execute all but the last pair (last kick handled separately)
    for (int pp = 1; pp <= (pairs - 1); pp++) {
        // Drift: update position using current velocity
        r_curr += midStep * v_curr;

        // Kick: update velocity using forces at new position
        force = gravitational_force(r_curr, i, npts, grav, g_active_halo_mass);
        dvdt = force + effective_angular_force(r_curr, ell);
        v_curr += midStep * dvdt;
    }

    // Final full drift
    r_curr += midStep * v_curr;

    // Final half-kick
    force = gravitational_force(r_curr, i, npts, grav, g_active_halo_mass);
    dvdt = force + effective_angular_force(r_curr, ell);
    v_curr += halfKick * dvdt;

    *r_out = r_curr;
    *v_out = v_curr;
}

/**
 * @brief Performs an adaptive full leapfrog step with error control over a physical timestep `h`.
 * @details This function implements an adaptive leapfrog algorithm to advance a particle's
 *          state \f$(r, v_{rad})\f$ over a full physical timestep `h` (\f$\Delta T_{phys}\f$). It iteratively
 *          refines the integration by comparing a "coarse" integration (using \f$2N+1\f$
 *          micro-steps via `doMicroLeapfrog`) with a "fine" integration (using \f$4N+1\f$
 *          micro-steps). The subdivision factor \f$N\f$ starts at 1 and is doubled if the
 *          relative differences in final radius and velocity between coarse and fine passes
 *          exceed `radius_tol` and `velocity_tol`, respectively. This process repeats up
 *          to a maximum subdivision factor `max_subdiv`.
 *          The final state \f$(r_{out}, v_{out})\f$ for the step `h` is chosen based on `out_type`
 *          (coarse, fine, or Richardson extrapolation) once convergence is met or
 *          `max_subdiv` is reached.
 *
 * @param i             [in] Particle index (0 to npts-1), for rank in gravitational force calculation.
 * @param npts          [in] Total number of particles in the simulation.
 * @param r_in          [in] Initial radial position (kpc) at the start of the step `h`.
 * @param v_in          [in] Initial radial velocity (kpc/Myr) at the start of `h`.
 * @param ell           [in] Angular momentum per unit mass (kpc^2/Myr), conserved.
 * @param radius_tol    [in] Relative convergence tolerance for radius comparison.
 * @param velocity_tol  [in] Relative convergence tolerance for velocity comparison.
 * @param max_subdiv    [in] Maximum allowed subdivision factor \f$N\f$ for micro-steps.
 * @param grav          [in] Gravitational constant G (simulation units, e.g., G_CONST).
 * @param out_type      [in] Result selection mode for converged integration:
 *                         0 for coarse result, 1 for fine result, 2 for Richardson extrapolation.
 * @param r_out         [out] Pointer to store the final radial position (kpc) after time `h`.
 * @param v_out         [out] Pointer to store the final radial velocity (kpc/Myr) after time `h`.
 */
void doAdaptiveFullLeap(int i, int npts, double r_in, double v_in, double ell, double radius_tol, double velocity_tol, int max_subdiv,
                        double grav, int out_type, double *r_out, double *v_out ) {
    int N = 1; // Start with N=1 micro-steps.

    // Initialize coarse/fine results (first coarse pass).
    double r_coarse = r_in, v_coarse = v_in;
    double r_fine = r_in, v_fine = v_in;

    while (N <= max_subdiv) {
        // Coarse pass (2N+1 steps)
        doMicroLeapfrog( i, npts, r_in, v_in, ell, N, (2 * N + 1), grav, &r_coarse, &v_coarse);

        // Fine pass (4N+1 steps)
        doMicroLeapfrog( i, npts, r_in, v_in, ell, N, (4 * N + 1), grav, &r_fine, &v_fine);

        // Compare radius, velocity.
        double radius_diff = fabs(r_fine - r_coarse) / (fabs(r_fine) + 1.0e-30);
        double velocity_diff = fabs(v_fine - v_coarse) / (fabs(v_fine) + 1.0e-30);

        if ((radius_diff < radius_tol) && (velocity_diff < velocity_tol)) {
            // Tolerances met => output according to out_type.
            if (out_type == 0) {
                *r_out = r_coarse;
                *v_out = v_coarse;
            } else if (out_type == 1) {
                *r_out = r_fine;
                *v_out = v_fine;
            } else {
                // Richardson extrapolation for higher-order result
                *r_out = 4.0 * r_fine - 3.0 * r_coarse;
                *v_out = 4.0 * v_fine - 3.0 * v_coarse;
            }
            return; // Done.
        } else // Not converged, increase refinement
            N *= 2;
    }

    // If we exit the while(N <= max_subdiv) loop, it means we never converged.
    // Return the last available result based on out_type.
    if (out_type == 0) {
        *r_out = r_coarse;
        *v_out = v_coarse;
    } else if (out_type == 1) {
        *r_out = r_fine;
        *v_out = v_fine;
    } else {
        *r_out = 4.0 * r_fine - 3.0 * r_coarse;
        *v_out = 4.0 * v_fine - 3.0 * v_coarse;
    }
}

// =========================================================================
// LEVI-CIVITA REGULARIZATION
// =========================================================================
//
// Physics-based regularized time integration method:
// - Transforms coordinates (r -> ρ = √r) to handle close encounters
// - Uses fictitious time τ to integrate equations of motion
// - Maps back to physical coordinates and time after integration
// - Provides enhanced stability for high-eccentricity orbits

/**
 * @brief Performs integration of particle motion over a physical time `dt` using Levi-Civita regularization.
 * @details This function implements a leapfrog-like integration scheme in Levi-Civita
 *          regularized coordinates \f$(ρ, v_{rad})\f$ and fictitious time \f$τ\f$. The physical
 *          coordinates are \f$r = ρ^2\f$, and physical time \f$t_{phys}\f$ is related to \f$τ\f$ by \f$dt_{phys} = ρ^2 dτ\f$.
 *          The integration proceeds by taking variable \f$Δτ\f$ steps (estimated based on `N_taumin`)
 *          until the accumulated physical time `t_phys` reaches or exceeds the target physical
 *          timestep `dt`. If a step overshoots `dt`, linear interpolation is used to obtain
 *          the state precisely at the target physical time.
 *          This method is particularly effective for handling close encounters where \f$r → 0\f$.
 *
 * @param i         [in] Particle index (0 to npts-1), for rank in force calculation.
 * @param npts      [in] Total number of particles.
 * @param r_in      [in] Initial physical radial position (kpc) at the start of the physical step `dt`.
 * @param v_in      [in] Initial physical radial velocity (kpc/Myr) at the start of `dt`.
 * @param ell       [in] Angular momentum per unit mass (kpc^2/Myr) for the force calculation.
 * @param N_taumin  [in] Target number of fictitious \f$τ\f$-steps within the `dt` interval; influences
 *                     the initial guess for \f$Δτ\f$.
 * @param grav      [in] Gravitational constant G (simulation units).
 * @param r_out     [out] Pointer to store the final physical radial position (kpc) after time `dt`.
 * @param v_out     [out] Pointer to store the final physical radial velocity (kpc/Myr) after time `dt`.
 */
void doLeviCivitaLeapfrog(int i, int npts, double r_in, double v_in, double ell, int N_taumin, double grav, double *r_out, double *v_out) {
    // Transform to Levi-Civita coordinates: rho = sqrt(r)
    double rho = sqrt(r_in);
    double v_rad = v_in;
    double tau = 0.0;
    double t_cur = 0.0;

    // Estimate initial tau step size based on N_taumin
    double deltaTau = 0.0;
    if (r_in > 1.0e-30 && N_taumin > 0)
        deltaTau = dt / (2.0 * r_in * N_taumin);
    else
        deltaTau = dt / 100.0;

    // Initialize integration state variables
    double rho_cur = rho;
    double v_cur = v_rad;
    double tau_cur = tau;
    double t_phys = t_cur;

    int stepCount = 0;
    int stepMax = 200000000; // Maximum step count to prevent infinite loops.

    while (1) {
        if (t_phys >= dt) // Exit loop when physical time reaches target
            break;

        // Leapfrog step 1: Evaluate force at current position
        double fval = forceLCfun(i, npts, g_active_halo_mass, grav, ell, rho_cur);

        // Leapfrog step 2: First half-kick for velocity
        double v_half = v_cur + 0.5 * deltaTau * fval;

        // Leapfrog step 3: Full drift for position
        double rho_next = rho_cur + deltaTau * dRhoDtaufun(rho_cur, v_half);

        // Leapfrog step 4: Second half-kick with force at new position
        double fval2 = forceLCfun(i, npts, g_active_halo_mass, grav, ell, rho_next);
        double v_next = v_half + 0.5 * deltaTau * fval2;

        // Update physical time using midpoint rho value
        double rho_mid = 0.5 * (rho_cur + rho_next);
        double t_next = t_phys + deltaTau * (rho_mid * rho_mid);

        double tau_next = tau_cur + deltaTau;

        // Handle case where step overshoots target time
        if (t_next >= dt) {
            double alpha = 0.0;
            if (fabs(t_next - t_phys) > 1.0e-30)
                alpha = (dt - t_phys) / (t_next - t_phys);
            else
                alpha = 1.0;

            // Interpolate state to exact target time
            double rho_f = rho_cur + alpha * (rho_next - rho_cur);
            double v_f = v_cur + alpha * (v_next - v_cur);

            // Transform back to physical coordinates
            double r_f = rho_f * rho_f;
            *r_out = r_f;
            *v_out = v_f;
            return;
        }

        // Update state for next iteration
        rho_cur = rho_next;
        v_cur = v_next;
        tau_cur = tau_next;
        t_phys = t_next;

        stepCount++;
        if (stepCount > stepMax) {
            // Safety exit if maximum iteration count exceeded
            *r_out = rho_cur * rho_cur;
            *v_out = v_cur;
            return;
        }
    }

    // Transform final state back to physical coordinates
    *r_out = rho_cur * rho_cur;
    *v_out = v_cur;

    return;
}

// =========================================================================
// ADAPTIVE FULL LEVI-CIVITA REGULARIZATION
// =========================================================================
//
// Enhanced regularization scheme that combines adaptive step sizing with
// Levi-Civita coordinate transformation for optimal performance near
// the coordinate origin.

/**
 * @brief Performs Levi-Civita regularized integration using a fixed number of micro-steps.
 * @details This helper function advances the particle state in regularized coordinates
 *          \f$(ρ, v_{rad}, t_{phys})\f$ over a total fictitious time interval `h_tau` (\f$Δτ_{total}\f$)
 *          by taking a specified number of `subSteps` fixed-size micro-steps (\f$δτ = Δτ_{total} / \text{subSteps}\f$).
 *          Each micro-step uses a leapfrog-like scheme (Kick-Drift-Kick for \f$ρ, v_{rad}\f$)
 *          and updates the accumulated physical time \f$t_{phys}\f$ using \f$dt_{phys} = δτ \cdot ρ_{mid}^2\f$.
 *          This function is called by `doSingleTauStepAdaptiveLeviCivita` for its coarse and fine passes.
 *
 * @param i         [in] Particle index (0 to npts-1), for rank in force calculation.
 * @param npts      [in] Total number of particles.
 * @param rho_in    [in] Initial \f$ρ = \sqrt{r}\f$ at the start of the `h_tau` interval.
 * @param v_in      [in] Initial radial velocity \f$v_{rad}\f$ at the start of `h_tau`.
 * @param t_in      [in] Initial accumulated physical time \f$t_{phys}\f$ at the start of `h_tau`.
 * @param subSteps  [in] Number of fixed micro-steps to perform over `h_tau`.
 * @param h_tau     [in] Total fictitious time interval \f$Δτ_{total}\f$ for this integration sequence.
 * @param grav      [in] Gravitational constant G (simulation units).
 * @param ell       [in] Angular momentum per unit mass (kpc^2/Myr).
 * @param rho_out   [out] Pointer to store the final \f$ρ\f$ after `h_tau`.
 * @param v_out     [out] Pointer to store the final \f$v_{rad}\f$ after `h_tau`.
 * @param t_out     [out] Pointer to store the final accumulated \f$t_{phys}\f$ after `h_tau`.
 */
void doMicroLeviCivita(int i, int npts, double rho_in, double v_in, double t_in, int subSteps, double h_tau, double grav, double ell, double *rho_out,
                       double *v_out, double *t_out) {
    double dtau = h_tau / (double)subSteps;
    double rho_curr = rho_in;
    double v_curr = v_in;
    double t_curr = t_in;

    for (int ss = 0; ss < subSteps; ss++) {
        // Half-kick.
        double fLC = forceLCfun(i, npts, g_active_halo_mass, grav, ell, rho_curr);
        double v_half = v_curr + 0.5 * dtau * fLC;

        // Drift for rho.
        double rho_next = rho_curr + dtau * (0.5 * rho_curr * v_half);

        // Second half-kick.
        double fLC2 = forceLCfun(i, npts, g_active_halo_mass, grav, ell, rho_next);
        double v_next = v_half + 0.5 * dtau * fLC2;

        // Accumulate physical time t(τ).
        // Simplest approach: use midpoint for rho =>  ρ_mid^2.
        double rho_mid = 0.5 * (rho_curr + rho_next);
        double dt_phys = dtau * (rho_mid * rho_mid);

        t_curr += dt_phys;
        rho_curr = rho_next;
        v_curr = v_next;
    }

    *rho_out = rho_curr;
    *v_out = v_curr;
    *t_out = t_curr;
}

/**
 * @brief Performs a single adaptive step in Levi-Civita coordinates over a proposed fictitious time `h_guess`.
 * @details This function integrates the equations of motion in regularized \f$(ρ, v_{rad}, t_{phys})\f$
 *          coordinates over a proposed fictitious time interval `h_guess` (\f$Δτ_{guess}\f$).
 *          It employs an adaptive refinement strategy by comparing a "coarse" integration
 *          (using \f$2N+1\f$ micro-steps via `doMicroLeviCivita`) with a "fine" integration
 *          (using \f$4N+1\f$ micro-steps). The subdivision factor \f$N\f$ starts at 1 and is
 *          doubled if the relative differences in \f$ρ\f$ and \f$v_{rad}\f$ between coarse and fine
 *          results exceed `radius_tol` and `velocity_tol`, respectively. This continues
 *          up to `max_subdiv`. The final state for the \f$Δτ_{guess}\f$ step is chosen based on
 *          `out_type` (coarse, fine, or Richardson extrapolation).
 *          The function outputs the final \f$ρ_{out}\f$, \f$v_{out}\f$ (radial velocity), and
 *          accumulated physical time \f$t_{out}\f$ corresponding to this adaptive \f$Δτ\f$ step.
 *
 * @param i             [in] Particle index (0 to npts-1), for rank in force calculation.
 * @param npts          [in] Total number of particles.
 * @param rho_in        [in] Initial regularized radial coordinate \f$ρ = \sqrt{r}\f$ at the start of \f$Δτ_{guess}\f$.
 * @param v_in          [in] Initial radial velocity \f$v_{rad}\f$ at the start of \f$Δτ_{guess}\f$.
 * @param t_in          [in] Initial accumulated physical time \f$t_{phys}\f$ at the start of \f$Δτ_{guess}\f$.
 * @param h_guess       [in] The proposed total fictitious time step \f$Δτ_{guess}\f$ for this adaptive step.
 * @param radius_tol    [in] Relative convergence tolerance for \f$ρ\f$ comparison.
 * @param velocity_tol  [in] Relative convergence tolerance for \f$v_{rad}\f$ comparison.
 * @param max_subdiv    [in] Maximum subdivision factor \f$N\f$ for micro-steps within `doMicroLeviCivita`.
 * @param grav          [in] Gravitational constant G (simulation units).
 * @param ell           [in] Angular momentum per unit mass (kpc^2/Myr).
 * @param out_type      [in] Result selection for converged micro-integration: 0=coarse, 1=fine, 2=Richardson.
 * @param rho_out       [out] Pointer to store the final \f$ρ\f$ after the adaptive \f$Δτ_{guess}\f$ step.
 * @param v_out         [out] Pointer to store the final \f$v_{rad}\f$ after the adaptive \f$Δτ_{guess}\f$ step.
 * @param t_out         [out] Pointer to store the final accumulated physical time \f$t_{phys}\f$ after this step.
 */
void doSingleTauStepAdaptiveLeviCivita(int i, int npts, double rho_in, double v_in, double t_in, double h_guess, double radius_tol, double velocity_tol,
                                       int max_subdiv, double grav, double ell, int out_type, double *rho_out, double *v_out, double *t_out) {
    // Initialize with the smallest subdivision factor N=1
    int N = 1;

    while (N <= max_subdiv) {
        // COARSE integration: use subSteps = 2N+1
        double rhoC, vC, tC;
        doMicroLeviCivita(i, npts, rho_in, v_in, t_in, (2 * N + 1), h_guess, grav, ell, &rhoC, &vC, &tC);

        // FINE integration: use subSteps = 4N+1
        double rhoF, vF, tF;
        doMicroLeviCivita(i, npts, rho_in, v_in, t_in, (4 * N + 1), h_guess, grav, ell, &rhoF, &vF, &tF);

        // Compare final radius and velocity for convergence
        // Calculate relative differences between coarse and fine solutions

        double rF = rhoF * rhoF;

        double rhodif = fabs(rhoF - rhoC) / (fabs(rF) + 1.0e-30); // Relative radial difference
        double vdif = fabs(vF - vC) / (fabs(vF) + 1.0e-30);       // Relative velocity difference

        if ((rhodif < radius_tol) && (vdif < velocity_tol)) {
            if (out_type == 0) {
                *rho_out = rhoC;
                *v_out = vC;
                *t_out = tC;
            } else if (out_type == 1) {
                *rho_out = rhoF;
                *v_out = vF;
                *t_out = tF;
            } else {
                // Apply Richardson extrapolation formula: result = 4*fine - 3*coarse
                // This provides a higher-order approximation by eliminating leading error terms
                double rho_rich = 4.0 * rhoF - 3.0 * rhoC; // Extrapolated ρ value
                double v_rich = 4.0 * vF - 3.0 * vC;       // Extrapolated velocity
                double t_rich = 4.0 * tF - 3.0 * tC;       // Extrapolated time value

                // Ensure non-negative radius
                if (rho_rich < 0.0)
                    rho_rich = 0.0;
                *rho_out = rho_rich;
                *v_out = v_rich;
                *t_out = t_rich;
            }
            return;
        } else // Not converged, double subdivision factor and try again
            N *= 2;
    }

    // Convergence not achieved within max_subdiv iterations
    // Use highest-resolution fine integration as fallback result
    double rhoF, vF, tF;
    doMicroLeviCivita(i, npts, rho_in, v_in, t_in, (4 * N + 1), h_guess, grav, ell, &rhoF, &vF, &tF);

    *rho_out = rhoF;
    *v_out = vF;
    *t_out = tF;
}

/**
 * @brief Performs a full adaptive leapfrog step using Levi-Civita regularization over a physical time interval `dt`.
 * @details This function integrates a particle's motion over a physical timestep `dt` (\f$\Delta T_{phys}\f$)
 *          by taking multiple adaptive steps in fictitious Levi-Civita time \f$τ\f$.
 *          It repeatedly calls `doSingleTauStepAdaptiveLeviCivita` to advance the state in
 *          \f$(ρ, v_{rad}, t_{phys})\f$ coordinates, where \f$ρ = \sqrt{r}\f$. Each call to
 *          `doSingleTauStepAdaptiveLeviCivita` takes an adaptive \f$Δτ\f$ step.
 *          The loop continues until the accumulated physical time `t_cur` (from summing \f$Δt_{phys}\f$
 *          corresponding to each \f$Δτ\f$) reaches or exceeds the target `dt`.
 *          If a \f$Δτ\f$ step overshoots `dt`, linear interpolation is used to find the
 *          state precisely at `t_cur = dt`. The final regularized state \f$(ρ_f, v_f)\f$ is then
 *          transformed back to physical coordinates \f$(r_{out}, v_{out})\f$.
 *          An initial guess for \f$Δτ\f$ is made based on `N_taumin` and the initial radius.
 *
 * @param i             [in] Particle index (0 to npts-1), used for rank in gravitational force calculation.
 * @param npts          [in] Total number of particles in the simulation.
 * @param r_in          [in] Initial physical radial position (kpc) at the start of the physical step `dt`.
 * @param v_in          [in] Initial physical radial velocity (kpc/Myr) at the start of `dt`.
 * @param ell           [in] Angular momentum per unit mass (kpc^2/Myr), conserved during integration.
 * @param N_taumin      [in] Target number of fictitious \f$τ\f$-steps within `dt`; influences the initial \f$Δτ\f$ guess.
 * @param radius_tol    [in] Relative convergence tolerance for \f$ρ\f$ comparison within each adaptive \f$τ\f$-step.
 * @param velocity_tol  [in] Relative convergence tolerance for velocity comparison within each adaptive \f$τ\f$-step.
 * @param max_subdiv    [in] Maximum allowed subdivision factor N for micro-steps within each adaptive \f$τ\f$-step.
 * @param grav          [in] Gravitational constant G (simulation units, e.g., G_CONST).
 * @param out_type      [in] Result selection mode for micro-steps within `doSingleTauStepAdaptiveLeviCivita`:
 *                         0 for coarse, 1 for fine, 2 for Richardson extrapolation.
 * @param r_out         [out] Pointer to store the final physical radial position (kpc) after time `dt`.
 * @param v_out         [out] Pointer to store the final physical radial velocity (kpc/Myr) after time `dt`.
 */
void doAdaptiveFullLeviCivita(int i, int npts, double r_in, double v_in, double ell, int N_taumin, double radius_tol, double velocity_tol,
                              int max_subdiv, double grav, int out_type, double *r_out, double *v_out) {
    // Handle near-zero radius edge case
    if (r_in < 1.0e-30) {
        *r_out = r_in;
        *v_out = v_in;
        return;
    }

    // Transform from physical to regularized coordinates
    double rho_current = sqrt(r_in);
    double v_current = v_in;

    double t_cur = 0.0; // Current physical time.

    // Estimate initial fictitious time step deltaTau
    double deltaTau = 0.0;
    if (r_in > 1.0e-30 && N_taumin > 0)
        deltaTau = dt / (2.0 * r_in * N_taumin);
    else
        deltaTau = dt / 100.0;

    int stepCount = 0;
    int stepMax = 100000000; // Safety limit on iteration count

    while (1) {
        if (t_cur >= dt)
            break; // Exit when physical time target is reached

        double rho_next, v_next, t_next;
        doSingleTauStepAdaptiveLeviCivita(i, npts, rho_current, v_current, t_cur, deltaTau, radius_tol, velocity_tol, max_subdiv, grav, ell, out_type,
                                          &rho_next, &v_next, &t_next);

        if (t_next > dt) {
            // Handle overshoot case with linear interpolation
            double alpha = 0.0;
            if (fabs(t_next - t_cur) > 1.0e-30)
                alpha = (dt - t_cur) / (t_next - t_cur);

            // Interpolate to exact target time dt
            double rho_final = rho_current + alpha * (rho_next - rho_current);
            double v_final = v_current + alpha * (v_next - v_current);

            double r_fin = rho_final * rho_final;
            *r_out = r_fin;
            *v_out = v_final;
            return;
        } else {
            // Update state variables for next iteration
            rho_current = rho_next;
            v_current = v_next;
            t_cur = t_next;
        }

        // Check for iteration limit
        stepCount++;
        if (stepCount > stepMax) {
            // Return best available result if maximum iterations reached
            double r_fin = rho_current * rho_current;
            *r_out = r_fin;
            *v_out = v_current;
            return;
        }
    }

    // Transform final regularized state back to physical coordinates
    double r_fin = rho_current * rho_current;
    *r_out = r_fin;
    *v_out = v_current;
}

// TRACKING


void update_trajectory_tacking(int current_step, int *inverse_map, int upper_npts_num_traj, double **trajectories, double **energies,
                               double **mu_arr, double **L_arr, double **E_arr, double **velocities_arr) {
    #pragma omp parallel for if (upper_npts_num_traj > 1000) schedule(static)
    for (int p = 0; p < upper_npts_num_traj; p++) {
        int idx = inverse_map[p];
        double rr = particles[0][idx];
        double vrad = particles[1][idx];
        double ell = particles[2][idx];
        double Psi_val = evaluatespline(splinePsi, Psiinterp, rr) * VEL_CONV_SQ;
        double vtot = to_velocity(vrad,ell,rr);
        double E_rel = Psi_val - 0.5 * sqr(vtot);
        trajectories[p][current_step] = rr;
        energies[p][current_step] = E_rel;
        E_arr[p][current_step] = E_rel;
        velocities_arr[p][current_step] = vrad;
        mu_arr[p][current_step] = vrad / vtot;
        L_arr[p][current_step] = ell;
    }
}

void update_inverse_map(int *inverse_map){
    #pragma omp parallel for default(shared) schedule(static)
    for (int idx = 0; idx < npts; idx++)
        inverse_map[(int)particles[3][idx]] = idx;
}

/**
 * @brief Performs euler_step method.
 */
void euler_step() {
    #pragma omp single
    sort_particles(particles, npts);
    #pragma omp barrier

    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++) {
        double r = particles[0][i];
        double vrad = particles[1][i];
        double ell = particles[2][i];
        double drdt = vrad;

        double force = gravitational_force(r, i, npts, G_CONST, g_active_halo_mass);
        double dvdt = force + effective_angular_force(r, ell);

        particles[0][i] += drdt * dt;
        particles[1][i] += dvdt * dt;
    }
}

/**
 * @brief Performs leapfrog method (position half step).
 */
void leapfrog_method_position_half_step() {
    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++)
        particles[0][i] += particles[1][i] * (dt / 2.0);
    #pragma omp single
    sort_particles(particles, npts);
    #pragma omp barrier

    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++) {
        double r = particles[0][i];
        double ell = particles[2][i];

        double force = gravitational_force(r, i, npts, G_CONST, g_active_halo_mass);
        double dvdt = force + effective_angular_force(r, ell);

        particles[1][i] += dvdt * dt;
        particles[0][i] += particles[1][i] * (dt / 2.0);
    }
}

/**
 * @brief Performs leapfrog method (velocity half step).
 */
void leapfrog_method_velocity_half_step() {
    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++) {
        double r = particles[0][i];
        double vrad = particles[1][i];
        double ell = particles[2][i];

        double force = gravitational_force(r, i, npts, G_CONST, g_active_halo_mass);
        double dvdt = force + effective_angular_force(r, ell);

        particles[1][i] = vrad + 0.5 * dvdt * dt;
    }

    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++)
        particles[0][i] += particles[1][i] * dt;

    #pragma omp single
    sort_particles(particles, npts);
    #pragma omp barrier

    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++) {
        double r = particles[0][i];
        double vrad = particles[1][i];
        double ell = particles[2][i];

        double force = gravitational_force(r, i, npts, G_CONST, g_active_halo_mass);
        double dvdt = force + effective_angular_force(r, ell);

        particles[1][i] = vrad + 0.5 * dvdt * dt;
    }
}

/**
 * @brief Full-step adaptive leapfrog integration.
 * @details Performs single step integration from (r_n, v_n) to (r_{n+1}, v_{n+1})
 *          using adaptive timestep control and the doAdaptiveFullLeap function.
 */
void leapfrog_method_full_step_adaptive() {
    double velocity_tol = 1.0e-5;
    double radius_tol = 1.0e-5;
    int max_subdiv = 1;

    #pragma omp single
    sort_particles(particles, npts);
    #pragma omp barrier

    /**
     * @brief Main integration loop - adaptive leapfrog update for each particle.
     * @details Each particle is advanced independently using adaptive timestepping.
     */
    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++) {
        double r = particles[0][i];
        double v = particles[1][i];
        double ell = particles[2][i];

        // One full step => h = dt.
        doAdaptiveFullLeap(i, npts, r, v, ell, radius_tol, velocity_tol, max_subdiv, G_CONST, 0, &particles[0][i], &particles[1][i]);
    }
}

/**
 * @brief Hybrid integration with adaptive method selection.
 * @details Uses Levi-Civita regularization for close encounters (r < r_crit)
 *          and standard leapfrog otherwise. Radius threshold r_crit is
 *          dynamically calculated for each particle.
 */
void hybrid_adaptive_method(){
    double velocity_tol = 1.0e-8;
    double radius_tol = 1.0e-8;
    int max_subdiv = 4096 * 4096;
    int N_taumin = 1000;
    double alpha_param = 0.05;

    #pragma omp single
    sort_particles(particles, npts);
    #pragma omp barrier

    /**
     * @brief Main integration loop - method selection based on orbital parameters.
     * @details Dynamically selects between standard leapfrog and Levi-Civita
     *          regularization based on particle's radius and angular momentum.
     */
    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++) {
        double r = particles[0][i];
        double v = particles[1][i];
        double ell = particles[2][i];
        double r_new, v_new;
        // Define critical radius r_crit for switching integration method:
        // r_crit = alpha_param * (ell^2) / (G * M(r))
        double r_crit = 0.0;
        if (ell != 0.0) {
            double M_enc = ((double)i / (double)npts) * g_active_halo_mass;
            double gravPart = VEL_CONV_SQ * G_CONST * M_enc;
            r_crit = sqr(ell) * alpha_param / gravPart;
        } else
            r_crit = 0.0;

        if ((r > 1.0e-30) && (r < r_crit))
            doLeviCivitaLeapfrog(i, npts, r, v, ell, N_taumin, G_CONST, &r_new, &v_new);
        else
            doAdaptiveFullLeap(i, npts, r, v, ell, radius_tol, velocity_tol,max_subdiv, G_CONST, 2, &r_new, &v_new);

        particles[0][i] = r_new;
        particles[1][i] = v_new;
    }
}

/**
 * @brief Performs adaptive leapfrog with adaptive Levi-Civita.
 */
void adaptive_leapfrog_adaptive_levi_civita() {
    double velocity_tol = 1.0e-7;
    double radius_tol = 1.0e-7;
    int max_subdiv = 4096 * 4096 * 16;
    int out_type = 2;
    int N_taumin = 10;
    double alpha_param = 0.05;

    #pragma omp single
    sort_particles(particles, npts);
    #pragma omp barrier

    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++) {
        double r = particles[0][i];
        double v = particles[1][i];
        double ell = particles[2][i];

        // Calculate r_crit, using special handling for particle i=0 (M_enc=0) to avoid division by zero.
        double r_crit = 0.0;
        if (fabs(ell) > 1.0e-30) {
            double M_enc;
            if (i == 0) // Special handling for i=0 to avoid zero mass in denominator
                M_enc = 0.1 * (1.0 / (double)npts) * g_active_halo_mass;
            else
                M_enc = ((double)i / (double)npts) * g_active_halo_mass;
            double gravPart = (VEL_CONV_SQ * G_CONST) * M_enc;
            r_crit = (ell * ell) * alpha_param / gravPart;
        }

        double r_new, v_new;
        if (r > 1.0e-30 && r < r_crit) // Switch based on critical radius
            doAdaptiveFullLeviCivita(i, npts, r, v, ell, N_taumin, radius_tol, velocity_tol, max_subdiv, G_CONST, out_type, &r_new, &v_new);
        else
            doAdaptiveFullLeap(i, npts, r, v, ell, radius_tol, velocity_tol, max_subdiv, G_CONST, out_type, &r_new, &v_new);
        particles[0][i] = r_new;
        particles[1][i] = v_new;
    }
}

/**
 * @brief Performs 4th-order Forest-Ruth-Yoshida integrator.
 */
void forest_ruth_yoshida_integration() {
    // Coefficients for 4th-order Forest-Ruth-Yoshida integrator (c1=c3).
    // Derived from: c1 = 1 / (2 - 2^(1/3)), c2 = 1 - 2*c1
    double c1 = 0.6756035959798289;
    double c2 = -0.3512071919596578; // = 1.0 - 2.0 * c1
    double c3 = c1;

    /**
     * @brief STEP 1: Kick by (c1 * dt/2).
     * @details Velocity update using the old position.
     */
    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++) {
        double r = particles[0][i];
        double vrad = particles[1][i];
        double ell = particles[2][i];

        double force = gravitational_force(r, i, npts, G_CONST, g_active_halo_mass);
        double dvdt = force + effective_angular_force(r, ell);
        particles[1][i] = vrad + 0.5 * c1 * dt * dvdt;
    }

    /**
     * @brief STEP 2: Drift by (c1 * dt).
     * @details Position update using the intermediate velocity v^*.
     */
    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++) {
        particles[0][i] += particles[1][i] * c1 * dt;
    }

    #pragma omp single
    sort_particles(particles, npts);
    #pragma omp barrier

    /**
     * @brief STEP 3: Kick by ((c1 + c2) * dt/2).
     * @details Velocity update using the new position after the first drift.
     */
    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++) {
        double r = particles[0][i];
        double vrad = particles[1][i];
        double ell = particles[2][i];

        double force = gravitational_force(r, i, npts, G_CONST, g_active_halo_mass);
        double dvdt = force + effective_angular_force(r, ell);
        double coeff = 0.5 * (c1 + c2);
        particles[1][i] = vrad + coeff * dt * dvdt;
    }

    /**
     * @brief STEP 4: Drift by (c2 * dt).
     * @details Position update using the intermediate velocity.
     */
    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++) {
        particles[0][i] += particles[1][i] * c2 * dt;
    }

    /**
     * @brief STEP 5: Kick by ((c2 + c3) * dt/2).
     * @details Velocity update using the new position after the second drift.
     */
    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++) {
        double r = particles[0][i];
        double vrad = particles[1][i];
        double ell = particles[2][i];

        // Recompute acceleration at new position.
        double force = gravitational_force(r, i, npts, G_CONST, g_active_halo_mass);
        double dvdt = force + effective_angular_force(r, ell);

        // Combination of the remaining half of c2 and half of c3.
        double coeff = 0.5 * (c2 + c3);
        particles[1][i] = vrad + coeff * dt * dvdt;
    }

    /**
     * @brief STEP 6: Drift by (c3 * dt).
     * @details Final position update in this integration step.
     */
    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++)
        particles[0][i] += particles[1][i] * c3 * dt;

    /**
     * @brief STEP 7: Kick by (c3 * dt/2).
     * @details Final velocity update (half-kick) to complete the integration step.
     */
    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++) {
        double r = particles[0][i];
        double vrad = particles[1][i];
        double ell = particles[2][i];

        double force = gravitational_force(r, i, npts, G_CONST, g_active_halo_mass);
        double dvdt = force + effective_angular_force(r, ell);
        particles[1][i] = vrad + 0.5 * c3 * dt * dvdt;
    }
}

/**
 * @brief Performs classic RK4 method.
 */
void rk4_method() {
    double *r_orig_by_id = (double *)malloc(npts * sizeof(double));
    double *v_orig_by_id = (double *)malloc(npts * sizeof(double));
    double *k1r_by_id = (double *)malloc(npts * sizeof(double));
    double *k1v_by_id = (double *)malloc(npts * sizeof(double));
    double *k2r_by_id = (double *)malloc(npts * sizeof(double));
    double *k2v_by_id = (double *)malloc(npts * sizeof(double));
    double *k3r_by_id = (double *)malloc(npts * sizeof(double));
    double *k3v_by_id = (double *)malloc(npts * sizeof(double));
    double *k4r_by_id = (double *)malloc(npts * sizeof(double));
    double *k4v_by_id = (double *)malloc(npts * sizeof(double));
    double h = dt;

    // Store original state by orig_id.
    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++) {
        int orig_id = (int)particles[3][i];
        r_orig_by_id[orig_id] = particles[0][i];
        v_orig_by_id[orig_id] = particles[1][i];
    }

    // K1 calculation.
    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++) {
        int orig_id = (int)particles[3][i];
        double r = particles[0][i];
        double vrad = particles[1][i];
        double ell = particles[2][i];

        double drdt = vrad;
        // Use the gravitational_force and effective_angular_force functions.
        double force = gravitational_force(r, i, npts, G_CONST, g_active_halo_mass);
        double dvdt = force + effective_angular_force(r, ell);

        k1r_by_id[orig_id] = drdt;
        k1v_by_id[orig_id] = dvdt;
    }

    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++) {
        int orig_id = (int)particles[3][i];
        double r_mid = particles[0][i];
        double v_mid = particles[1][i];
        double ell = particles[2][i];

        double drdt = v_mid;
        // Use the gravitational_force and effective_angular_force functions.
        double force = gravitational_force(r_mid, i, npts, G_CONST, g_active_halo_mass);
        double dvdt = force + effective_angular_force(r_mid, ell);

        k2r_by_id[orig_id] = drdt;
        k2v_by_id[orig_id] = dvdt;
    }

    #pragma omp single
    sort_particles(particles, npts);
    #pragma omp barrier

    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++) {
        int orig_id = (int)particles[3][i];
        double r_mid = particles[0][i];
        double v_mid = particles[1][i];
        double ell = particles[2][i];

        double drdt = v_mid;
        // Use the gravitational_force and effective_angular_force functions.
        double force = gravitational_force(r_mid, i, npts, G_CONST, g_active_halo_mass);
        double dvdt = force + effective_angular_force(r_mid, ell);

        k3r_by_id[orig_id] = drdt;
        k3v_by_id[orig_id] = dvdt;
    }

    // K4 calculation.
    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++) {
        int orig_id = (int)particles[3][i];
        double r_end = particles[0][i];
        double v_end = particles[1][i];
        double ell = particles[2][i];

        double drdt = v_end;
        // Use the gravitational_force and effective_angular_force functions.
        double force = gravitational_force(r_end, i, npts, G_CONST, g_active_halo_mass);
        double dvdt = force + effective_angular_force(r_end, ell);

        k4r_by_id[orig_id] = drdt;
        k4v_by_id[orig_id] = dvdt;
    }

    // Final RK4 combination.
    #pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < npts; i++) {
        int orig_id = (int)particles[3][i];

        double r_new = r_orig_by_id[orig_id] + (h / 6.0) * (k1r_by_id[orig_id] + 2.0 * k2r_by_id[orig_id] + 2.0 * k3r_by_id[orig_id] + k4r_by_id[orig_id]);
        double v_new = v_orig_by_id[orig_id] + (h / 6.0) * (k1v_by_id[orig_id] + 2.0 * k2v_by_id[orig_id] + 2.0 * k3v_by_id[orig_id] + k4v_by_id[orig_id]);

        particles[0][i] = r_new;
        particles[1][i] = v_new;
    }

    // Free RK4 arrays.
    free(r_orig_by_id);
    free(v_orig_by_id);
    free(k1r_by_id);
    free(k1v_by_id);
    free(k2r_by_id);
    free(k2v_by_id);
    free(k3r_by_id);
    free(k3v_by_id);
    free(k4r_by_id);
    free(k4v_by_id);
}

/**
 * @brief Make a dynamic step using the selected method.
 */
void make_dynamic_step() {
    if (method_select == 1)
        adaptive_leapfrog_adaptive_levi_civita();
    else if (method_select == 2)
        hybrid_adaptive_method();
    else if (method_select == 3)
        leapfrog_method_full_step_adaptive();
    else if (method_select == 4)
        forest_ruth_yoshida_integration();
    else if (method_select == 6)
        leapfrog_method_velocity_half_step();
    else if (method_select == 7)
        leapfrog_method_position_half_step();
    else if (method_select == 8)
        rk4_method();
    else if (method_select == 9)
        euler_step();
}
