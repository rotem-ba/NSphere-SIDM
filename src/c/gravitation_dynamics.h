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

#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "globals.h"
#include "utils.h"

/**
 * @brief Calculates gravitational acceleration at a given radius.
 *
 * Parameters
 * ----------
 * r : double
 *     Radius in kpc.
 * current_rank : int
 *     Particle rank (0 to npts-1), used for M(r) approximation.
 * npts : int
 *     Total number of particles.
 * G_value : double
 *     Gravitational constant value (e.g., G_CONST).
 * halo_mass_value : double
 *     Total halo mass (e.g., HALO_MASS).
 *
 * Returns
 * -------
 * double
 *     Gravitational acceleration (force per unit mass) in simulation units (kpc/Myr^2).
 *
 * @note Returns 0.0 if `use_identity_gravity` is set to 1.
 * @note M(r) is approximated as `(current_rank / npts) * halo_mass_value`.
 */
inline double gravitational_force(double r, int current_rank, int npts, double G_value, double halo_mass_value) {
    if (use_identity_gravity) // Testing mode: no gravitational force
        return 0.0;
    else // Calculate gravitational force: F = -G * M(r) / r², where M(r) is proportional to particle rank
        return -(VEL_CONV_SQ * G_value) * ((double)current_rank / (double)npts) * halo_mass_value / sqr(r);
}

/**
 * @brief Computes the effective centrifugal acceleration due to angular momentum.
 *
 * Parameters
 * ----------
 * r : double
 *     Radius (kpc).
 * ell : double
 *     Angular momentum per unit mass (kpc^2/Myr).
 *
 * Returns
 * -------
 * double
 *     Centrifugal acceleration: L²/r³ (kpc/Myr^2).
 */
inline double effective_angular_force(double r, double ell) {
    return sqr(ell) / cube(r);
}

/**
 * @brief Alternative gravitational force calculation using transformed coordinates.
 * @details Used in the Levi-Civita regularization scheme. Calculates F/m in rho coordinates.
 *
 * Parameters
 * ----------
 * rho : double
 *     Transformed radial coordinate (sqrt(r)). Units: sqrt(kpc).
 * current_rank : int
 *     Particle rank (0 to npts-1).
 * npts : int
 *     Total number of particles.
 * G_value : double
 *     Gravitational constant value (e.g., G_CONST).
 * halo_mass_value : double
 *     Total halo mass (e.g., HALO_MASS).
 *
 * Returns
 * -------
 * double
 *     Gravitational acceleration in transformed coordinates (units related to kpc^(3/2)/Myr^2).
 *
 * @note Returns 0.0 if `use_identity_gravity` is set to 1.
 * @see gravitational_force
 * @see doLeviCivitaLeapfrog
 */
inline double gravitational_force_rho_v(double rho, int current_rank, int npts, double G_value, double halo_mass_value)
{
    if (use_identity_gravity) // Testing mode: no gravitational force
        return 0.0;
    else // Gravitational force in transformed coordinates
        return -(VEL_CONV_SQ * G_value) * ((double)current_rank / (double)npts) * halo_mass_value / sqr(rho);
}

/**
 * @brief Computes effective centrifugal acceleration in transformed coordinates.
 * @details Used in the Levi-Civita regularization scheme. Calculates L^2/r^3 in rho coordinates.
 *
 * Parameters
 * ----------
 * rho : double
 *     Transformed radial coordinate (sqrt(r)). Units: sqrt(kpc).
 * ell : double
 *     Angular momentum per unit mass (kpc^2/Myr).
 *
 * Returns
 * -------
 * double
 *     Centrifugal acceleration in transformed coordinates (units related to kpc^(3/2)/Myr^2).
 *
 * @see effective_angular_force
 * @see doLeviCivitaLeapfrog
 */
inline double effective_angular_force_rho_v(double rho, double ell) {
    return sqr(ell) / (rho * rho * rho * rho);
}

/**
 * @brief Calculates \f$d\rho/d\tau\f$, the derivative of the regularized coordinate \f$\rho\f$ with respect to fictitious time \f$\tau\f$.
 * @details In Levi-Civita regularization, \f$d\rho/d\tau = \frac{1}{2} \rho v_{rad}\f$, where \f$\rho = \sqrt{r}\f$
 *          and \f$v_{rad}\f$ is the radial velocity in physical units (though often represented as \f$v\f$ or \f$v_{\rho}\f$
 *          in transformed equations of motion depending on the specific formulation).
 *          This function implements this relationship.
 *
 * @param rhoVal [in] The current value of the regularized radial coordinate \f$\rho = \sqrt{r}\f$.
 * @param vVal   [in] The current radial velocity \f$v_{rad}\f$ (kpc/Myr).
 * @return double The value of \f$d\rho/d\tau\f$.
 */
inline double dRhoDtaufun(double rhoVal, double vVal) {
    // dρ/dτ = 0.5 * ρ * v
    return 0.5 * rhoVal * vVal;
}

/**
 * @brief Calculates the total effective force per unit mass in Levi-Civita transformed coordinates.
 * @details This function computes \f$F_{\rho}/m = (F_{grav,\rho} + F_{centrifugal,\rho})/m\f$,
 *          where \f$F_{grav,\rho}\f$ is the gravitational force and \f$F_{centrifugal,\rho}\f$ is the
 *          effective centrifugal force, both expressed in the regularized radial coordinate \f$\rho = \sqrt{r}\f$.
 *          It calls `gravitational_force_rho_v` and `effective_angular_force_rho_v`.
 *          This combined force is used in the equations of motion for Levi-Civita regularization.
 *
 * @param i          [in] Particle index (0 to npts-1), for rank in gravitational force calculation.
 * @param npts       [in] Total number of particles.
 * @param totalmass  [in] Total halo mass of the system (Msun) used for gravitational force.
 * @param grav       [in] Gravitational constant G (simulation units).
 * @param ell        [in] Angular momentum per unit mass (kpc^2/Myr).
 * @param rhoVal     [in] Current value of the regularized radial coordinate \f$\rho = \sqrt{r}\f$.
 * @return double    The total transformed force per unit mass \f$F_{\rho}/m\f$.
 */
inline double forceLCfun(int i, int npts, double totalmass, double grav, double ell, double rhoVal) {
    double gravPart = gravitational_force_rho_v(rhoVal, i, npts, grav, totalmass);
    double angPart = effective_angular_force_rho_v(rhoVal, ell);
    return gravPart + angPart;
}

void doMicroLeapfrog(int i, int npts,double r_in, double v_in, double ell, int N, int subSteps, double grav, double *r_out, double *v_out);
void doAdaptiveFullLeap(int i, int npts, double r_in, double v_in, double ell, double radius_tol, double velocity_tol, int max_subdiv,
                        double grav, int out_type, double *r_out, double *v_out);
void doLeviCivitaLeapfrog(int i, int npts, double r_in, double v_in, double ell, int N_taumin, double grav, double *r_out, double *v_out);
void doMicroLeviCivita(int i, int npts, double rho_in, double v_in, double t_in, int subSteps, double h_tau, double grav, double ell, double *rho_out,
                       double *v_out, double *t_out);
void doSingleTauStepAdaptiveLeviCivita(int i, int npts, double rho_in, double v_in, double t_in, double h_guess, double radius_tol, double velocity_tol,
                                       int max_subdiv, double grav, double ell, int out_type, double *rho_out, double *v_out, double *t_out);
void doAdaptiveFullLeviCivita(int i, int npts, double r_in, double v_in, double ell, int N_taumin, double radius_tol, double velocity_tol,
                              int max_subdiv, double grav, int out_type, double *r_out, double *v_out);

void update_trajectory_tacking(int current_step, int *inverse_map, int upper_npts_num_traj, double **trajectories, double **energies,
                               double **mu_arr, double **L_arr, double **E_arr, double **velocities_arr);
void update_inverse_map(int *inverse_map);
void euler_step();
void leapfrog_method_position_half_step();
void leapfrog_method_velocity_half_step();
void leapfrog_method_full_step_adaptive();
void hybrid_adaptive_method();
void adaptive_leapfrog_adaptive_levi_civita();
void forest_ruth_yoshida_integration();
void rk4_method();
void make_dynamic_step();
#endif // DYNAMICS_H
