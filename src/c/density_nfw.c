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

 #include "density_nfw.h"
 #include "globals.h"
 #include <math.h>

 /**
  * @brief Calculates \f$d\rho/dr\f$ for the NFW-like profile with a power-law cutoff.
  * @details The NFW-like density profile used is:
  *          \f$\rho(r) = \text{nt_nfw_scaler} \times [ (r_s + \epsilon)(1+r_s)^2 (1 + (r_s/C)^N) ]^{-1}\f$
  *          where \f$r_s = r / \text{rc_param}\f$, \f$\epsilon\f$ is a softening parameter (0.01),
  *          \f$C\f$ is the `falloff_C_param`, and \f$N\f$ is a power-law index (10.0).
  *          This function computes the analytical derivative of this \f$\rho(r)\f$ with respect to \f$r\f$.
  *
  * @param r               [in] Radial coordinate (kpc) at which to evaluate the derivative.
  * @param rc_param        [in] Scale radius (RC) of the NFW-like profile (kpc).
  * @param nt_nfw_scaler   [in] Density normalization constant (nt_nfw) for the profile.
  * @param falloff_C_param [in] Falloff transition factor \f$C\f$ for the power-law cutoff.
  * @return double The value of \f$d\rho/dr\f$ at radius `r`. Returns 0.0 if `rc_param` is non-positive.
  */
 double drhodr_profile_nfwcutoff(double r, double rc_param, double nt_nfw_scaler, double falloff_C_param) { // Added falloff_C_param

     // Parameters for the NFW-like profile shape
     const double epsilon_softening = 0.01;
     double C_cutoff_factor = falloff_C_param;
     if (C_cutoff_factor <= 0) C_cutoff_factor = 19.0; // Safety default
     const double N_cutoff_power = 10.0;

     if (rc_param <= 0) return 0.0;

     double rs = r / rc_param;

     //rho_shape(rs) = 1.0 / ( (rs+eps) * (1+rs)^2 * (1+(rs/C)^N) )
     double term_s = rs + epsilon_softening;
     if (term_s <= 1e-9) term_s = 1e-9; // Avoid division by zero for denominator
     double term_n = (1.0 + rs) * (1.0 + rs);
     double term_c_base = rs / C_cutoff_factor;
     double term_c_pow_N = pow(term_c_base, N_cutoff_power);
     double term_c = 1.0 + term_c_pow_N;

     double density_shape_val;
     if (term_s < 1e-9 || term_n < 1e-9 || term_c < 1e-9) {
         density_shape_val = 0.0; // Or very large if r is tiny
         if (r < 1e-6 && term_s < 1e-3) {
             density_shape_val = 1.0 / (term_s * term_n * term_c);
         }
     } else {
         density_shape_val = 1.0 / (term_s * term_n * term_c);
     }

     // Derivative of each term in the denominator w.r.t rs (d/d(rs)):
     // d/drs (rs+eps) = 1
     // d/drs (1+rs)^2 = 2*(1+rs)
     // d/drs (1+(rs/C)^N) = N * (rs/C)^(N-1) * (1/C)

     double d_log_term_s_d_rs = 1.0 / term_s;
     double d_log_term_n_d_rs = 2.0 / (1.0 + rs);
     double d_log_term_c_d_rs = (N_cutoff_power / C_cutoff_factor) * pow(term_c_base, N_cutoff_power - 1.0) / term_c;
     if (!isfinite(term_c_base) || (term_c_base < 1e-9 && N_cutoff_power -1 < 0)) { // Avoid pow(small_negative_base)
          d_log_term_c_d_rs = 0; // If rs/C is zero and N-1 is negative
     }

     // d(rho_shape)/dr = d(rho_shape)/d(rs) * d(rs)/dr = d(rho_shape)/d(rs) * (1/rc_param)
     // d(log(rho_shape))/d(rs) = - (d_log_term_s_d_rs + d_log_term_n_d_rs + d_log_term_c_d_rs)
     // d(rho_shape)/d(rs) = rho_shape * d(log(rho_shape))/d(rs)

     double d_rho_shape_d_rs = -density_shape_val * (d_log_term_s_d_rs + d_log_term_n_d_rs + d_log_term_c_d_rs);
     double drho_dr = nt_nfw_scaler * d_rho_shape_d_rs / rc_param;


     return drho_dr;
 }

 /**
  * @brief Integrand for calculating the I(E) component of the NFW distribution function.
  * @details This function is integrated with respect to t_integration_var = sqrt(E_shell - Psi_true_at_r).
  *          It computes -2 * (d(rho)/d(Psi_true)) to ensure a positive integrand,
  *          leading to a monotonically increasing I(E). Psi_true_at_r is the potential
  *          corresponding to the radius r reached when energy E_shell has been reduced by t_integration_var^2.
  *          The radius r is determined from Psi_true_at_r using a spline.
  *          Derivatives d(rho)/dr and d(Psi_true)/dr are then calculated at this r.
  *
  * @param t_integration_var [in] The integration variable, t_prime = sqrt(E_shell - Psi_true_at_r).
  * @param params            [in] Void pointer to a `fE_integrand_params_NFW_t` structure. This struct
  *                               contains the current energy shell E_shell, splines for r(Psi) and M(r),
  *                               physical constants (G), NFW profile-specific parameters (scale radius,
  *                               density normalization, falloff factor C), and global Psimin/Psimax
  *                               for physical range validation.
  * @return double The value of the integrand -2 * (d(rho)/d(Psi_true)). Returns 0.0 if Psi_true_at_r
  *                is outside the physical range, if the derived radius is non-physical,
  *                if dPsi/dr is too small (or zero), or if the result is non-finite.
  */
 double fEintegrand_nfw(double t_integration_var, void *params) {
     fE_integrand_params_NFW_t *p_nfw = (fE_integrand_params_NFW_t *)params;

     double E_shell = p_nfw->E_current_shell; // Energy of the current shell for I(E)

     // Calculate Psi_true from t_integration_var: Psi_true = E_shell - t_integration_var^2
     // This means as t_integration_var goes from 0 to sqrt(E_shell - Psimin_global),
     // Psi_true goes from E_shell down to Psimin_global.
     double Psi_true_at_r = E_shell - t_integration_var * t_integration_var;


     // Ensure Psi_true_at_r is within the valid physical range [Psimin_global, Psimax_global]
     if (Psi_true_at_r < p_nfw->Psimin_global - 1e-7*fabs(p_nfw->Psimin_global) || Psi_true_at_r > p_nfw->Psimax_global + 1e-7*fabs(p_nfw->Psimax_global)) {
          return 0.0;
     }

     // Get radius r from Psi_true_at_r. Spline p->spline_r_of_Psi expects -Psi_true as input.
     double r_val;
     double spline_x_input_rPsi = -Psi_true_at_r; // Input for r_of_Psi spline
     double spline_rPsi_x_min = p_nfw->spline_r_of_Psi->x[0];
     double spline_rPsi_x_max = p_nfw->spline_r_of_Psi->x[p_nfw->spline_r_of_Psi->size - 1];

     if (spline_x_input_rPsi < spline_rPsi_x_min) spline_x_input_rPsi = spline_rPsi_x_min;
     if (spline_x_input_rPsi > spline_rPsi_x_max) spline_x_input_rPsi = spline_rPsi_x_max;

     r_val = gsl_spline_eval(p_nfw->spline_r_of_Psi, spline_x_input_rPsi, p_nfw->accel_r_of_Psi);


     // If radius is non-physical (negative or zero), the integrand is ill-defined or zero.
     if (r_val <= 1e-10) { // Using a slightly larger epsilon than machine precision for safety
         return 0.0;
     }
     // No more flooring of r_val here; use it as is if positive, or return 0 if not.

     // Calculate drho/dr at r_val
     double drho_dr_val = drhodr_profile_nfwcutoff(r_val, p_nfw->profile_rc_const, p_nfw->profile_nt_norm_const, p_nfw->profile_falloff_C_const);

     // Calculate dPsi_true/dr = G * M(r_val) / r_val^2 (magnitude)
     double M_at_r_val = gsl_spline_eval(p_nfw->spline_M_of_r, r_val, p_nfw->accel_M_of_r);
     if (M_at_r_val < 0) M_at_r_val = 0; // Mass must be non-negative

     // Calculate dPsi_true/dr = G * M(r_val) / r_val^2 (magnitude)
     // Since r_val > 1e-10, division by r_val^2 is safe
     double dPsi_dr_mag = p_nfw->const_G_universal * M_at_r_val / (r_val * r_val);

     if (fabs(dPsi_dr_mag) < 1e-30) { // If dPsi/dr is effectively zero (e.g. M(r)=0 at r=0)
         return 0.0; // drho/dPsi would be undefined or infinite
     }

     // drho/dPsi = (drho/dr) / (dPsi/dr)
     // Sign convention: Assume Psi is defined such that dPsi/dr is positive (potential less negative further out).
     // drho/dr is negative. So drho/dPsi is negative.
     // The quantity 2 * drho/dPsi is typically what appears in one form of Eddington's.
     double drho_dPsi_val = drho_dr_val / dPsi_dr_mag;
     // As per ANFIS.1 and subsequent findings, for I(E) to be increasing,
     // the integrand 2*d(rho)/d(Psi_true) needs to be positive.
     // Since drho_dPsi_val = (drho/dr) / (dPsi/dr_mag) is (negative/positive) = negative,
     // we need to flip the sign.
     double integrand_value = -2.0 * drho_dPsi_val;


     if (!isfinite(integrand_value)) {
         if (g_doDebug) fprintf(stderr, "Warning: NFW fEintegrand (refactored) returning non-finite value for t_in=%.3e, E_shell=%.3e\n", t_integration_var, E_shell);
         return 0.0; // Return 0 for non-finite cases
     }

     return integrand_value;
 }

 /**
  * @brief GSL integrand \f$r^2 \rho(r)\f$ for NFW-like profile mass calculation.
  * @details Computes \f$r^2 \rho(r)\f$ where \f$\rho(r)\f$ is an NFW-like density profile
  *          with an inner softening term and an outer power-law cutoff.
  *          The density \f$\rho(r) = p[2] \times [ (r_s + \epsilon)(1+r_s)^2 (1 + (r_s/C)^N) ]^{-1}\f$,
  *          where \f$r_s = r/p[0]\f$, \f$\epsilon=0.01\f$, \f$C=p[3]\f$, \f$N=10.0\f$.
  *          This integrand is used in GSL routines to calculate the normalization factor
  *          or enclosed mass for the NFW-like profile.
  *
  * @param r      [in] Radial coordinate \f$r\f$ (kpc).
  * @param params [in] Void pointer to a `double` array `p` of size 4:
  *                    - `p[0]` (rc_param): Scale radius (kpc).
  *                    - `p[1]` (halo_mass): Target total halo mass (Msun) - used to derive nt_nfw_scaler.
  *                    - `p[2]` (nt_nfw_scaler): Density normalization constant \f$nt_{NFW}\f$.
  *                    - `p[3]` (falloff_C_param): Falloff transition factor \f$C\f$.
  * @return double The value of the mass integrand \f$r^2 \rho(r)\f$. Returns 0.0 if `rc_param` (p[0]) is non-positive.
  */
 double massintegrand_profile_nfwcutoff(double r, void *params) {
     double *p = (double *)params;
     double rc_param = p[0];                 // Scale radius RC from parameters
     // p[1] is current_profile_halo_mass, not used directly in density formula
     double nt_nfw_scaler = p[2];            // Density scaling factor nt_nfw

     // Parameters for the NFW-like profile shape
     const double epsilon_softening = 0.01;  // Softening parameter for r/rc term
     double C_cutoff_factor = p[3];         // Falloff factor from params
     if (C_cutoff_factor <= 0) C_cutoff_factor = 19.0; // Safety default if param is bad
     const double N_cutoff_power = 10.0;     // Power for the cutoff term

     if (rc_param <= 0) { // Avoid division by zero if rc is invalid
         return 0.0;
     }

     double rs = r / rc_param; // r normalized by scale radius

     // Calculate the structural part of the density profile (unscaled by nt_nfw)
     double term_softening = rs + epsilon_softening;
     if (term_softening <= 1e-9) term_softening = 1e-9; // Avoid division by zero from softening

     double term_nfw_slope = (1.0 + rs) * (1.0 + rs); // (1 + r/rc)^2

     double cutoff_rs = rs / C_cutoff_factor;
     double term_cutoff = 1.0 + pow(cutoff_rs, N_cutoff_power);

     double density_shape;
     if (term_softening < 1e-9 || term_nfw_slope < 1e-9 || term_cutoff < 1e-9) { // Denominator terms too small
          density_shape = 0.0; // Or handle as very large if r is very small
          if (r < 1e-6 && term_softening < 1e-3) { // special handling for very small r to match NFW cusp
              density_shape = 1.0 / (term_softening * term_nfw_slope * term_cutoff); // Let it be large
          }
     } else {
          density_shape = 1.0 / (term_softening * term_nfw_slope * term_cutoff);
     }

     // Apply the overall density scaling factor
     double physical_density = nt_nfw_scaler * density_shape;

     return r * r * physical_density;
 }
