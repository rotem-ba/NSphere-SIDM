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

 #include "density.h"
 #include "globals.h"
 #include "utils.h"
 #include "exit.h"
 #include <math.h>

 double fEintegrand(double t, void *params)
 {
     fEintegrand_params *p = (fEintegrand_params *)params;
     gsl_spline *splinePsi = p->splinePsi;
     gsl_spline *splinemass = p->splinemass;
     gsl_interp_accel *rofPsiarray = p->rofPsiarray;
     gsl_interp_accel *massarray = p->massarray;
     double E = p->E;

     double Psi = E - t * t;
     double r = evaluatespline(splinePsi, rofPsiarray, -Psi);
     double drhodr(double r1);
     double drhodpsi = -(g_cored_profile_halo_mass / normalization) * drhodr(r) / (G_CONST * evaluatespline(splinemass, massarray, r) / (r * r));
     return 2.0 * drhodpsi;
 }

 /**
  * @brief GSL integrand for calculating the gravitational potential \f$\Psi(r)\f$.
  * @details This function computes the integrand \f$M_{enc}(r')/r'\f$ or \f$r' \rho_{shape}(r')\f$
  *          (depending on the exact formulation of \f$\Psi\f$) needed for the integral part of
  *          the potential calculation: \f$\Psi(r) = G M(<r)/r + G \int_r^{\infty} \text{integrand_val } dr'\f$.
  *          It uses a function pointer (`p_psi->massintegrand_func`) passed via the `params`
  *          argument (a `Psiintegrand_params` struct) to call the profile-specific
  *          mass integrand (which itself returns \f$r'^2 \rho_{shape}(r')\f$ or similar).
  *          The function then typically divides by `rp` to get \f$r' \rho_{shape}(r')\f$ if `massintegrand_func`
  *          returned \f$r'^2 \rho_{shape}(r')\f$.
  *          It includes a safety check to exit if `rp` (radius prime) is non-positive.
  *
  * @param rp     [in] The radial integration variable \f$r'\f$ (kpc).
  * @param params [in] Void pointer to a `Psiintegrand_params` struct. This struct contains
  *                    a function pointer to the profile-specific mass integrand and its parameters.
  * @return double The value of the potential integrand at `rp`.
  *
  * @see Psiintegrand_params
  * @see massintegrand
  * @see massintegrand_profile_nfwcutoff
  */
 double Psiintegrand(double rp, void *params)
 {
     Psiintegrand_params *p_psi = (Psiintegrand_params *)params;
     if (rp <= 0.0)
     {
         printf("rp out of range\n");
         CLEAN_EXIT(1);
     }
     // Call the profile-specific mass integrand function
     return p_psi->massintegrand_func(rp, p_psi->params_for_massintegrand) / rp;
 }

 /**
  * @brief Calculates \f$d\rho/dr\f$ for the Cored Plummer-like density profile.
  * @details The density profile is \f$\rho(r) \propto (1 + (r/RC)^2)^{-3}\f$.
  *          This function computes its analytical derivative with respect to \f$r\f$.
  *          It directly uses the `RC` macro for the scale radius.
  *
  * @param r [in] Radial coordinate (kpc) at which to evaluate the derivative.
  * @return double The value of \f$d\rho/dr\f$ at radius `r`.
  */
 double drhodr(double r)
 {
     return -6.0 * r / (g_cored_profile_rc * g_cored_profile_rc) / pow(1.0 + sqr(r / g_cored_profile_rc), 4.0);
 }

 /**
  * @brief GSL integrand \f$r^2 \rho_{shape}(r)\f$ for Cored Plummer-like profile mass calculation.
  * @details Computes the term \f$r^2 \rho_{shape}(r)\f$ for the Cored Plummer-like density profile,
  *          where \f$\rho_{shape}(r) = (1 + (r/RC)^2)^{-3}\f$. This integrand is used in
  *          GSL numerical integration routines (e.g., `gsl_integration_qag`) to calculate
  *          the normalization factor or the enclosed mass \f$M(<r) = 4\pi \int_0^r r'^2 \rho_{physical}(r') dr'\f$.
  *          It directly uses the `RC` macro for the scale radius.
  *
  * @param r      [in] Radial coordinate \f$r\f$ (kpc).
  * @param params [in] Void pointer to parameters (unused in this version, hence `__attribute__((unused))`).
  * @return double The value of the mass integrand \f$r^2 \rho_{shape}(r)\f$ at radius `r`.
  */
 double massintegrand(double r, void *params __attribute__((unused)))
 {
     double startingprofile = 1.0 / cube((1.0 + sqr(r / g_cored_profile_rc)));
     return r * r * startingprofile;
 }
