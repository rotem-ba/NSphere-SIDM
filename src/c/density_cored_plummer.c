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

 #include "density_cored_plummer.h"
 #include "globals.h"
 #include "utils.h"
 #include <math.h>


 /**
  * @brief Calculates \f$d\rho/dr\f$ for the Cored Plummer-like density profile.
  * @details The density profile is \f$\rho(r) \propto (1 + (r/RC)^2)^{-3}\f$.
  *          This function computes its analytical derivative with respect to \f$r\f$.
  *          It directly uses the `RC` macro for the scale radius.
  *
  * @param r [in] Radial coordinate (kpc) at which to evaluate the derivative.
  * @return double The value of \f$d\rho/dr\f$ at radius `r`.
  */
 double drhodr_cored_plummer(double r)
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
 double massintegrand_cored_plummer(double r, void *params __attribute__((unused)))
 {
     double startingprofile = 1.0 / cube((1.0 + sqr(r / g_cored_profile_rc)));
     return r * r * startingprofile;
 }
