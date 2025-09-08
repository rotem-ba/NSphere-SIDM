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

 #ifndef DENSITY_H
 #define DENSITY_H

 #include <gsl/gsl_spline.h>
 #include <gsl/gsl_interp.h>

 // =========================================================================
 // ENERGY CALCULATION AND INTEGRATION STRUCTURES
 // =========================================================================

 /**
  * @brief Parameters for energy integration calculations.
  * @details Used as `void* params` argument in GSL integration routines,
  *          specifically for distribution function calculations (`fEintegrand`).
  */
  struct fEintegrand_params
  {
      double E;                      ///< Energy value (relative energy).
      gsl_spline *splinePsi;         ///< Interpolation spline for potential Psi(r).
      gsl_spline *splinemass;        ///< Interpolation spline for enclosed mass M(r).
      gsl_interp_accel *rofPsiarray; ///< Accelerator for radius lookup from potential r(Psi).
      gsl_interp_accel *massarray;   ///< Accelerator for mass lookups M(r).
  };
  typedef struct fEintegrand_params fEintegrand_params;

  /**
   * @brief Parameters for Psiintegrand to support profile-specific mass integrands.
   * @details Allows Psiintegrand to call the appropriate mass integrand function
   *          based on the selected density profile.
   */
  struct Psiintegrand_params{
      double (*massintegrand_func)(double, void *); ///< Function pointer to profile-specific mass integrand
      void *params_for_massintegrand;               ///< Parameters for the mass integrand function
  };
  typedef struct Psiintegrand_params Psiintegrand_params;

  double fEintegrand(double t, void *params);
  double Psiintegrand(double rp, void *params);

  double drhodr(double r);
  double massintegrand(double r, void *params __attribute__((unused)));

  void set_global_density_params();

#endif // DENSITY_H
