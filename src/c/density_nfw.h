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

 #ifndef DENSITY_NFW_H
 #define DENSITY_NFW_H

 #include "density.h"
 #include <gsl/gsl_spline.h>
 #include <gsl/gsl_interp.h>

 // =========================================================================
 // ENERGY CALCULATION AND INTEGRATION STRUCTURES
 // =========================================================================


 struct fE_integrand_params_NFW_t{
     double E_current_shell;         ///< Energy E of the current shell for which I(E) is being computed.
     gsl_spline *spline_r_of_Psi;    ///< Spline for r(-Psi_true), i.e., radius as a function of negated true potential.
     gsl_interp_accel *accel_r_of_Psi; ///< Accelerator for the r(-Psi_true) spline.
     gsl_spline *spline_M_of_r;      ///< Spline for M(r), enclosed mass as a function of radius.
     gsl_interp_accel *accel_M_of_r;   ///< Accelerator for the M(r) spline.
     double const_G_universal;       ///< Universal gravitational constant G.
     double profile_rc_const;        ///< Scale radius (rc) of the NFW profile.
     double profile_nt_norm_const;   ///< Density normalization constant (nt_nfw) for the NFW profile.
     double profile_falloff_C_const; ///< Falloff transition factor C for power-law cutoff in NFW profile.
     double Psimin_global;           ///< Minimum potential value for physical range validation.
     double Psimax_global;           ///< Maximum potential value for physical range validation.
 };
 typedef struct fE_integrand_params_NFW_t fE_integrand_params_NFW_t;

 double drhodr_profile_nfwcutoff(double r, double rc_param, double nt_nfw_scaler, double falloff_C_param);
 double fEintegrand_nfw(double t_integration_var, void *params);
 double massintegrand_profile_nfwcutoff(double r, void *params);

#endif // DENSITY_NFW_H
