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
 #include "logging.h"
 #include <math.h>
 #include <string.h>

 // =========================================================================
 //  Shared data variables
 // =========================================================================
 // Common data arrays
 double *radius = NULL;
 double *mass = NULL;
 double *Psivalues = NULL;
 double *nPsivalues = NULL;
 double *Evalues = NULL;
 double *innerintegrandvalues = NULL;
 double *radius_monotonic_grid_nfw = NULL;

 // Key scalar values
 double Psimin = 0.0;
 double Psimax = 0.0;
 double rmax = 0.0;
 int num_points = 0;

 // Common spline objects and accelerators
gsl_spline *splinemass = NULL;
gsl_interp_accel *enclosedmass = NULL;
gsl_spline *splinePsi = NULL;
gsl_interp_accel *Psiinterp = NULL;
gsl_spline *splinerofPsi = NULL;
gsl_interp_accel *rofPsiinterp = NULL;
gsl_interp *g_main_fofEinterp = NULL;
gsl_interp_accel *g_main_fofEacc = NULL;


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

 /**
  * @brief Determine active profile type (NFW is default).
  * @details Conforms the global tags to match the desired density.
  */
void set_global_density_params(){
    if (g_profile_type_str_provided) {
        if (strcmp(g_profile_type_str, "nfw") == 0) {
            g_use_nfw_profile = 1;
        } else if (strcmp(g_profile_type_str, "cored") == 0) {
            g_use_nfw_profile = 0;
        } else {
            // Should have been caught by parser, but as a safeguard:
            log_message("WARNING", "Unknown profile type '%s', defaulting to NFW.", g_profile_type_str);
            g_use_nfw_profile = 1;
        }
    } else {
        // Default to NFW if --profile flag was not provided
        g_use_nfw_profile = 1;
        strcpy(g_profile_type_str, "nfw"); // Update string for consistency in printouts
    }

    // Set up profile-specific parameters based on generalized flags and profile defaults
    if (g_use_nfw_profile) {
        // NFW Profile Path
        // Halo Mass for NFW
        if (g_halo_mass_param_provided) { // --halo-mass overrides NFW default
            g_nfw_profile_halo_mass = g_halo_mass_param;
        } else { // No --halo-mass, NFW uses its own default
            g_nfw_profile_halo_mass = HALO_MASS_NFW;
            g_halo_mass_param = g_nfw_profile_halo_mass; // Update general param to reflect NFW's choice
        }
        // Scale Radius for NFW
        if (g_scale_radius_param_provided) { // --scale-radius overrides NFW default
            g_nfw_profile_rc = g_scale_radius_param;
        } else { // No --scale-radius, NFW uses its own default
            g_nfw_profile_rc = RC_NFW_DEFAULT;
            g_scale_radius_param = g_nfw_profile_rc; // Update general param to reflect NFW's choice
        }
        // Cutoff Factor for NFW
        if (g_cutoff_factor_param_provided) { // --cutoff-factor overrides NFW default
            g_nfw_profile_rmax_norm_factor = g_cutoff_factor_param;
        } else { // No --cutoff-factor, NFW uses its own default
            g_nfw_profile_rmax_norm_factor = CUTOFF_FACTOR_NFW_DEFAULT;
            // g_cutoff_factor_param is NOT updated here by NFW default; it keeps its own (Cored's) default or user value.
        }
        // Falloff Factor for NFW
        if (g_falloff_factor_param_provided) { // --falloff-factor overrides NFW default
            g_nfw_profile_falloff_factor = g_falloff_factor_param;
        } else { // No --falloff-factor, NFW uses its own default
            g_nfw_profile_falloff_factor = FALLOFF_FACTOR_NFW_DEFAULT;
            // Optionally, update g_falloff_factor_param if NFW is the overall default and no flag given
            // For now, let g_falloff_factor_param keep its own default unless explicitly set by user
        }
    } else {
        // Cored Profile Path
        // Halo Mass for Cored (already defaults to HALO_MASS or takes from --halo-mass via g_halo_mass_param)
        g_cored_profile_halo_mass = g_halo_mass_param;
        // Scale Radius for Cored (already defaults to RC or takes from --scale-radius via g_scale_radius_param)
        g_cored_profile_rc = g_scale_radius_param;
        // Cutoff Factor for Cored (directly uses generalized or its (Cored's) default)
        g_cored_profile_rmax_factor = g_cutoff_factor_param;
    }

    // Set the single g_active_halo_mass for N-body forces and tdyn from the finalized g_halo_mass_param
    g_active_halo_mass = g_halo_mass_param;

}

/**
 * @def free spline and accelerator objects.
 */
void free_splines_accelerators() {
    gsl_spline_free(splinemass);
    gsl_spline_free(splinePsi);
    gsl_spline_free(splinerofPsi);
    gsl_interp_accel_free(enclosedmass);
    gsl_interp_accel_free(Psiinterp);
    gsl_interp_accel_free(rofPsiinterp);
    gsl_interp_free(g_main_fofEinterp);
    gsl_interp_accel_free(g_main_fofEacc);
}

/**
 * @def free density data arrays.
 */
void free_density_data_arrays() {
    free(mass);
    free(radius);
    if (radius_monotonic_grid_nfw != NULL) {
        free(radius_monotonic_grid_nfw);
        radius_monotonic_grid_nfw = NULL;
    }
    free(Psivalues);
    free(nPsivalues);
    free(innerintegrandvalues);
    free(Evalues);
}
