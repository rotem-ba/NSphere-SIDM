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
inline double gravitational_force(double r, int current_rank, int npts, double G_value, double halo_mass_value)
{
    if (use_identity_gravity)
    {
        // Testing mode: no gravitational force
        return 0.0;
    }
    else
    {
        // Calculate gravitational force: F = -G * M(r) / r²
        // where M(r) is proportional to particle rank
        return -(VEL_CONV_SQ * G_value) * ((double)current_rank / (double)npts) * halo_mass_value / (r * r);
    }
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
inline double effective_angular_force(double r, double ell)
{
    return (ell * ell) / (r * r * r);
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
    if (use_identity_gravity)
    {
        // Testing mode: no gravitational force
        return 0.0;
    }
    else
    {
        // Gravitational force in transformed coordinates
        return -(VEL_CONV_SQ * G_value) * ((double)current_rank / (double)npts) * halo_mass_value / (rho * rho);
    }
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
inline double effective_angular_force_rho_v(double rho, double ell)
{
    return (ell * ell) / (rho * rho * rho * rho);
}

#endif // DYNAMICS_H
