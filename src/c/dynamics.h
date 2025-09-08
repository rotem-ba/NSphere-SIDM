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

extern inline double gravitational_force(double r, int current_rank, int npts, double G_value, double halo_mass_value);
extern inline double effective_angular_force(double r, double ell);
extern inline double gravitational_force_rho_v(double rho, int current_rank, int npts, double G_value, double halo_mass_value);
extern inline double effective_angular_force_rho_v(double rho, double ell);

#endif // DYNAMICS_H
