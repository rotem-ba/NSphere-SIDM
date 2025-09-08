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

 #ifndef DEBUG_H
 #define DEBUG_H

 #define DEBUG_PARTICLE_ID 4    // Particle ID to track for debugging
 #define DEBUG_MAX_STEPS 100000 // Maximum number of debug energy snapshots

 /** @brief Arrays for tracking energy components through simulation for debugging. */
 extern double dbg_approxE[DEBUG_MAX_STEPS]; ///< Theoretical model energy (per unit mass).
 extern double dbg_dynE[DEBUG_MAX_STEPS];    ///< Actual dynamical energy (per unit mass).
 extern double dbg_kinE[DEBUG_MAX_STEPS];    ///< Kinetic energy component (per unit mass).
 extern double dbg_potE[DEBUG_MAX_STEPS];    ///< Potential energy component (per unit mass).
 extern double dbg_time[DEBUG_MAX_STEPS];    ///< Simulation time at each snapshot (Myr).
 extern double dbg_radius[DEBUG_MAX_STEPS];  ///< Particle radius at each snapshot (kpc).
 extern int dbg_count;                       ///< Number of debug snapshots recorded.

 void store_debug_approxE(int snapIndex, double E_value, double time_val);
 void store_debug_dynE_components(int snapIndex, double totalE, double kinE, double potE, double time_val, double radius_val);
 void finalize_debug_energy_output(void);

 #endif // DEBUG_H
