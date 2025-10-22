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


#ifndef TIMESTEPS_H
#define TIMESTEPS_H

#include "globals.h"
#include <math.h>

int adjust_ntimesteps(int Ntimes_initial, int nout, int dtwrite);
void adjust_Ntimes();

void initialize_simulation_time();

inline int get_print_step(int k) {return (int)floor(k * (double)(log_every_p()) / 100 * Ntimes);}

#endif // TIMESTEPS_H
