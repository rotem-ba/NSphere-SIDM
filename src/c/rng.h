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

#ifndef RNG_H
#define RNG_H

#include <stdio.h>

void init_rng();
void init_rng_per_thread();

void set_sidm_seed(unsigned long current_time_pid_seed, char *seed_filepath, FILE *fp_seed);
void set_seed(unsigned long *seed, char *seed_type, int seed_provided, const char *seed_filename_base, unsigned long current_time_pid_seed);

#endif // RNG_H
