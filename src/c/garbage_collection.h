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

#ifndef GARBAGE_COLLECTION_H
#define GARBAGE_COLLECTION_H

#include <gsl/gsl_interp.h>

void free_accelerators(gsl_interp_accel **accelerators, int max_threads);
void free_double_array(double **array, int n);

#endif //GARBAGE_COLLECTION_H
