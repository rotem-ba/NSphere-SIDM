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

#include <gsl/gsl_interp.h>

/** @def Garbage collection for accelerators */
void free_accelerators(gsl_interp_accel **accelerators, int max_threads) {
    for (int i = 0; i < max_threads; i++)
        gsl_interp_accel_free(accelerators[i]);
    free(accelerators);
}

/** @def Garbage collection for double arrays */
void free_double_array(double **array, int n) {
    for (int i = 0; i < n; i++)
        free(array[i]);
    free(array);
}
