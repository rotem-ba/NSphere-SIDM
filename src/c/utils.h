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

#ifndef UTILS_H
#define UTILS_H
#include <gsl/gsl_spline.h>

/** @def imin(a, b) Minimum of two integer values. */
#define imin(a, b) ((a) < (b) ? (a) : (b))
/** @def sqr(x) Calculates the square of a value. */
#define sqr(x) ((x) * (x))
/** @def cube(x) Calculates the cube of a value. */
#define cube(x) ((x) * (x) * (x))
/** @def in_range(x) Checks if the value x is in the range [gte,lte]. */
#define in_range(x,gte,lte) (x >= gte && x <= lte)
/** @def if x<low return low */
#define at_least(x,low) ((x < low) ? low : x)
/** @def if x>high return high */
#define at_most(x,high) ((x > high) ? high : x)
/** @def if x<low return low, if x>high, return high, else return x */
#define clip(x, low, high) (at_most(at_least(x,low),high))

int isInteger(const char *str);
int isFloat(const char *str);
double evaluatespline(gsl_spline *spline, gsl_interp_accel *acc, double value);
void fill_geomspace(double *r_grid, double *log_r_grid, int grid_size, double max_r, double min_r, int snap);

#endif // UTILS_H
