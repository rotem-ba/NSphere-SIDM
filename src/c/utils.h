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

/**
 * @brief Three-dimensional vector structure for SIDM calculations.
 * @details Used for representing velocity vectors and performing vector operations
 *          in SIDM scattering calculations.
 */
struct threevector {
    double x; ///< x-component of the vector.
    double y; ///< y-component of the vector.
    double z; ///< z-component of the vector.
};
typedef struct threevector threevector;

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
/** @def Constructs a three-dimensional vector from its Cartesian components */
#define to_vector(x, y, z) ((threevector){x, y, z})
/** @def Calculates the difference of 2 threevectors */
#define vector_diff(v1, v2) ((threevector){v1.x - v2.x, v1.y - v2.y, v1.z - v2.z})
/** @def Calculates the sum of 2 threevectors */
#define vector_sum(v1, v2) ((threevector){v1.x + v2.x, v1.y + v2.y, v1.z + v2.z})
/** @def Calculates the product of a threevector with a scalar */
#define vector_scalar(v1, scalar) ((threevector){v1.x * scalar, v1.y * scalar, v1.z * scalar})

/**
 * @brief Computes the scalar dot product of two three-dimensional vectors.
 * @details Calculates \f$X \cdot Y = X_x Y_x + X_y Y_y + X_z Y_z\f$.
 *          The dot product is a measure of the projection of one vector onto another
 *          and is used in various physics calculations, such as determining the
 *          magnitude squared of a vector (\f$V \cdot V = |V|^2\f$) or the angle between vectors.
 *
 * @param X [in] The first threevector operand.
 * @param Y [in] The second threevector operand.
 * @return double The scalar result of the dot product \f$X \cdot Y\f$.
 */
inline double dotproduct(threevector X, threevector Y) {
    return X.x * Y.x + X.y * Y.y + X.z * Y.z;
}

/**
 * @brief Computes the vector cross product of two three-dimensional vectors.
 * @details Calculates \f$Z = X \times Y\f$, where \f$X = (X_x, X_y, X_z)\f$ and \f$Y = (Y_x, Y_y, Y_z)\f$.
 *          The components of the resulting vector \f$Z\f$ are determined by:
 *          \f$Z_x = X_y Y_z - X_z Y_y\f$
 *          \f$Z_y = X_z Y_x - X_x Y_z\f$
 *          \f$Z_z = X_x Y_y - X_y Y_x\f$
 *          This follows the standard right-hand rule for vector cross products.
 *
 * @param X [in] The first threevector operand.
 * @param Y [in] The second threevector operand.
 * @return threevector The resulting vector \f$Z = X \times Y\f$.
 */
inline threevector crossproduct(threevector X, threevector Y) {
    threevector Z;
    Z.x = X.y * Y.z - X.z * Y.y;
    Z.y = X.z * Y.x - X.x * Y.z;
    Z.z = X.x * Y.y - X.y * Y.x;
    return Z;
}

int isInteger(const char *str);
int isFloat(const char *str);
double evaluatespline(gsl_spline *spline, gsl_interp_accel *acc, double value);
void fill_geomspace(double *r_grid, double *log_r_grid, int grid_size, double max_r, double min_r, int snap);

/** @def calculate total velocity norm from L, r, and vrad */
#define to_velocity(vrad,L,r) (sqrt(sqr(vrad) + sqr(L) / sqr(r)))

#endif // UTILS_H
