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

#include <ctype.h>
#include <gsl/gsl_spline.h>
#include <stdio.h>

/**
 * @brief Checks if a given string represents a valid integer.
 * @details Allows an optional leading '+' or '-' sign. Validates that all
 *          subsequent characters in the string are digits. Returns 0 (false)
 *          for empty strings, strings containing only a sign, or strings with
 *          non-digit characters after the optional sign.
 *
 * @param str [in] The null-terminated string to check.
 * @return int 1 if the string is a valid integer, 0 otherwise.
 */
int isInteger(const char *str)
{
    if (*str == '-' || *str == '+')
    {
        str++;
    }
    if (!*str)
    {
        return 0;
    }
    while (*str)
    {
        if (!isdigit((unsigned char)*str))
        {
            return 0;
        }
        str++;
    }
    return 1;
}

/**
 * @brief Checks if a given string represents a valid floating-point number.
 * @details Validates if the input string conforms to common floating-point number
 *          formats, including an optional leading sign ('+' or '-'), digits,
 *          at most one decimal point (if not in exponent part), and an optional
 *          exponent part (e.g., "e+10", "E-5").
 *          The function requires at least one digit to be present for a number to be
 *          considered valid (e.g., "." or "+." are not valid floats).
 *
 * @param str [in] The null-terminated string to check.
 * @return int 1 if the string is a valid float, 0 otherwise.
 */
int isFloat(const char *str)
{
    if (*str == '-' || *str == '+')
    {
        str++;
    }
    if (!*str)
    {
        return 0;
    }

    int has_digit = 0;
    int has_decimal = 0;
    int has_exponent = 0;

    while (*str)
    {
        if (isdigit((unsigned char)*str))
        {
            has_digit = 1;
        }
        else if (*str == '.' && !has_decimal && !has_exponent)
        {
            has_decimal = 1;
        }
        else if ((*str == 'e' || *str == 'E') && !has_exponent && has_digit)
        {
            has_exponent = 1;
            str++;
            // Check for optional sign after exponent
            if (*str == '-' || *str == '+')
            {
                str++;
            }
            if (!*str || !isdigit((unsigned char)*str))
            {
                return 0; // Exponent must have at least one digit
            }
            // Don't reset has_digit - we already have valid digits before exponent
        }
        else
        {
            return 0;
        }
        str++;
    }

    return has_digit;
}

/**
 * @brief Safely evaluates a GSL spline at a given value with robust bounds checking.
 * @details This function evaluates the provided GSL spline at the specified `value`.
 *          It includes critical safety checks:
 *          1. It verifies that the `spline` and accelerator `acc` pointers are not NULL.
 *          2. It checks if the `value` is outside the defined range of the spline's x-values.
 *             If `value` is out of bounds, it clamps `value` to the nearest valid boundary
 *             (plus/minus a small MARGIN) before evaluation to prevent GSL domain errors.
 *          This robust approach ensures that spline evaluations do not cause crashes due to
 *          out-of-range inputs, which can occur due to floating-point inaccuracies or
 *          unexpected data.
 *
 * @param spline [in] Pointer to the initialized GSL spline object.
 * @param acc    [in] Pointer to the GSL interpolation accelerator associated with the spline.
 * @param value  [in] The x-coordinate at which to evaluate the spline.
 * @return double The interpolated y-value from the spline. Returns the boundary spline value
 *                if `value` was clamped. Returns 0.0 if `spline` or `acc` is NULL (error logged).
 */
double evaluatespline(gsl_spline *spline, gsl_interp_accel *acc, double value)
{
    // NULL pointer safety check.
    if (spline == NULL || acc == NULL)
    {
        fprintf(stderr, "Error: NULL pointer passed to evaluatespline (spline=%p, acc=%p)\n",
                (void *)spline, (void *)acc);
        return 0.0; // Return a default value instead of crashing.
    }

    // Get the actual min and max ranges of the spline from its data directly.
    double x_min = spline->x[0];
    double x_max = spline->x[spline->size - 1];

    // Ensure the value is within the valid interpolation range with a small safety margin.
    const double MARGIN = 1e-10; // Small safety margin.

    if (value < x_min)
    {
#ifdef DEBUG_SPLINE
        fprintf(stderr, "Warning: Spline interpolation value %g below minimum %g, clamping\n",
                value, x_min);
#endif
        // Clamp to minimum with a tiny margin to stay inside the valid range.
        return gsl_spline_eval(spline, x_min + MARGIN, acc);
    }
    else if (value > x_max)
    {
#ifdef DEBUG_SPLINE
        fprintf(stderr, "Warning: Spline interpolation value %g above maximum %g, clamping\n",
                value, x_max);
#endif
        // Clamp to maximum with a tiny margin to stay inside the valid range.
        return gsl_spline_eval(spline, x_max - MARGIN, acc);
    }

    // Normal case - value is within range.
    return gsl_spline_eval(spline, value, acc);
}
