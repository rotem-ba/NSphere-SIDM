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
