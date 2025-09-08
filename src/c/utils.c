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
