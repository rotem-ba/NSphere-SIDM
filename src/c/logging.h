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

#ifndef LOGGING_H
#define LOGGING_H

void log_message(const char *level, const char *format, ...);
void print_input_parameters();
void print_density_params();
void print_logging_status();
void print_SIDM_OpenMP_status();
void warn_no_OpenMP();
void log_scattering();
void log_disk_space(long long size);
void raise_insufficient_memory(long long total_disk_space, long long available_space);
void warn_low_memory(long long total_disk_space, long long available_space, double usage_after);

#endif // LOGGING_H
