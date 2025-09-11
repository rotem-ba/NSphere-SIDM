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

#ifndef CLI_H
#define CLI_H

void printUsage(const char *prog);
static int validate_int_input(int *i, int argc, char **argv, const char *param_name);
static float validate_float_input(int *i, int argc, char **argv, const char *param_name);
static void raise_incompatible_choice_error(const char *main_param, const char *param2, int condition);
int prompt_yes_no(const char *prompt);
void prompt_quit(const char *prompt, const char *exit_message);
void parseSaveArgs(int argc, char *argv[], int *pIndex);

void exit_on_help(int argc, char *argv[]);
void read_user_arguments(int argc, char *argv[]);
void parse_user_arguments(int argc, char *argv[]);
void set_method_name();


#endif // CLI_H
