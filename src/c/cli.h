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
void errorAndExit(const char *msg, const char *arg, const char *prog);
int prompt_yes_no(const char *prompt);
void parseSaveArgs(int argc, char *argv[], int *pIndex);

#endif // CLI_H
