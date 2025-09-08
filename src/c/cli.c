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

 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include "globals.h"
 #include "exit.h"

 // =========================================================================
 // COMMAND LINE ARGUMENT PROCESSING
 // =========================================================================

 /**
  * @brief Displays detailed usage information for command-line arguments.
  * @details Prints a comprehensive help message to stderr showing all available
  *          command-line options, their default values, and brief descriptions.
  *          Includes information about integration methods, sorting algorithms,
  *          data saving modes, and basic usage examples.
  *
  * Parameters
  * ----------
  * prog : const char*
  *     The program name to display in the usage message (typically argv[0]).
  *
  * Returns
  * -------
  * None (prints to stderr).
  *
  * @note Called when the user specifies `--help` or when errors occur during
  *       argument parsing.
  */

void printUsage(const char *prog)
 {
     fprintf(stderr,
             "Usage: %s [options]\n"
             "  --help                        Show this usage message.\n"
             "  --log                         [Default Off] Enable writing detailed logs to log/nsphere.log\n"
             "  --tag <string>                [Default Off] Add a custom tag to output filenames (e.g. \"run1\")\n"
             "  --save <subargs>              [Default all] Enable various data-saving modes (may combine any):\n"
             "                                          all           => output everything\n"
             "                                          raw-data      => only raw particle data\n"
             "                                          psi-snaps     => plus Psi snapshots\n"
             "                                          full-snaps    => plus full data snapshots\n"
             "                                          debug-energy  => plus energy diagnostics\n"
             "                                          If multiple subargs are given, the highest priority\n"
             "                                          one overrides the lower: raw-data < psi-snaps <\n"
             "                                          full-snaps < all/debug-energy.\n"
             "  --restart                     [Default Off] Restart processing from the last written snapshot\n"
             "\n"
             "  --nparticles <int>            [Default 100000] Number of particles\n"
             "  --ntimesteps <int>            [Default 10000] Requested total timesteps\n"
             "                                     Note: Ntimes will be adjusted to the minimum value that\n"
             "                                     satisfies the constraint (Ntimes - 1) = k*(dtwrite)*(nout)\n"
             "  --dtwrite <int>               [Default 100] Low level diskwrite interval in timesteps\n"
             "  --nout <int>                  [Default 100] Number of post-processing output data snapshot times\n"
             "  --tfinal <int>                [Default 5] Final simulation time in units of the dynamical time\n"
             "\n"
             "  --readinit <file>             [Default Off] Read initial conditions from <file> (binary)\n"
             "  --writeinit <file>            [Default Off] Write initial conditions to <file> (binary)\n"
             "  --master-seed <int>           [Default Random] Set a master seed to derive other seeds.\n"
             "  --load-seeds                  [Default Off] Load seeds from previous run's output files\n"
             "  --init-cond-seed <int>        [Default Random/Master] Set seed for IC generation.\n"
             "                                     Overrides derivation from master-seed.\n"
             "\n"
             "  --method <int>                [Default 1] Integration method (1..9):\n"
             "                                          1   Adaptive Leapfrog with Adaptive Levi-Civita\n"
             "                                          2   Full-step adaptive Leapfrog + Levi-Civita\n"
             "                                          3   Full-step adaptive Leapfrog\n"
             "                                          4   Yoshida 4th-order\n"
             "                                          5   Adams-Bashforth 3\n"
             "                                          6   Leapfrog (vel half-step)\n"
             "                                          7   Leapfrog (pos half-step)\n"
             "                                          8   Classic RK4\n"
             "                                          9   Euler\n"
             "  --methodtag                   [Default Off] Include method string in output filenames\n"
             "  --sort <int>                  [Default 1] Sorting algorithm (1..4):\n"
             "                                          1   Parallel Quadsort\n"
             "                                          2   Sequential Quadsort\n"
             "                                          3   Parallel Insertion Sort\n"
             "                                          4   Sequential Insertion Sort\n"
             "\n"
             "  --halo-mass <float>           [Default 1.15e9] Total halo mass in M☉ for the selected profile.\n"
             "  --profile <type>              [Default nfw] Profile type for ICs: 'nfw' or 'cored'.\n"
             "  --scale-radius <float>        [Default 1.18] Scale radius in kpc for the selected profile.\n"
             "  --cutoff-factor <float>       [Default 85.0] Absolute r_max in units of scale radius.\n"
             "  --falloff-factor <float>      [Default 19.0] NFW concentration parameter transition factor.\n"
             "  --ftidal <float>              [Default 0.0] Set tidal fraction outer stripping value (0.0 to 1.0)\n"
             "\n"
             "  --sidm                        [Default Off] Enable self-interacting scattering physics\n"
             "  --sidm-seed <int>             [Default Random/Master] Set seed for SIDM calculations.\n"
             "                                     Overrides derivation from master-seed.\n"
             "  --sidm-mode <serial|parallel> [Default parallel] Select SIDM execution mode.\n"
             "                                     Parallel mode requires OpenMP.\n"
             "  --sidm-kappa <float>          [Default 50.0] SIDM opacity kappa in cm^2/g.\n"
             "\n"
             "Example:\n"
             "  %s --nparticles 50000 --ntimesteps 20000 --tfinal 5 \\\n"
             "     --nout 100 --dtwrite 10 --method 6 --sort 2 --readinit initial_conditions.bin\n",
             prog, prog);
 }

 /**
  * @brief Displays an error message, suggests `--help`, and terminates the program.
  * @details Formats and prints an error message to stderr, includes a suggestion
  *          to use the `--help` flag for usage information, performs necessary
  *          cleanup of allocated resources, and then exits with a non-zero status.
  *
  * Parameters
  * ----------
  * msg : const char*
  *     The error message to display.
  * arg : const char*
  *     Optional argument value that caused the error (shown in quotes),
  *     or NULL to omit this part of the message.
  * prog : const char*
  *     The program name to display in the `--help` usage suggestion.
  *
  * Returns
  * -------
  * None (calls exit(1) and never returns).
  *
  * @note Calls cleanup_all_particle_data() before exiting to free allocated memory.
  * @warning This function does not return; execution is terminated.
  */
void errorAndExit(const char *msg, const char *arg, const char *prog)
 {
     if (arg != NULL)
     {
         fprintf(stderr, "Error: %s '%s'\n", msg, arg);
     }
     else
     {
         fprintf(stderr, "Error: %s\n", msg);
     }
     fprintf(stderr, "Use '%s --help' for usage information.\n", prog ? prog : "./nsphere");
     cleanup_all_particle_data();
     exit(1);
 }

 /**
  * @brief Prompts the user with a yes/no question and reads their response from stdin.
  * @details Displays the given `prompt` string followed by "[y/N]: ".
  *          Reads a line of input from the user.
  *          - Returns 1 (yes) if the first character of the input is 'y' or 'Y'.
  *          - Returns 0 (no) if the first character is 'n', 'N', or if the input is
  *            an empty line (user just pressed Enter, defaulting to No).
  *          - If any other input is received, the prompt is repeated.
  *          Handles potential EOF or read errors by defaulting to No.
  *
  * @param prompt [in] The question/prompt message to display to the user.
  * @return int 1 if the user confirms (yes), 0 otherwise (no/default).
  */
 int prompt_yes_no(const char *prompt) {
     int response;

     while (1) {
         printf("%s [y/N]: ", prompt);
         fflush(stdout);

         // Read entire line
         char buffer[256];
         if (fgets(buffer, sizeof(buffer), stdin) == NULL) {
             // EOF or error - treat as 'N'
             return 0;
         }

         // Check first character
         response = buffer[0];

         if (response == 'y' || response == 'Y') {
             return 1;
         } else if (response == 'n' || response == 'N' || response == '\n') {
             // Empty line (just Enter) or explicit 'n'/'N'
             return 0;
         }
         // Any other input repeats the prompt
     }
 }

 /**
  * @brief Parses sub-arguments for the `--save` command-line option.
  * @details This function is called when the `--save` option is encountered during
  *          command-line argument parsing. It reads subsequent arguments (until another
  *          option starting with '-' is found, or arguments end) which specify the
  *          level or type of data to save. It then sets the corresponding global data
  *          output flags (`g_doDebug`, `g_doDynPsi`, `g_doDynRank`, `g_doAllParticleData`)
  *          based on the highest priority valid sub-argument encountered.
  *          Valid sub-arguments and their priority (lowest to highest):
  *          - "raw-data": Enables `g_doAllParticleData`.
  *          - "psi-snaps": Enables `g_doAllParticleData`, `g_doDynPsi`.
  *          - "full-snaps": Enables `g_doAllParticleData`, `g_doDynPsi`, `g_doDynRank`.
  *          - "all" or "debug-energy": Enables all flags (`g_doDebug`, `g_doDynPsi`, `g_doDynRank`, `g_doAllParticleData`).
  *
  * @param argc   [in] The total argument count from `main()`.
  * @param argv   [in] The argument array from `main()`.
  * @param pIndex [in,out] Pointer to the current index in `argv`. On input, it points to the
  *                       `--save` option. On output, it is updated to point to the last
  *                       sub-argument consumed by this function.
  * @note Exits the program with an error message if an unknown sub-argument to `--save` is found.
  */
 void parseSaveArgs(int argc, char *argv[], int *pIndex)
 {
     // Use an integer priority to track the highest level of saving requested
     static int savePriority = 0; // 0=none, 1=raw, 2=psi, 3=full, 4=all/debug

     // Start checking arguments after "--save"
     int i = *pIndex + 1;

     // Process arguments until the end or another option (starting with '-') is found
     while (i < argc && argv[i][0] != '-')
     {
         const char *subarg = argv[i];

         if (strcmp(subarg, "all") == 0 || strcmp(subarg, "debug-energy") == 0)
         {
             savePriority = 4; // Highest priority
         }
         else if (strcmp(subarg, "full-snaps") == 0)
         {
             if (savePriority < 3)
                 savePriority = 3;
         }
         else if (strcmp(subarg, "psi-snaps") == 0)
         {
             if (savePriority < 2)
                 savePriority = 2;
         }
         else if (strcmp(subarg, "raw-data") == 0)
         {
             if (savePriority < 1)
                 savePriority = 1;
         }
         else
         {
             fprintf(stderr, "Error: unknown argument to --save '%s'\n", subarg);
             exit(1);
         }

         i++;
     }

     // Set global flags based on the highest priority encountered
     switch (savePriority)
     {
     case 4: // All or debug-energy
         g_doDebug = 1;
         g_doDynPsi = 1;
         g_doDynRank = 1;
         g_doAllParticleData = 1;
         break;
     case 3: // Full-snaps.
         g_doDebug = 0;
         g_doDynPsi = 1;
         g_doDynRank = 1;
         g_doAllParticleData = 1;
         break;
     case 2: // Psi-snaps.
         g_doDebug = 0;
         g_doDynPsi = 1;
         g_doDynRank = 0;
         g_doAllParticleData = 1;
         break;
     case 1: // Raw-data.
         g_doDebug = 0;
         g_doDynPsi = 0;
         g_doDynRank = 0;
         g_doAllParticleData = 1;
         break;
     default: // No saving option specified, all flags remain 0
         break;
     }

     // Update index to point to the last processed argument
     *pIndex = i - 1;
 }
