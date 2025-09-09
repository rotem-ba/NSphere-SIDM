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
 #include "utils.h"
 #include "logging.h"

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
  * @brief Prompts the user with a yes/no question, and quits on a "yes".
  * @param prompt [in] The question/prompt message to display to the user.
  * @param exit_message [in] The exit message displayd if the user elects to quit.
  */
 void prompt_quit(const char *prompt, const char *exit_message) {
     if (!prompt_yes_no(prompt)) {
         printf("%s", exit_message);
         CLEAN_EXIT(0);
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

 /**
  * @brief Checks if the user used the "help" prompt.
  * @details If one of the keywords is the "help" prompt, print the help menu and exit.
  *
  * @param argc [in] The number of input arguments.
  * @param argv [in] The input arguments.
  * @return int 1 if the used the "help" prompt, 0 otherwise.
  */
int exit_on_help(int argc, char *argv[]) {
     for (int i = 1; i < argc; i++)
     {
         if (strcmp(argv[i], "--help") == 0)
         {
             printUsage(argv[0]);
             return 1;
         }
     }
     return 0;
 }


 /**
  * @brief Parse user arguments.
  * @details Reads out the user input arguments and update global fields.
  *
  * @param argc [in] The number of input arguments.
  * @param argv [in] The input arguments.
  * @return None
  */
 void read_user_arguments(int argc, char *argv[]) {
     /** @note Handle the case where no command-line arguments are provided. */
     if (argc == 1)
     {
         printf("No command-line arguments given. Using default parameters.\n");
         printf("Run `%s --help` to learn how to adjust parameters.\n\n", argv[0]);
     }

     /** @note Main command-line argument parsing loop (strict parsing). */
     for (int i = 1; i < argc; i++)
     {
         /** @note Forbid using '=' within options; require space separation. */
         // Check for "--option=value" format, which is disallowed.
         if (strncmp(argv[i], "--", 2) == 0 && strstr(argv[i], "=") != NULL)
         {
             errorAndExit("use space, not '=' after option", argv[i], argv[0]);
         }

         if (strcmp(argv[i], "--nparticles") == 0)
         {
             if (i + 1 >= argc)
             {
                 errorAndExit("--nparticles requires an integer argument", NULL, argv[0]);
             }
             if (!isInteger(argv[i + 1]))
             {
                 errorAndExit("invalid integer for --nparticles", argv[i + 1], argv[0]);
             }
             npts = atoi(argv[++i]);
         }
         else if (strcmp(argv[i], "--ntimesteps") == 0)
         {
             if (i + 1 >= argc)
             {
                 errorAndExit("--ntimesteps requires an integer argument", NULL, argv[0]);
             }
             if (!isInteger(argv[i + 1]))
             {
                 errorAndExit("invalid integer for --ntimesteps", argv[i + 1], argv[0]);
             }
             Ntimes = atoi(argv[++i]);
         }
         else if (strcmp(argv[i], "--tfinal") == 0)
         {
             if (i + 1 >= argc)
             {
                 errorAndExit("--tfinal requires an integer argument", NULL, argv[0]);
             }
             if (!isInteger(argv[i + 1]))
             {
                 errorAndExit("invalid integer for --tfinal", argv[i + 1], argv[0]);
             }
             tfinal_factor = atoi(argv[++i]);
         }
         else if (strcmp(argv[i], "--nout") == 0)
         {
             if (i + 1 >= argc)
             {
                 errorAndExit("--nout requires an integer argument", NULL, argv[0]);
             }
             if (!isInteger(argv[i + 1]))
             {
                 errorAndExit("invalid integer for --nout", argv[i + 1], argv[0]);
             }
             nout = atoi(argv[++i]);
         }
         else if (strcmp(argv[i], "--dtwrite") == 0)
         {
             if (i + 1 >= argc)
             {
                 errorAndExit("--dtwrite requires an integer argument", NULL, argv[0]);
             }
             if (!isInteger(argv[i + 1]))
             {
                 errorAndExit("invalid integer for --dtwrite", argv[i + 1], argv[0]);
             }
             dtwrite = atoi(argv[++i]);
         }
         else if (strcmp(argv[i], "--tag") == 0)
         {
             if (i + 1 >= argc)
             {
                 errorAndExit("--tag requires a string argument", NULL, argv[0]);
             }

             strncpy(custom_tag, argv[++i], 255);
             custom_tag[255] = '\0'; // Ensure null-termination.
         }
         else if (strcmp(argv[i], "--method") == 0)
         {
             if (i + 1 >= argc)
             {
                 errorAndExit("--method requires an integer argument", NULL, argv[0]);
             }
             if (!isInteger(argv[i + 1]))
             {
                 errorAndExit("invalid integer for --method", argv[i + 1], argv[0]);
             }
             method_select = atoi(argv[++i]);

             if (method_select < 1 || method_select > 9)
             {
                 char buf[256];
                 snprintf(buf, sizeof(buf), "method must be in [1..9]");
                 errorAndExit(buf, NULL, argv[0]);
             }
         }
         else if (strcmp(argv[i], "--sort") == 0)
         {
             if (i + 1 >= argc)
             {
                 errorAndExit("--sort requires an integer argument", NULL, argv[0]);
             }
             if (!isInteger(argv[i + 1]))
             {
                 errorAndExit("invalid integer for --sort", argv[i + 1], argv[0]);
             }
             int sort_val = atoi(argv[++i]);

             if (sort_val < 1 || sort_val > 4)
             {
                 char buf[256];
                 snprintf(buf, sizeof(buf), "sort must be in [1..4]");
                 errorAndExit(buf, NULL, argv[0]);
             }

             display_sort = sort_val;

             switch (sort_val)
             {
             case 1:
                 g_defaultSortAlg = "quadsort_parallel"; // Formerly case 1.
                 break;
             case 2:
                 g_defaultSortAlg = "quadsort"; // Formerly case 0.
                 break;
             case 3:
                 g_defaultSortAlg = "insertion_parallel"; // Formerly case 2.
                 break;
             case 4:
                 g_defaultSortAlg = "insertion"; // Formerly case 3.
                 break;
             }
         }
         else if (strcmp(argv[i], "--readinit") == 0)
         {
             if (i + 1 >= argc) // Check if filename argument exists
             {
                 errorAndExit("--readinit requires a file argument", NULL, argv[0]);
             }

             /** @warning Check for incompatibility with `--restart` and `--writeinit`. */
             if (g_doRestart)
             {
                 errorAndExit("--readinit is incompatible with --restart. These options cannot be used together. Use either --restart OR --readinit, not both.", NULL, argv[0]);
             }
             if (doWriteInit)
             {
                 errorAndExit("--readinit is incompatible with --writeinit. These options cannot be used together. Use either --readinit OR --writeinit, not both.", NULL, argv[0]);
             }

             const char* user_filename = argv[++i]; // consume next arg
             // Prefix the path with "init/" directory
             static char prefixed_read_path[512]; // Static buffer for the path
             snprintf(prefixed_read_path, sizeof(prefixed_read_path), "init/%s", user_filename);
             readInitFilename = prefixed_read_path; // Assign the prefixed path
             doReadInit = 1;
         }
         else if (strcmp(argv[i], "--writeinit") == 0)
         {
             if (i + 1 >= argc) // Check if filename argument exists
             {
                 errorAndExit("--writeinit requires a file argument", NULL, argv[0]);
             }

             /** @warning Check for incompatibility with `--restart` and `--readinit`. */
             if (g_doRestart)
             {
                 errorAndExit("--writeinit is incompatible with --restart. These options cannot be used together. Use either --restart OR --writeinit, not both.", NULL, argv[0]);
             }
             if (doReadInit)
             {
                 errorAndExit("--writeinit is incompatible with --readinit. These options cannot be used together. Use either --writeinit OR --readinit, not both.", NULL, argv[0]);
             }

             const char* user_filename = argv[++i]; // consume next arg
             // Prefix the path with "init/" directory
             static char prefixed_write_path[512]; // Static buffer for the path
             snprintf(prefixed_write_path, sizeof(prefixed_write_path), "init/%s", user_filename);
             writeInitFilename = prefixed_write_path; // Assign the prefixed path
             doWriteInit = 1;
         }
         else if (strcmp(argv[i], "--restart") == 0)
         {
             /** @warning Check for incompatibility with `--readinit` and `--writeinit`. */
             if (doReadInit)
             {
                 errorAndExit("--restart is incompatible with --readinit. These options cannot be used together. Use either --restart OR --readinit, not both.", NULL, argv[0]);
             }
             if (doWriteInit)
             {
                 errorAndExit("--restart is incompatible with --writeinit. These options cannot be used together. Use either --restart OR --writeinit, not both.", NULL, argv[0]);
             }

             g_doRestart = 1;
             printf("Restart mode enabled. Will look for existing data products to resume processing.\n\n");
         }
         else if (strcmp(argv[i], "--save") == 0)
         {
             /** @note Delegate parsing of subsequent arguments to parseSaveArgs. */
             parseSaveArgs(argc, argv, &i);
             // The main loop continues from the updated index 'i'.
         }
         else if (strcmp(argv[i], "--ftidal") == 0)
         {
             if (i + 1 >= argc)
             {
                 errorAndExit("--ftidal requires a float argument", NULL, argv[0]);
             }
             if (!isFloat(argv[i + 1]))
             {
                 errorAndExit("invalid float for --ftidal", argv[i + 1], argv[0]);
             }
             tidal_fraction = atof(argv[++i]);

             /** @warning Check range [0.0, 1.0]. */
             if (tidal_fraction < 0.0 || tidal_fraction > 1.0)
             {
                 char buf[256];
                 snprintf(buf, sizeof(buf), "tidal_fraction must be in [0.0..1.0]");
                 errorAndExit(buf, NULL, argv[0]);
             }
         }
         else if (strcmp(argv[i], "--methodtag") == 0)
         {
             /** @note Flag to include integration method string in output filename suffix. */
             include_method_in_suffix = 1;
         }
         else if (strcmp(argv[i], "--log") == 0)
         {
             /** @note Flag to enable logging to log/nsphere.log. */
             g_enable_logging = 1;
         }
         else if (strcmp(argv[i], "--sidm") == 0)
         {
             /** @note Flag to enable Self-Interacting Dark Matter physics. */
             g_enable_sidm_scattering = 1;
             // This flag does not take a value, so 'i' is not incremented further.
         }
         else if (strcmp(argv[i], "--sidm-mode") == 0)
         {
             if (i + 1 >= argc) {
                 errorAndExit("--sidm-mode requires an argument (serial or parallel)", NULL, argv[0]);
             }
             char* mode_arg = argv[++i];
             if (strcmp(mode_arg, "serial") == 0) {
                 g_sidm_execution_mode = 0;
             } else if (strcmp(mode_arg, "parallel") == 0) {
                 #ifndef _OPENMP
                     printf("Warning: OpenMP is not enabled in this build. SIDM will run serially despite '--sidm-mode parallel'.\n");
                     log_message("WARNING", "OpenMP not enabled, SIDM forced to serial despite --sidm-mode parallel request.");
                     g_sidm_execution_mode = 0; // Force serial if no OpenMP
                 #else
                     g_sidm_execution_mode = 1;
                 #endif
             } else {
                 errorAndExit("Invalid argument for --sidm-mode. Use 'serial' or 'parallel'.", mode_arg, argv[0]);
             }
         }
         else if (strcmp(argv[i], "--sidm-kappa") == 0) {
             if (i + 1 >= argc || !isFloat(argv[i + 1])) {
                 errorAndExit("--sidm-kappa requires a float argument", argv[i + 1], argv[0]);
             }
             g_sidm_kappa = atof(argv[++i]);
             if (g_sidm_kappa < 0) { // Kappa can be 0 (no interaction) but not negative
                 errorAndExit("--sidm-kappa must be non-negative", NULL, argv[0]);
             }
             g_sidm_kappa_provided = 1;
         }
         else if (strcmp(argv[i], "--master-seed") == 0) {
             if (i + 1 >= argc || !isInteger(argv[i + 1])) {
                 errorAndExit("--master-seed requires an integer argument", argv[i + 1], argv[0]);
             }
             g_master_seed = strtoul(argv[++i], NULL, 10);
             g_master_seed_provided = 1;
         } else if (strcmp(argv[i], "--init-cond-seed") == 0) {
             if (i + 1 >= argc || !isInteger(argv[i + 1])) {
                 errorAndExit("--init-cond-seed requires an integer argument", argv[i + 1], argv[0]);
             }
             g_initial_cond_seed = strtoul(argv[++i], NULL, 10);
             g_initial_cond_seed_provided = 1;
         } else if (strcmp(argv[i], "--sidm-seed") == 0) {
             if (i + 1 >= argc || !isInteger(argv[i + 1])) {
                 errorAndExit("--sidm-seed requires an integer argument", argv[i + 1], argv[0]);
             }
             g_sidm_seed = strtoul(argv[++i], NULL, 10);
             g_sidm_seed_provided = 1;
         } else if (strcmp(argv[i], "--load-seeds") == 0) {
             g_attempt_load_seeds = 1;
         } else if (strcmp(argv[i], "--profile") == 0) {
             if (i + 1 >= argc) {
                 errorAndExit("--profile requires a type argument ('nfw' or 'cored')", NULL, argv[0]);
             }
             strncpy(g_profile_type_str, argv[++i], sizeof(g_profile_type_str) - 1);
             g_profile_type_str[sizeof(g_profile_type_str) - 1] = '\0'; // Ensure null termination
             if (strcmp(g_profile_type_str, "nfw") != 0 && strcmp(g_profile_type_str, "cored") != 0) {
                 errorAndExit("Invalid argument for --profile. Use 'nfw' or 'cored'.", g_profile_type_str, argv[0]);
             }
             g_profile_type_str_provided = 1;
         } else if (strcmp(argv[i], "--scale-radius") == 0) {
             if (i + 1 >= argc || !isFloat(argv[i + 1])) {
                 errorAndExit("--scale-radius requires a float argument", argv[i + 1], argv[0]);
             }
             g_scale_radius_param = atof(argv[++i]);
             if (g_scale_radius_param <= 0) errorAndExit("--scale-radius must be positive", NULL, argv[0]);
             g_scale_radius_param_provided = 1;
         } else if (strcmp(argv[i], "--halo-mass") == 0) {
             if (i + 1 >= argc || !isFloat(argv[i + 1])) { // Ensure isFloat is robust for scientific notation
                 errorAndExit("--halo-mass requires a float argument", argv[i + 1], argv[0]);
             }
             g_halo_mass_param = atof(argv[++i]);
             if (g_halo_mass_param <= 0) errorAndExit("--halo-mass must be positive", NULL, argv[0]);
             g_halo_mass_param_provided = 1;
         } else if (strcmp(argv[i], "--cutoff-factor") == 0) {
             if (i + 1 >= argc || !isFloat(argv[i + 1])) {
                 errorAndExit("--cutoff-factor requires a float argument", argv[i + 1], argv[0]);
             }
             g_cutoff_factor_param = atof(argv[++i]);
             if (g_cutoff_factor_param <= 0) errorAndExit("--cutoff-factor must be positive", NULL, argv[0]);
             g_cutoff_factor_param_provided = 1;
         } else if (strcmp(argv[i], "--falloff-factor") == 0) {
             if (i + 1 >= argc || !isFloat(argv[i + 1])) {
                 errorAndExit("--falloff-factor requires a float argument", argv[i + 1], argv[0]);
             }
             g_falloff_factor_param = atof(argv[++i]);
             if (g_falloff_factor_param <= 0) errorAndExit("--falloff-factor must be positive", NULL, argv[0]);
             g_falloff_factor_param_provided = 1;
         }
         else if (strncmp(argv[i], "--", 2) == 0)
         {
             errorAndExit("unrecognized option", argv[i], argv[0]);
         }
         else
         {
             errorAndExit("unrecognized argument", argv[i], argv[0]);
         }
     }
 }

 /**
  * @brief set method name.
  * @details updates the method_name variable based on the method_select value.
  *
  * @param None
  * @return None
  */
 void set_method_name() {
     switch (method_select)
     {
     case 1:
         method_name = "Adaptive Leapfrog with Adaptive Levi-Civita";
         strcpy(method_filename, "adp.leap.adp.levi");
         break;
     case 2:
         method_name = "Full-Step Adaptive Leapfrog + Levi-Civita";
         strcpy(method_filename, "adp.leap.levi");
         break;
     case 3:
         method_name = "Full-Step Adaptive Leapfrog";
         strcpy(method_filename, "adp.leap");
         break;
     case 4:
         method_name = "Yoshida 4th-Order";
         strcpy(method_filename, "fr4.yoshi");
         break;
     case 5:
         method_name = "Adams-Bashforth 3rd-Order";
         strcpy(method_name, "ab3");
         break;
     case 6:
         method_name = "Leapfrog (Vel Half-Step)";
         strcpy(method_filename, "vel.leap");
         break;
     case 7:
         method_name = "Leapfrog (Pos Half-Step)";
         strcpy(method_filename, "pos.leap");
         break;
     case 8:
         method_name = "Classic RK4";
         strcpy(method_filename, "rk4");
         break;
     case 9:
         method_name = "Euler";
         strcpy(method_filename, "euler");
         break;
     default:
         method_name = "Unknown Method";
         strcpy(method_filename, "unknown");
         break;
     }
 }

 /**
  * @brief Parse user arguments.
  * @details Reads out the user input arguments and update global fields.
  *
  * @param argc [in] The number of input arguments.
  * @param argv [in] The input arguments.
  * @return 1 if the "help" prompt was requested (leave main immedietly) or 0 otherwise.
  */
 int parse_user_arguments(int argc, char *argv[]) {
     /** @note Check for the `--help` argument first before parsing other options. */
     if (exit_on_help(argc, argv)) {
         return 1;
     }
     parse_user_arguments(argc, argv);
     set_method_name();
     return 0;
}
