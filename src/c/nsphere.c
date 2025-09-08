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
#include <stdarg.h>
#include <stdint.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <unistd.h>    /* For getpid */
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h> // For usleep function
#include "globals.h"
#include "particle_data.h"
#include "sidm.h"
#include "logging.h"
#include "particle_array_ops.h"
#include "exit.h"
#include "utils.h"
#include "signal_processing.h"
#include "gravitation_dynamics.h"
#include "density.h"
#include "density_nfw.h"
// #include "density_cored_plummer.h"
#include "debug.h"
#include "cli.h"
#include "io.h"
#include <float.h> // For DBL_MAX
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/statvfs.h>
#endif
#ifndef MYBINIO_H
#define MYBINIO_H
#endif

/**
 * @brief Adjusts the total number of timesteps to align with desired output snapshot intervals.
 * @details This function calculates an adjusted number of total simulation timesteps, \f$N'_{times}\f$,
 *          such that it is greater than or equal to the initially requested `Ntimes_initial` (\f$N\f$)
 *          and satisfies the constraint: \f$(N'_{times} - 1)\f$ must be an integer multiple of
 *          \f$(M - 1) \times p\f$. Here, \f$M\f$ is `nout` (number of desired output snapshot points,
 *          which means \f$M-1\f$ intervals) and \f$p\f$ is `dtwrite` (the low-level write interval
 *          in terms of simulation timesteps).
 *          This alignment ensures that exactly `nout` snapshots can be produced at intervals
 *          that are multiples of `dtwrite` and that also evenly span the total adjusted simulation duration.
 *
 * @param Ntimes_initial [in] Initially requested total number of simulation timesteps (\f$N\f$).
 * @param nout           [in] Number of desired output snapshot points (\f$M\f$). Must be >= 2 for adjustment to apply.
 * @param dtwrite        [in] The interval (in timesteps) at which low-level data is potentially written (\f$p\f$). Must be >= 1.
 * @return int The adjusted total number of timesteps (\f$N'_{times}\f$). Returns `Ntimes_initial`
 *             if `nout < 2` or `dtwrite < 1` or other edge cases where the constraint cannot be met.
 */
static int adjust_ntimesteps(int Ntimes_initial, int nout, int dtwrite)
{
    // Notation: M = nout, p = dtwrite, N = Ntimes_initial.
    // Find the smallest N' >= N such that (N' - 1) is a multiple of (M - 1) * p.

    int M = nout;
    int p = dtwrite;
    int N = Ntimes_initial;

    // Edge cases:
    if (M < 2)
    {
        // If only 0 or 1 snapshot requested, no interval constraint applies.
        return N;
    }
    if (p < 1)
    {
        // Invalid write interval.
        return N;
    }

    // The total number of intervals between M snapshots is (M - 1).
    // The total number of steps spanning these intervals must be a multiple of p.
    // Therefore, the total number of steps (N' - 1) must be a multiple of (M - 1) * p.
    // Find the smallest integer k >= 1 such that (M - 1) * k * p >= (N - 1).
    double required_steps = (double)(N - 1);
    double steps_per_output_cycle = (M - 1) * (double)p;

    // Handle case where denominator is zero (e.g., M=1 or p=0, caught above but added safety)
    if (steps_per_output_cycle <= 0) {
        return N; // Cannot satisfy constraint
    }

    double ratio = required_steps / steps_per_output_cycle;
    int k = (int)ceil(ratio);
    if (k < 1)
    {
        k = 1; // Ensure at least one full output cycle.
    }

    int Nprime_minus_1 = (M - 1) * k * p;
    int Nprime = Nprime_minus_1 + 1;

    return Nprime;
}

/**
 * @brief Main entry point for the n-sphere dark matter simulation program.
 * @details Orchestrates the overall simulation workflow:
 *          1. Parses command-line arguments.
 *          2. Sets up global parameters and logging.
 *          3. Initializes random number generators.
 *          4. Generates or loads initial conditions (ICs) for either NFW or Cored Plummer-like profiles.
 *             - Includes theoretical calculations for density, mass, potential, and f(E) splines.
 *             - Includes a diagnostic loop to test IC generation with varied numerical parameters.
 *             - Performs particle sampling based on the derived distribution function.
 *          5. Optionally performs tidal stripping and re-assigns particle IDs.
 *          6. Converts particle velocities to physical simulation units.
 *          7. Executes the main N-body timestepping loop using a selected integration method.
 *             - Performs gravitational updates.
 *             - Optionally performs SIDM scattering via `handle_sidm_step`.
 *             - Tracks particle trajectories and energies.
 *             - Periodically writes simulation data and progress.
 *          8. If `g_doAllParticleData` is enabled, processes all particle data to generate
 *             snapshot files for Rank/Mass/Radius/Velocity/Potential/Energy/Density.
 *          9. Writes final summary plots and theoretical profiles.
 *          10. Cleans up allocated resources.
 *          Handles restart/resume functionality by checking for existing data products.
 *
 * @param argc [in] Standard argument count from the command line.
 * @param argv [in] Standard array of argument strings from the command line.
 * @return int Exit code: 0 for successful execution, non-zero for errors.
 *
 * @note This application supports OpenMP for parallelization in various sections.
 * @warning Large particle counts or long simulations can be memory and CPU intensive.
 *          Disk space requirements for full data output can also be significant.
 */
int main(int argc, char *argv[])
{
// Print the tool header first, regardless of OpenMP status
printf("\n===================================================================================================\n");
printf("NSphere Simulation Tool\n");
printf("===================================================================================================\n");
printf("  \n");

#ifdef _OPENMP
    /** @note OpenMP section: Configures parallel execution environment when compiled with OpenMP. */
    // OpenMP is available - configure parallel execution environment
    int max_threads = omp_get_max_threads();
    int num_processors = omp_get_num_procs();

    // Enable nested parallelism with max_active_levels (replacing deprecated omp_set_nested)
    omp_set_max_active_levels(10); // Allow up to 10 levels of nested parallelism

    // Set OpenMP to use maximum available thread parallelism
    omp_set_num_threads(max_threads);

    printf("OpenMP Status: ENABLED (%d logical processors, using %d threads)\n",
           num_processors, max_threads);
    printf("\n\n");

    log_message("INFO", "Using maximum thread parallelism: %d threads", max_threads);
#else
    /** @warning OpenMP section: Warns user when compiled without OpenMP support. */
    // OpenMP is not available - warn user about performance implications
    printf("WARNING: OpenMP is NOT ENABLED in this build!\n");
    printf("This will result in significantly reduced performance.\n");
    printf("For better performance, please install OpenMP and recompile with -fopenmp flag.\n\n\n");

    log_message("WARNING", "OpenMP not available - running in single-threaded mode");

    // Check if --help flag is used (don't delay in that case)
    int help_requested = 0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            help_requested = 1;
            break;
        }
    }

    // Add a delay to ensure the warning is noticed, but only when not showing help
    if (!help_requested) {
        printf("Continuing in single-threaded mode");
        fflush(stdout);
        for (int i = 0; i < 5; i++) {
            usleep(500000); // 500ms * 5 = 2.5 seconds
            printf(".");
            fflush(stdout);
        }
        printf("\n\n");
    }

    // Define a single thread variable for the code below to use
    int max_threads __attribute__((unused)) = 1;
#endif

#ifdef _OPENMP
    /** @note Initializes FFTW thread support if compiled with OpenMP. */
    // Initialize FFTW thread support (only effective if FFTW was compiled with threading)
    fftw_init_threads();
    fftw_plan_with_nthreads(max_threads);
#endif

    int npts = 100000;
    int Ntimes = 10000;
    int tfinal_factor = 5;
    int nout = 100;
    int dtwrite = 100;
    double tidal_fraction = 0.0;
    int noutsnaps;
    g_total_sidm_scatters = 0; // Initialize global SIDM scatter counter

    int method_select = 1;            // Default: option 1 (Adaptive Leapfrog with Adaptive Levi-Civita)
    int display_sort = 1;             // Default: option 1 (Parallel Quadsort)
    int include_method_in_suffix = 0; // Default: exclude method from filenames.
    char custom_tag[256] = {0};       // Default: no custom tag.

    /** @note Check for the `--help` argument first before parsing other options. */
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--help") == 0)
        {
            printUsage(argv[0]);
            return 0;
        }
    }

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

    /** @note Convert user-facing method number (1-9) to internal identifier (0-8) and get description string. */
    int display_method = method_select; // Store original user input for display
    char *method_verbose_name;

    // Convert method number to internal index and set descriptive name.
    switch (method_select)
    {
    case 1:
        method_select = 5;
        method_verbose_name = "Adaptive Leapfrog with Adaptive Levi-Civita";
        break;
    case 2:
        method_select = 4;
        method_verbose_name = "Full-Step Adaptive Leapfrog + Levi-Civita";
        break;
    case 3:
        method_select = 3;
        method_verbose_name = "Full-Step Adaptive Leapfrog";
        break;
    case 4:
        method_select = 6;
        method_verbose_name = "Yoshida 4th-Order";
        break;
    case 5:
        method_select = 8;
        method_verbose_name = "Adams-Bashforth 3rd-Order";
        break;
    case 6:
        method_select = 2;
        method_verbose_name = "Leapfrog (Vel Half-Step)";
        break;
    case 7:
        method_select = 1;
        method_verbose_name = "Leapfrog (Pos Half-Step)";
        break;
    case 8:
        method_select = 7;
        method_verbose_name = "Classic RK4";
        break;
    case 9:
        method_select = 0;
        method_verbose_name = "Euler";
        break;
    default:
        method_verbose_name = "Unknown Method";
        break;
    }

    /** @note Generate internal method string identifier for filenames. */
    char method_str[32];
    switch (method_select)
    {
    case 0:
        strcpy(method_str, "euler");
        break;
    case 1:
        strcpy(method_str, "pos.leap");
        break;
    case 2:
        strcpy(method_str, "vel.leap");
        break;
    case 3:
        strcpy(method_str, "adp.leap");
        break;
    case 4:
        strcpy(method_str, "adp.leap.levi");
        break;
    case 5:
        strcpy(method_str, "adp.leap.adp.levi");
        break;
    case 6:
        strcpy(method_str, "fr4.yoshi");
        break;
    case 7:
        strcpy(method_str, "rk4");
        break;
    case 8:
        strcpy(method_str, "ab3");
        break;
    default:
        strcpy(method_str, "unknown");
        break;
    }

    // Determine active profile type (NFW is default)
    if (g_profile_type_str_provided) {
        if (strcmp(g_profile_type_str, "nfw") == 0) {
            g_use_nfw_profile = 1;
        } else if (strcmp(g_profile_type_str, "cored") == 0) {
            g_use_nfw_profile = 0;
        } else {
            // Should have been caught by parser, but as a safeguard:
            log_message("WARNING", "Unknown profile type '%s', defaulting to NFW.", g_profile_type_str);
            g_use_nfw_profile = 1;
        }
    } else {
        // Default to NFW if --profile flag was not provided
        g_use_nfw_profile = 1;
        strcpy(g_profile_type_str, "nfw"); // Update string for consistency in printouts
    }

    // Set up profile-specific parameters based on generalized flags and profile defaults
    if (g_use_nfw_profile) {
        // NFW Profile Path
        // Halo Mass for NFW
        if (g_halo_mass_param_provided) { // --halo-mass overrides NFW default
            g_nfw_profile_halo_mass = g_halo_mass_param;
        } else { // No --halo-mass, NFW uses its own default
            g_nfw_profile_halo_mass = HALO_MASS_NFW;
            g_halo_mass_param = g_nfw_profile_halo_mass; // Update general param to reflect NFW's choice
        }
        // Scale Radius for NFW
        if (g_scale_radius_param_provided) { // --scale-radius overrides NFW default
            g_nfw_profile_rc = g_scale_radius_param;
        } else { // No --scale-radius, NFW uses its own default
            g_nfw_profile_rc = RC_NFW_DEFAULT;
            g_scale_radius_param = g_nfw_profile_rc; // Update general param to reflect NFW's choice
        }
        // Cutoff Factor for NFW
        if (g_cutoff_factor_param_provided) { // --cutoff-factor overrides NFW default
            g_nfw_profile_rmax_norm_factor = g_cutoff_factor_param;
        } else { // No --cutoff-factor, NFW uses its own default
            g_nfw_profile_rmax_norm_factor = CUTOFF_FACTOR_NFW_DEFAULT;
            // g_cutoff_factor_param is NOT updated here by NFW default; it keeps its own (Cored's) default or user value.
        }
        // Falloff Factor for NFW
        if (g_falloff_factor_param_provided) { // --falloff-factor overrides NFW default
            g_nfw_profile_falloff_factor = g_falloff_factor_param;
        } else { // No --falloff-factor, NFW uses its own default
            g_nfw_profile_falloff_factor = FALLOFF_FACTOR_NFW_DEFAULT;
            // Optionally, update g_falloff_factor_param if NFW is the overall default and no flag given
            // For now, let g_falloff_factor_param keep its own default unless explicitly set by user
        }

    } else {
        // Cored Profile Path
        // Halo Mass for Cored (already defaults to HALO_MASS or takes from --halo-mass via g_halo_mass_param)
        g_cored_profile_halo_mass = g_halo_mass_param;
        // Scale Radius for Cored (already defaults to RC or takes from --scale-radius via g_scale_radius_param)
        g_cored_profile_rc = g_scale_radius_param;
        // Cutoff Factor for Cored (directly uses generalized or its (Cored's) default)
        g_cored_profile_rmax_factor = g_cutoff_factor_param;
    }

    // Set the single g_active_halo_mass for N-body forces and tdyn from the finalized g_halo_mass_param
    g_active_halo_mass = g_halo_mass_param;

    /** @note Display final parameter values used for the simulation run. */
    printf("Parameter values requested:\n\n");

    printf("  Number of Particles:          %d\n", npts);
    printf("  Number of Time Steps:         %d\n", Ntimes);
    printf("  Number of Dynamical Times:    %d\n", tfinal_factor);
    printf("  Number of Output Snapshots:   %d\n", nout);
    printf("  Steps Between Writes:         %d\n", dtwrite);
    printf("  Tidal Stripping Fraction:     %.5f\n", tidal_fraction);
    printf("  Integration Method:           %d (%s)\n", display_method, method_verbose_name);
    printf("  Sorting Algorithm:            %d (%s)\n", display_sort, get_sort_description(g_defaultSortAlg));
    /** @note Build the filename tag string based on options for display purposes. */
    char filename_tag[512] = "";

    if (custom_tag[0] != '\0')
    {
        strcat(filename_tag, custom_tag);
    }

    if (include_method_in_suffix)
    {
        if (filename_tag[0] != '\0')
        {
            strcat(filename_tag, "_");
        }
        strcat(filename_tag, method_str);
    }

    printf("  Filename Tag:                 %s\n", filename_tag[0] ? filename_tag : "[none]");
    printf("  SIDM Scattering:              %s\n", g_enable_sidm_scattering ? "Enabled via --sidm" : "Disabled (Default)");
    printf("  SIDM Execution Mode:          %s\n", g_sidm_execution_mode == 1 ? "Parallel (Default)" : "Serial");
    printf("  SIDM Opacity Kappa:           %.1f cm^2/g (Default: 50.0, User set: %s)\n", g_sidm_kappa, g_sidm_kappa_provided ? "Yes" : "No");

    printf("  Initial Conditions Profile:   %s\n", g_use_nfw_profile ? "NFW-like with Cutoff" : "Cored Plummer-like");
    if (g_use_nfw_profile) {
        printf("    NFW Profile Scale Radius (IC): %.3f kpc (NFW Default: %.2f, User set via --scale-radius: %s)\n", g_nfw_profile_rc, RC_NFW_DEFAULT, g_scale_radius_param_provided ? "Yes" : "No");
    } else {
        printf("    Cored Profile Scale Radius (IC): %.3f kpc (Cored Default: %.2f, User set via --scale-radius: %s)\n", g_cored_profile_rc, RC, g_scale_radius_param_provided ? "Yes" : "No");
    }
    if (g_use_nfw_profile) {
        printf("    NFW Profile Halo Mass (IC): %.3e Msun (NFW Default: %.2e, User set via --halo-mass: %s)\n", g_nfw_profile_halo_mass, HALO_MASS_NFW, g_halo_mass_param_provided ? "Yes" : "No");
    } else {
        printf("    Cored Profile Halo Mass (IC): %.3e Msun (Cored Default: %.2e, User set via --halo-mass: %s)\n", g_cored_profile_halo_mass, HALO_MASS, g_halo_mass_param_provided ? "Yes" : "No");
    }
    printf("    Profile Cutoff Factor:      %.1f (CmdLine/Default: %.1f, User set: %s)\n", g_cutoff_factor_param, (g_use_nfw_profile ? CUTOFF_FACTOR_NFW_DEFAULT : CUTOFF_FACTOR_CORED_DEFAULT), g_cutoff_factor_param_provided ? "Yes" : "No");

    if (g_use_nfw_profile) {
        printf("    NFW Profile Falloff Factor (C): %.1f (NFW Default: %.1f, User set via --falloff-factor: %s)\n", g_nfw_profile_falloff_factor, FALLOFF_FACTOR_NFW_DEFAULT, g_falloff_factor_param_provided ? "Yes" : "No");
    }
    // This g_active_halo_mass is now correctly set from g_halo_mass_param which reflects the chosen profile's mass
    printf("  N-body Active Halo Mass (tdyn): %.3e Msun\n", g_active_halo_mass);

    /** @note Display logging status based on g_enable_logging flag. */
    if (g_enable_logging)
    {
        printf("  Logging:                      Enabled (log/nsphere.log)\n\n");
        log_message("INFO", "Simulation started with %d particles, %d timesteps, %d dynamical times",
                    npts, Ntimes, tfinal_factor);
    }
    else
    {
        printf("  Logging:                      Disabled\n\n");
    }

    // Check for SIDM + parallel mode without OpenMP
    if (g_enable_sidm_scattering && g_sidm_execution_mode == 1) {
        #ifndef _OPENMP
            printf("Warning: SIDM parallel mode is enabled by default, but OpenMP is not available in this build.\n");
            printf("         SIDM will run serially. Use '--sidm-mode serial' to suppress this warning.\n\n");
            log_message("WARNING", "SIDM parallel mode requested but OpenMP not available, will run serially.");
            g_sidm_execution_mode = 0; // Force serial if no OpenMP
        #endif
    }

    // Validation occurs earlier in the argument parsing loop.

    // ... proceed with simulation ...

    /**
     * @brief Oversample initial conditions based on tidal fraction.
     * @details Calculate the number of initial particles (`npts_initial`) needed
     *          before tidal stripping to ensure `npts` particles remain afterwards.
     *          If `tidal_fraction` is 0, `npts_initial` equals `npts`.
     */
    int npts_initial;
    if (tidal_fraction > 0.0)
    {
        // Use ceiling to ensure enough particles remain after stripping
        npts_initial = ceil(npts / (1.0 - tidal_fraction));
    }
    else
    {
        npts_initial = npts;
    }

    /** @note Calculate number of snapshots (`noutsnaps`) based on desired intervals (`nout`). */
    nout = nout + 1; // nout specifies intervals, noutsnaps is number of points (intervals + 1)
    noutsnaps = nout;

    /**
     * @brief Adjust Ntimes using adjust_ntimesteps to align with output schedule.
     * @details Ensures (Ntimes - 1) is a multiple of (noutsnaps - 1) * dtwrite.
     * @see adjust_ntimesteps
     */
    int oldN = Ntimes;
    Ntimes = adjust_ntimesteps(Ntimes, noutsnaps, dtwrite); // Use noutsnaps here
    if (Ntimes != oldN)
    {
        printf("Adjusted Number of Time Steps to %d to satisfy parameter constraints.\n", Ntimes);
    }
    /** @brief Calculate total number of write events and steps between major snapshots. */
    int total_writes = ((Ntimes - 1) / dtwrite) + 1; // Total potential write points
    int stepBetweenSnaps = (int)floor(
        (double)(total_writes - 1) / (double)(noutsnaps - 1) + 0.5); // Steps between major snapshots
    int ext_Ntimes; ///< Extended time steps potentially needed for trajectory arrays bounds.
    ext_Ntimes = Ntimes + dtwrite; // Allocate trajectory arrays slightly larger

    /** @brief Initialize the global file suffix string `g_file_suffix` based on cmd line args. */
    g_file_suffix[0] = '\0';

    /** @brief Add custom tag to suffix if provided via `--tag`. */
    if (custom_tag[0] != '\0')
    {
        snprintf(g_file_suffix, sizeof(g_file_suffix), "_%s", custom_tag);
    }

    /** @brief Add method/parameter tag to suffix. */
    char temp[256];
    if (include_method_in_suffix)
    {
        snprintf(temp, sizeof(temp), "_%s_%d_%d_%d", method_str, npts, Ntimes, tfinal_factor);
    }
    else
    {
        snprintf(temp, sizeof(temp), "_%d_%d_%d", npts, Ntimes, tfinal_factor);
    }

    /** @brief Append the parameter tag to the global suffix. */
    strcat(g_file_suffix, temp); // Append temp to g_file_suffix

    /** @brief Set flag `skip_file_writes=1` if in restart mode (`g_doRestart`). */
    if (g_doRestart)
    {
        printf("Restart mode: Skipping simulation phase and file writes, proceeding directly to data product generation.\n\n");
        skip_file_writes = 1;
    }

    /** @brief Create the output 'data' directory if it doesn't exist. */
    {
        struct stat st = {0};
        if (stat("data", &st) == -1)
        {
            mkdir("data", 0755); // POSIX standard, works on most systems including MinGW/Cygwin
        }
    }

    /** @brief Create the 'init' directory if it doesn't exist. */
    {
        struct stat st_init = {0};
        if (stat("init", &st_init) == -1)
        {
            #if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
                if (mkdir("init") != 0) {
                     perror("Error creating init directory");
                     // Decide if this is fatal - perhaps not if only writing
                } else {
                     log_message("INFO", "Created init/ directory.");
                }
            #else
                if (mkdir("init", 0755) != 0) { // POSIX standard
                     perror("Error creating init directory");
                     // Decide if this is fatal
                } else {
                     log_message("INFO", "Created init/ directory.");
                }
            #endif
        }
    }

    /** @brief Write current run parameters to `data/lastparams<suffix>.dat` and create a standard link `data/lastparams.dat`. */
    {
        char file_tag[512] = "";
        if (custom_tag[0] != '\0')
        {
            strcat(file_tag, custom_tag);
            if (include_method_in_suffix)
            {
                strcat(file_tag, "_");
                strcat(file_tag, method_str);
            }
        }
        else if (include_method_in_suffix)
        {
            strcat(file_tag, method_str);
        }

        /** @note Create filename with suffix for the specific run parameters file. */
        char filename[512]; // Holds suffixed filename, e.g., data/lastparams_run1_100k_10k_5.dat
        get_suffixed_filename("data/lastparams.dat", 1, filename, sizeof(filename));
        printf("Saving parameters: %s\n", filename);

        FILE *fp_params = fopen(filename, "w"); // Text mode for regular fprintf
        if (!fp_params)
        {
            printf("Error: cannot open %s\n", filename);
            return 1;
        }

        fprintf(fp_params, "%d %d %d %s\n", npts, Ntimes, tfinal_factor, file_tag);
        fclose(fp_params);

        /** @note Create a standard-named link `data/lastparams.dat` pointing to the
                 suffixed version for compatibility with scripts. Uses copy on Windows. */
        char linkname[512] = "data/lastparams.dat"; // Standard name

        /* Platform detection using standard predefined macros */
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
        /** @note Windows: Copy file content as symlinks can be unreliable or require special privileges. */
        // Windows or Windows-like environment: create a direct file copy.
            // Ensure compatibility with various Windows environments.

        // Use lower-level file operations instead of system commands for better compatibility.
        FILE *source, *dest;
        source = fopen(filename, "rb");
        if (!source)
        {
            printf("Warning: Failed to open source file %s for copying\n", filename);
        }
        else
        {
            dest = fopen(linkname, "wb");
            if (!dest)
            {
                printf("Warning: Failed to create destination file %s\n", linkname);
                fclose(source);
            }
            else
            {
                // Copy file content.
                char buffer[4096];
                size_t bytes_read;

                while ((bytes_read = fread(buffer, 1, sizeof(buffer), source)) > 0)
                {
                    fwrite(buffer, 1, bytes_read, dest);
                }

                fclose(dest);
                fclose(source);

                // Extract the basename for display purposes
                const char *basename = strrchr(filename, '/');
                basename = basename ? basename + 1 : filename; // Skip the '/' or use full name if no '/'

                printf("Created link: %s -> %s\n\n", basename, linkname);
            }
        }
#else
        /** @note Unix: Create symbolic link from basename(filename) to 'linkname', fallback to copy. */
            // Unix-like systems (Linux, macOS, etc.) and fallback for other platforms: use symbolic links.
        char command[1024];

        // First, remove any existing link or file.
        snprintf(command, sizeof(command), "rm -f \"%s\" 2>/dev/null", linkname);
        system(command);

        // Then create the symbolic link - use the basename of the file, not the full path
        // Extract the basename from filename
        const char *basename = strrchr(filename, '/');
        basename = basename ? basename + 1 : filename; // Skip the '/' or use full name if no '/'

        snprintf(command, sizeof(command), "ln -s \"%s\" \"%s\"", basename, linkname);
        if (system(command) != 0)
        {
            // If symbolic link fails, fall back to copying the file.
            snprintf(command, sizeof(command), "cp \"%s\" \"%s\"", filename, linkname);
            if (system(command) != 0)
            {
                printf("Warning: Failed to create link or copy %s to %s\n", filename, linkname);
            }
            else
            {
                printf("Created link: %s -> %s\n\n", filename, linkname);
            }
        }
        else
        {
            printf("Created link: %s -> %s\n\n", filename, linkname);
        }
#endif
    }

    /**
     * @brief Common IC generation variables shared between profile pathways.
     * @details These variables are declared before the profile selection block
     *          and will be populated by whichever profile pathway is chosen.
     */
    double **particles = NULL;           ///< Main particle data array
    int i = 0;                          ///< Loop counter
    double result, error;               ///< GSL integration results
    double calE;                        ///< Energy value for calculations
    gsl_integration_workspace *w = NULL; ///< GSL integration workspace

    // Common spline objects and accelerators
    gsl_spline *splinemass = NULL;      ///< Spline for mass profile M(r)
    gsl_interp_accel *enclosedmass = NULL; ///< Accelerator for mass spline
    gsl_spline *splinePsi = NULL;       ///< Spline for potential profile Psi(r)
    gsl_interp_accel *Psiinterp = NULL; ///< Accelerator for potential spline
    gsl_spline *splinerofPsi = NULL;    ///< Spline for inverse potential r(Psi)
    gsl_interp_accel *rofPsiinterp = NULL; ///< Accelerator for r(Psi) spline
    gsl_interp *g_main_fofEinterp = NULL;  ///< Main f(E) interpolator
    gsl_interp_accel *g_main_fofEacc = NULL; ///< Accelerator for f(E)

    // Common data arrays
    double *radius = NULL;              ///< Radial grid points
    double *mass = NULL;                ///< Mass values at radial points
    double *Psivalues = NULL;           ///< Potential values at radial points
    double *nPsivalues = NULL;          ///< Negative potential values (for r(Psi) spline)
    double *Evalues = NULL;             ///< Energy grid points
    double *innerintegrandvalues = NULL; ///< f(E) integrand values
    double *radius_monotonic_grid_nfw = NULL; ///< Monotonic radial grid for NFW calculations

    // Key scalar values
    double Psimin = 0.0;                ///< Minimum potential (at rmax)
    double Psimax = 0.0;                ///< Maximum potential (at r=0)
    double rmax = 0.0;                  ///< Maximum radius for profile calculations
    int num_points = 0;                 ///< Number of points for spline interpolation

    // File handling
    char fname[256];                    ///< Buffer for file names
    FILE *fp;                           ///< File pointer for data output

    /**
     * @brief Allocate main particle data array before profile selection.
     * @details This ensures both NFW and Cored pathways use the same particles array.
     */
    particles = (double **)malloc(5 * sizeof(double *));
    if (particles == NULL) {
        fprintf(stderr, "ERROR: Memory allocation failed for particle array pointer\n");
        CLEAN_EXIT(1);
    }
    for (i = 0; i < 5; i++) {
        particles[i] = (double *)malloc(npts_initial * sizeof(double));
        if (particles[i] == NULL) {
            fprintf(stderr, "ERROR: Memory allocation failed for particles[%d]\n", i);
            CLEAN_EXIT(1);
        }
    }

    /**
     * @brief Initialize random number generator for particle generation.
     * @details Sets up the GSL Random Number Generator environment and allocates
     *          the global GSL RNG state used throughout the simulation.
     *          Needed for the Sample Generator if not reading ICs from file.
     */
    gsl_rng_env_setup();                           // Setup GSL RNG environment
    const gsl_rng_type * T_rng = gsl_rng_default;  // Use default RNG type
    g_rng = gsl_rng_alloc(T_rng);                  // Allocate global GSL RNG state
    if (g_rng == NULL) {                           // Check allocation success
        fprintf(stderr, "Error allocating GSL RNG.\n");
        CLEAN_EXIT(1);
    }

    // Seed the global g_rng (used for ICs and Serial SIDM)
    gsl_rng_set(g_rng, g_initial_cond_seed);       // Use the determined IC seed for g_rng
    log_message("INFO", "Global g_rng (intended primarily for IC generation) seeded with %lu", g_initial_cond_seed);

    // Initialize per-thread GSL RNGs if OpenMP is enabled
    #ifdef _OPENMP
        g_max_omp_threads_for_rng = omp_get_max_threads();
        if (g_max_omp_threads_for_rng <= 0) g_max_omp_threads_for_rng = 1; // Safety
    #else
        g_max_omp_threads_for_rng = 1;
    #endif

    g_rng_per_thread = (gsl_rng **)malloc(g_max_omp_threads_for_rng * sizeof(gsl_rng *));
    if (g_rng_per_thread == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for per-thread RNG array.\n");
        CLEAN_EXIT(1);
    }

    const gsl_rng_type *T_rng_thread = gsl_rng_default;

    for (int i_rng = 0; i_rng < g_max_omp_threads_for_rng; ++i_rng) {
        g_rng_per_thread[i_rng] = gsl_rng_alloc(T_rng_thread);
        if (g_rng_per_thread[i_rng] == NULL) {
            fprintf(stderr, "Error: Failed to allocate GSL RNG for thread %d.\n", i_rng);
            for (int k_rng = 0; k_rng < i_rng; ++k_rng) gsl_rng_free(g_rng_per_thread[k_rng]);
            free(g_rng_per_thread);
            CLEAN_EXIT(1);
        }
        gsl_rng_set(g_rng_per_thread[i_rng], g_sidm_seed + (unsigned long int)i_rng);
    }
    log_message("INFO", "Initialized %d per-thread GSL RNGs (for SIDM) using base SIDM seed %lu", g_max_omp_threads_for_rng, g_sidm_seed);

    if (g_use_nfw_profile) {
        log_message("INFO", "Starting IC generation using NFW-like profile pathway.");
        log_message("INFO", "Generating Initial Conditions using NFW-like profile with its specific numerics...");

        // Diagnostic loop for NFW (similar to Cored's diagnostic)
        if (g_doDebug) {
            log_message("INFO", "NFW DIAGNOSTIC LOOP: Starting convergence tests for NFW profile.");
            int diag_integration_points_array[2] = {1000, 10000};
            int diag_spline_points_array[2] = {1000, 10000};

            for (int ii_ip_nfw = 0; ii_ip_nfw < 2; ii_ip_nfw++) {
                for (int ii_sp_nfw = 0; ii_sp_nfw < 2; ii_sp_nfw++) {
                    int Nintegration_diag = diag_integration_points_array[ii_ip_nfw];
                    int Nspline_diag_base = diag_spline_points_array[ii_sp_nfw];
                    int num_points_diag = Nspline_diag_base * 10;

                    log_message("DEBUG", "Diagnostic iteration %d, spline_base=%d (points=%d)",
                                Nintegration_diag, Nspline_diag_base, num_points_diag);

                    // Use the main NFW profile parameters for this diagnostic run
                    double current_diag_rc = g_nfw_profile_rc;
                    double current_diag_halo_mass = g_nfw_profile_halo_mass;
                    double current_diag_rmax_factor = g_nfw_profile_rmax_norm_factor;
                    double current_diag_falloff_C = g_nfw_profile_falloff_factor;
                    double rmax_diag = current_diag_rmax_factor * current_diag_rc;

                    // Local GSL workspace and variables for this diagnostic iteration
                    gsl_integration_workspace *w_diag = gsl_integration_workspace_alloc(Nintegration_diag);
                    if (!w_diag) {
                        log_message("ERROR", "Failed to allocate GSL workspace for NFW diagnostic");
                        continue;
                    }

                    double result_diag, error_diag;
                    double normalization_diag;

                    // Declare all splines and accelerators locally for the diagnostic loop
                    gsl_spline *splinemass_diag = NULL;
                    gsl_interp_accel *enclosedmass_diag = NULL;
                    gsl_spline *splinePsi_diag = NULL;
                    gsl_interp_accel *Psiinterp_diag = NULL;
                    gsl_spline *splinerofPsi_diag = NULL;
                    gsl_interp_accel *rofPsiinterp_diag = NULL;
                    gsl_interp *fofEinterp_diag = NULL;
                    gsl_interp_accel *fofEacc_diag = NULL;

                    double *mass_diag_arr = NULL;
                    double *radius_diag_arr = NULL;
                    double *radius_for_rofPsi_diag_arr = NULL;
                    double *Psivalues_diag_arr = NULL;
                    double *nPsivalues_diag_arr = NULL;
                    double *Evalues_diag_arr = NULL;
                    double *innerintegrandvalues_diag_arr = NULL;

                    // Prepare NFW params for this diagnostic iteration's integrands
                    double nfw_params_diag[4];
                    nfw_params_diag[0] = current_diag_rc;
                    nfw_params_diag[1] = current_diag_halo_mass;
                    nfw_params_diag[2] = 1.0; // Initial nt_nfw guess
                    nfw_params_diag[3] = current_diag_falloff_C;

                    gsl_function F_nfw_diag;
                    F_nfw_diag.function = &massintegrand_profile_nfwcutoff;
                    F_nfw_diag.params = nfw_params_diag;

                    // --- 1. Normalization for NFW Diagnostic ---
                    gsl_integration_qag(&F_nfw_diag, 0.0, rmax_diag, 1e-12, 1e-12, Nintegration_diag, GSL_INTEG_GAUSS51,
                                        w_diag, &result_diag, &error_diag);
                    normalization_diag = result_diag;
                    if (normalization_diag <= 1e-30) {
                        log_message("ERROR", "NFW diagnostic normalization too small: %e", normalization_diag);
                        gsl_integration_workspace_free(w_diag);
                        continue;
                    }
                    nfw_params_diag[2] = current_diag_halo_mass / (4.0 * M_PI * normalization_diag);

                    // --- 2. M(r) spline for NFW Diagnostic ---
                    mass_diag_arr = (double *)malloc(num_points_diag * sizeof(double));
                    radius_diag_arr = (double *)malloc(num_points_diag * sizeof(double));
                    if (!mass_diag_arr || !radius_diag_arr) {
                        log_message("ERROR", "Failed to allocate arrays for NFW diagnostic M(r)");
                        gsl_integration_workspace_free(w_diag);
                        free(mass_diag_arr);
                        free(radius_diag_arr);
                        continue;
                    }

                    mass_diag_arr[0] = 0.0;
                    radius_diag_arr[0] = 0.0;
                    for (int k = 1; k < num_points_diag; k++) {
                        double r_k = (double)k * rmax_diag / (num_points_diag - 1.0);
                        if (k == num_points_diag - 1) r_k = rmax_diag;
                        radius_diag_arr[k] = r_k;
                        gsl_integration_qag(&F_nfw_diag, 0.0, r_k, 1e-12, 1e-12, Nintegration_diag, GSL_INTEG_GAUSS51,
                                            w_diag, &result_diag, &error_diag);
                        mass_diag_arr[k] = 4.0 * M_PI * result_diag;
                    }
                    enclosedmass_diag = gsl_interp_accel_alloc();
                    splinemass_diag = gsl_spline_alloc(gsl_interp_cspline, num_points_diag);
                    gsl_spline_init(splinemass_diag, radius_diag_arr, mass_diag_arr, num_points_diag);


                    // Write mass profile diagnostic file
                    char diag_fname_mass[256];
                    char diag_base_mass[128];
                    snprintf(diag_base_mass, sizeof(diag_base_mass), "data/massprofile_Ni%d_Ns%d.dat", Nintegration_diag, Nspline_diag_base);
                    get_suffixed_filename(diag_base_mass, 1, diag_fname_mass, sizeof(diag_fname_mass));
                    FILE *fp_diag_mass = fopen(diag_fname_mass, "wb");
                    if (fp_diag_mass) {
                        for (double r_write_diag = 0.0; r_write_diag < radius_diag_arr[num_points_diag - 1]; r_write_diag += rmax_diag / 900.0) {
                            if (r_write_diag >= radius_diag_arr[0]) {
                                 fprintf_bin(fp_diag_mass, "%f %f\n", r_write_diag, gsl_spline_eval(splinemass_diag, r_write_diag, enclosedmass_diag));
                            }
                        }
                        fclose(fp_diag_mass);
                        log_message("DEBUG", "Wrote diagnostic file: %s", diag_fname_mass);
                    } else {
                        log_message("ERROR", "Failed to open diagnostic file: %s", diag_fname_mass);
                    }

                    // --- 3. Psi(r) spline for NFW Diagnostic ---
                    Psivalues_diag_arr = (double *)malloc(num_points_diag * sizeof(double));
                    nPsivalues_diag_arr = (double *)malloc(num_points_diag * sizeof(double));
                    radius_for_rofPsi_diag_arr = (double *)malloc(num_points_diag * sizeof(double));
                    if (!Psivalues_diag_arr || !nPsivalues_diag_arr || !radius_for_rofPsi_diag_arr) {
                        log_message("ERROR", "Failed to allocate arrays for NFW diagnostic Psi(r)");
                        goto cleanup_diag_iteration;
                    }

                    Psiintegrand_params psi_params_nfw_diag;
                    psi_params_nfw_diag.massintegrand_func = &massintegrand_profile_nfwcutoff;
                    psi_params_nfw_diag.params_for_massintegrand = nfw_params_diag;
                    gsl_function F_psi_nfw_diag;
                    F_psi_nfw_diag.function = &Psiintegrand;
                    F_psi_nfw_diag.params = &psi_params_nfw_diag;

                    for (int k = 0; k < num_points_diag; k++) {
                        double r_k = radius_diag_arr[k];
                        double r1_psi_k = fmax(r_k, current_diag_rc / 1000000.0);
                        gsl_integration_qagiu(&F_psi_nfw_diag, r1_psi_k, 1e-12, 1e-12, Nintegration_diag,
                                              w_diag, &result_diag, &error_diag);
                        double M_at_r1_k = gsl_spline_eval(splinemass_diag, r1_psi_k, enclosedmass_diag);
                        double first_term_psi_k = (r1_psi_k > 1e-9) ? (G_CONST * M_at_r1_k / r1_psi_k) : 0.0;
                        double second_term_psi_k = G_CONST * 4.0 * M_PI * result_diag;
                        Psivalues_diag_arr[k] = first_term_psi_k + second_term_psi_k;
                        nPsivalues_diag_arr[k] = -Psivalues_diag_arr[k];
                        radius_for_rofPsi_diag_arr[k] = r_k;
                    }
                    Psiinterp_diag = gsl_interp_accel_alloc();
                    splinePsi_diag = gsl_spline_alloc(gsl_interp_cspline, num_points_diag);
                    gsl_spline_init(splinePsi_diag, radius_diag_arr, Psivalues_diag_arr, num_points_diag);


                    // Write potential profile diagnostic file
                    char diag_fname_psi[256];
                    char diag_base_psi[128];
                    snprintf(diag_base_psi, sizeof(diag_base_psi), "data/Psiprofile_Ni%d_Ns%d.dat", Nintegration_diag, Nspline_diag_base);
                    get_suffixed_filename(diag_base_psi, 1, diag_fname_psi, sizeof(diag_fname_psi));
                    FILE *fp_diag_psi = fopen(diag_fname_psi, "wb");
                    if (fp_diag_psi) {
                        for (double r_write_diag = 0.0; r_write_diag < radius_diag_arr[num_points_diag - 1]; r_write_diag += rmax_diag / 900.0) {
                             if (r_write_diag >= radius_diag_arr[0]) {
                                fprintf_bin(fp_diag_psi, "%f %f\n", r_write_diag, evaluatespline(splinePsi_diag, Psiinterp_diag, r_write_diag));
                             }
                        }
                        fclose(fp_diag_psi);
                        log_message("DEBUG", "Wrote diagnostic file: %s", diag_fname_psi);
                    } else {
                        log_message("ERROR", "Failed to open diagnostic file: %s", diag_fname_psi);
                    }

                    // --- 4. r(Psi) spline for NFW Diagnostic ---
                    struct RrPsiPair *temp_pairs_npsi_nfw_diag = (struct RrPsiPair *)malloc(num_points_diag * sizeof(struct RrPsiPair));
                    if(!temp_pairs_npsi_nfw_diag) {
                        log_message("ERROR", "Failed to allocate sorting pairs for NFW diagnostic r(Psi)");
                        goto cleanup_diag_iteration;
                    }
                    for(int k_sort = 0; k_sort < num_points_diag; ++k_sort) {
                        temp_pairs_npsi_nfw_diag[k_sort].rr = nPsivalues_diag_arr[k_sort];
                        temp_pairs_npsi_nfw_diag[k_sort].psi = radius_for_rofPsi_diag_arr[k_sort];
                    }
                    qsort(temp_pairs_npsi_nfw_diag, num_points_diag, sizeof(struct RrPsiPair), compare_by_rr);
                    for(int k_sort = 0; k_sort < num_points_diag; ++k_sort) {
                        nPsivalues_diag_arr[k_sort] = temp_pairs_npsi_nfw_diag[k_sort].rr;
                        radius_for_rofPsi_diag_arr[k_sort] = temp_pairs_npsi_nfw_diag[k_sort].psi;
                    }
                    free(temp_pairs_npsi_nfw_diag);

                    rofPsiinterp_diag = gsl_interp_accel_alloc();
                    splinerofPsi_diag = gsl_spline_alloc(gsl_interp_cspline, num_points_diag);
                    if(!check_strict_monotonicity(nPsivalues_diag_arr, num_points_diag, "nPsivalues_diag_arr (NFW_DIAG)")) {
                        log_message("ERROR", "NFW diagnostic nPsivalues not monotonic");
                        goto cleanup_diag_iteration;
                    }
                    gsl_spline_init(splinerofPsi_diag, nPsivalues_diag_arr, radius_for_rofPsi_diag_arr, num_points_diag);


                    // --- 5. I(E) spline (fofEinterp_diag) for NFW Diagnostic ---
                    double Psimin_diag = Psivalues_diag_arr[num_points_diag - 1];
                    double Psimax_diag = Psivalues_diag_arr[0];
                    if (Psimax_diag <= Psimin_diag) {
                        log_message("ERROR", "NFW diagnostic potential not monotonic: Psimax=%e <= Psimin=%e", Psimax_diag, Psimin_diag);
                        goto cleanup_diag_iteration;
                    }

                    Evalues_diag_arr = (double *)malloc((num_points_diag + 1) * sizeof(double));
                    innerintegrandvalues_diag_arr = (double *)malloc((num_points_diag + 1) * sizeof(double));
                    if(!Evalues_diag_arr || !innerintegrandvalues_diag_arr) {
                        log_message("ERROR", "Failed to allocate arrays for NFW diagnostic I(E)");
                        goto cleanup_diag_iteration;
                    }

                    Evalues_diag_arr[0] = Psimin_diag;
                    innerintegrandvalues_diag_arr[0] = 0.0;
                    gsl_function F_fE_nfw_diag;
                    F_fE_nfw_diag.function = &fEintegrand_nfw;

                    for (int k = 1; k <= num_points_diag; k++) {
                        double calE_diag = Psimin_diag + (Psimax_diag - Psimin_diag) * ((double)k) / ((double)num_points_diag);
                        if (k > 0 && calE_diag <= Evalues_diag_arr[k-1]) {
                            calE_diag = Evalues_diag_arr[k-1] + DBL_EPSILON * fabs(Evalues_diag_arr[k-1]) + DBL_MIN;
                        }
                        Evalues_diag_arr[k] = calE_diag;

                        fE_integrand_params_NFW_t params_fE_nfw_diag = {
                            calE_diag, splinerofPsi_diag, rofPsiinterp_diag,
                            splinemass_diag, enclosedmass_diag, G_CONST,
                            current_diag_rc, nfw_params_diag[2], current_diag_falloff_C,
                            Psimin_diag, Psimax_diag
                        };
                        F_fE_nfw_diag.params = &params_fE_nfw_diag;
                        double t_upper_diag = sqrt(fmax(0.0, calE_diag - Psimin_diag));
                        double t_lower_diag = (t_upper_diag > 1e-9) ? t_upper_diag / 1.0e4 : 0.0;
                        if (t_lower_diag >= t_upper_diag - 1e-12) {
                            result_diag = 0.0;
                        } else {
                            gsl_integration_qag(&F_fE_nfw_diag, t_lower_diag, t_upper_diag, 1e-8, 1e-8,
                                                Nintegration_diag, GSL_INTEG_GAUSS61, w_diag, &result_diag, &error_diag);
                        }
                        innerintegrandvalues_diag_arr[k] = result_diag;
                    }

                    fofEacc_diag = gsl_interp_accel_alloc();
                    fofEinterp_diag = gsl_interp_alloc(gsl_interp_linear, num_points_diag + 1);
                    if(!check_strict_monotonicity(Evalues_diag_arr, num_points_diag + 1, "Evalues_diag_arr (NFW_DIAG)")) {
                        log_message("ERROR", "NFW diagnostic Evalues not monotonic");
                        goto cleanup_diag_iteration;
                    }
                    gsl_interp_init(fofEinterp_diag, Evalues_diag_arr, innerintegrandvalues_diag_arr, num_points_diag + 1);

                    // Write f(E) diagnostic file
                    char diag_fname_fofe[256];
                    char diag_base_fofe[128];
                    snprintf(diag_base_fofe, sizeof(diag_base_fofe), "data/f_of_E_Ni%d_Ns%d.dat", Nintegration_diag, Nspline_diag_base);
                    get_suffixed_filename(diag_base_fofe, 1, diag_fname_fofe, sizeof(diag_fname_fofe));
                    FILE *fp_diag_fofe = fopen(diag_fname_fofe, "wb");
                    if (fp_diag_fofe) {
                        for (int k = 0; k <= num_points_diag; k++) {
                            double E_diag = Evalues_diag_arr[k];
                            double deriv_diag = 0.0;
                            if (E_diag > Evalues_diag_arr[0] && E_diag < Evalues_diag_arr[num_points_diag]) {
                                 deriv_diag = gsl_interp_eval_deriv(fofEinterp_diag, Evalues_diag_arr, innerintegrandvalues_diag_arr, E_diag, fofEacc_diag);
                            }
                            double fE_val_diag = fabs(deriv_diag) / (sqrt(8.0) * PI * PI);
                            if (!isfinite(fE_val_diag)) fE_val_diag = 0.0;
                            fprintf_bin(fp_diag_fofe, "%f %f\n", E_diag, fE_val_diag);
                        }
                        fclose(fp_diag_fofe);
                        log_message("DEBUG", "Wrote diagnostic file: %s", diag_fname_fofe);
                    } else {
                        log_message("ERROR", "Failed to open diagnostic file: %s", diag_fname_fofe);
                    }

                    // Write NFW diagnostic integrand file
                    char diag_fname_integrand[256];
                    char diag_base_integrand[128];
                    snprintf(diag_base_integrand, sizeof(diag_base_integrand), "data/integrand_Ni%d_Ns%d.dat", Nintegration_diag, Nspline_diag_base);
                    get_suffixed_filename(diag_base_integrand, 1, diag_fname_integrand, sizeof(diag_fname_integrand));
                    FILE *fp_diag_int = fopen(diag_fname_integrand, "wb");
                    if (fp_diag_int) {
                        // Simple integrand convergence test (like Cored profile)
                        double calE_for_integrand = Psimax_diag;
                        fE_integrand_params_NFW_t params_int_nfw = {
                            calE_for_integrand, splinerofPsi_diag, rofPsiinterp_diag,
                            splinemass_diag, enclosedmass_diag, G_CONST,
                            current_diag_rc, nfw_params_diag[2], current_diag_falloff_C,
                            Psimin_diag, Psimax_diag
                        };

                        for (int k_int = 0; k_int < num_points_diag; k_int++) {
                            double t_diag_int = sqrt(fmax(0.0, calE_for_integrand - Psimin_diag)) * ((double)k_int) / ((double)num_points_diag);
                            fprintf_bin(fp_diag_int, "%f %f\n", t_diag_int, fEintegrand_nfw(t_diag_int, &params_int_nfw));
                        }

                        fclose(fp_diag_int);
                        log_message("DEBUG", "Wrote diagnostic file: %s", diag_fname_integrand);
                    } else {
                        log_message("ERROR", "Failed to open diagnostic file: %s", diag_fname_integrand);
                    }

                    // Write NFW diagnostic density profile
                    char diag_fname_dens[256];
                    char diag_base_dens[128];
                    snprintf(diag_base_dens, sizeof(diag_base_dens), "data/density_profile_Ni%d_Ns%d.dat", Nintegration_diag, Nspline_diag_base);
                    get_suffixed_filename(diag_base_dens, 1, diag_fname_dens, sizeof(diag_fname_dens));
                    FILE *fp_diag_dens = fopen(diag_fname_dens, "wb");
                    if (fp_diag_dens) {
                        for (int k = 0; k < num_points_diag; k++) {
                            double rr_k = radius_diag_arr[k];
                            double rs_k = rr_k / current_diag_rc;
                            double term_s_k = rs_k + 0.01;
                            if (term_s_k <= 1e-9) term_s_k = 1e-9;
                            double term_n_k = (1.0 + rs_k) * (1.0 + rs_k);
                            double term_c_base_k = rs_k / current_diag_falloff_C;
                            double term_c_k = 1.0 + pow(term_c_base_k, 10.0);
                            double rho_shape_k = (term_s_k < 1e-9 || term_n_k < 1e-9 || term_c_k < 1e-9) ? 0.0 : (1.0 / (term_s_k * term_n_k * term_c_k));
                            if (rr_k < 1e-6 && term_s_k < 1e-3 && rho_shape_k == 0.0) {
                               rho_shape_k = 1.0 / (term_s_k * term_n_k * term_c_k);
                            }
                            double rho_r_k = nfw_params_diag[2] * rho_shape_k;
                            fprintf_bin(fp_diag_dens, "%f %f\n", rr_k, rho_r_k);
                        }
                        fclose(fp_diag_dens);
                        log_message("DEBUG", "Wrote diagnostic file: %s", diag_fname_dens);
                    } else {
                        log_message("ERROR", "Failed to open diagnostic file: %s", diag_fname_dens);
                    }

                    // Write NFW diagnostic dPsi/dr
                    char diag_fname_dpsi[256];
                    char diag_base_dpsi[128];
                    snprintf(diag_base_dpsi, sizeof(diag_base_dpsi), "data/dpsi_dr_Ni%d_Ns%d.dat", Nintegration_diag, Nspline_diag_base);
                    get_suffixed_filename(diag_base_dpsi, 1, diag_fname_dpsi, sizeof(diag_fname_dpsi));
                    FILE *fp_diag_dpsi = fopen(diag_fname_dpsi, "wb");
                    if (fp_diag_dpsi) {
                        for (int k = 0; k < num_points_diag; k++) {
                            double rr_k = radius_diag_arr[k];
                            if (rr_k > 1e-9) {
                                double Menc_k = gsl_spline_eval(splinemass_diag, rr_k, enclosedmass_diag);
                                double dpsidr_k = -(G_CONST * Menc_k) / (rr_k * rr_k);
                                fprintf_bin(fp_diag_dpsi, "%f %f\n", rr_k, dpsidr_k);
                            }
                        }
                        fclose(fp_diag_dpsi);
                        log_message("DEBUG", "Wrote diagnostic file: %s", diag_fname_dpsi);
                    } else {
                        log_message("ERROR", "Failed to open diagnostic file: %s", diag_fname_dpsi);
                    }

                    // Write NFW diagnostic drho/dPsi
                    char diag_fname_drhodpsi[256];
                    char diag_base_drhodpsi[128];
                    snprintf(diag_base_drhodpsi, sizeof(diag_base_drhodpsi), "data/drho_dpsi_Ni%d_Ns%d.dat", Nintegration_diag, Nspline_diag_base);
                    get_suffixed_filename(diag_base_drhodpsi, 1, diag_fname_drhodpsi, sizeof(diag_fname_drhodpsi));
                    FILE *fp_diag_drhodpsi = fopen(diag_fname_drhodpsi, "wb");
                    if (fp_diag_drhodpsi) {
                        for (int k = 1; k < num_points_diag - 1; k++) {
                            double rr_k = radius_diag_arr[k];
                            if (rr_k <= 1e-9) continue;

                            double drhodr_val_k = drhodr_profile_nfwcutoff(rr_k, current_diag_rc, nfw_params_diag[2], current_diag_falloff_C);
                            double Menc_k = gsl_spline_eval(splinemass_diag, rr_k, enclosedmass_diag);
                            double dPsidr_mag_k = (G_CONST * Menc_k) / (rr_k * rr_k);

                            if (fabs(dPsidr_mag_k) > 1e-30) {
                                double Psi_val_k = evaluatespline(splinePsi_diag, Psiinterp_diag, rr_k);
                                double drho_dPsi_val_k = drhodr_val_k / dPsidr_mag_k;
                                fprintf_bin(fp_diag_drhodpsi, "%f %f\n", Psi_val_k, drho_dPsi_val_k);
                            }
                        }
                        fclose(fp_diag_drhodpsi);
                        log_message("DEBUG", "Wrote diagnostic file: %s", diag_fname_drhodpsi);
                    } else {
                        log_message("ERROR", "Failed to open diagnostic file: %s", diag_fname_drhodpsi);
                    }

                    log_message("DEBUG", "Finished diagnostic iteration %d, spline_base=%d", Nintegration_diag, Nspline_diag_base);

                    // --- Cleanup for NFW Diagnostic Iteration ---
cleanup_diag_iteration:
                    gsl_integration_workspace_free(w_diag);
                    if(splinemass_diag) gsl_spline_free(splinemass_diag);
                    if(enclosedmass_diag) gsl_interp_accel_free(enclosedmass_diag);
                    if(splinePsi_diag) gsl_spline_free(splinePsi_diag);
                    if(Psiinterp_diag) gsl_interp_accel_free(Psiinterp_diag);
                    if(splinerofPsi_diag) gsl_spline_free(splinerofPsi_diag);
                    if(rofPsiinterp_diag) gsl_interp_accel_free(rofPsiinterp_diag);
                    if(fofEinterp_diag) gsl_interp_free(fofEinterp_diag);
                    if(fofEacc_diag) gsl_interp_accel_free(fofEacc_diag);
                    free(mass_diag_arr);
                    free(radius_diag_arr);
                    free(radius_for_rofPsi_diag_arr);
                    free(Psivalues_diag_arr);
                    free(nPsivalues_diag_arr);
                    free(Evalues_diag_arr);
                    free(innerintegrandvalues_diag_arr);
                } // end Nspline_diag_base loop
            } // end Nintegration_diag loop
            log_message("INFO", "NFW DIAGNOSTIC LOOP: Completed convergence tests for NFW profile.");
        } // end if(g_doDebug) for NFW diagnostic loop

        // NFW PROFILE IC GENERATION PATHWAY

        /**
         * @brief NFW-specific theoretical calculation for initial conditions.
         * @details Calculates mass profile, potential, and distribution function
         *          for the NFW-like profile with power-law cutoff.
         */

        // Re-establish main NFW calc parameters for the main calculation
        double current_profile_rc = g_nfw_profile_rc;
        double current_profile_halo_mass = g_nfw_profile_halo_mass;
        double current_profile_rmax_norm_factor = g_nfw_profile_rmax_norm_factor;
        double current_profile_falloff_C = g_nfw_profile_falloff_factor;

        num_points = 100000; // Main NFW calculation's num_points
        int num_maxv2f = 1000; // Resolution for velocity envelope calculation only
        rmax = current_profile_rmax_norm_factor * current_profile_rc; // Main NFW rmax

        // Local variables for NFW pathway
        double nfw_result, nfw_error;
        double nfw_calE;
        int i_nfw;

        // Set NFW-specific parameters from generalized profile parameters
        g_nfw_profile_rc = g_scale_radius_param;
        g_nfw_profile_halo_mass = g_halo_mass_param;
        g_nfw_profile_rmax_norm_factor = g_cutoff_factor_param;

        // Update current_profile_* variables with new values
        current_profile_rc = g_nfw_profile_rc;
        current_profile_halo_mass = g_nfw_profile_halo_mass;
        current_profile_rmax_norm_factor = g_nfw_profile_rmax_norm_factor;
        current_profile_falloff_C = g_nfw_profile_falloff_factor;

        // Set numerical parameters for NFW
        num_points = 100000;  // Conservative default for NFW
        rmax = current_profile_rmax_norm_factor * current_profile_rc;

        // Allocate GSL workspace
        w = gsl_integration_workspace_alloc(1000);
        if (!w) {
            fprintf(stderr, "NFW_PATH: Failed to allocate GSL workspace\n");
            CLEAN_EXIT(1);
        }

        // Prepare mass integrand function
        gsl_function F_nfw_calc;
        F_nfw_calc.function = &massintegrand_profile_nfwcutoff;

        // Prepare parameters for NFW mass integrand: [rc, halo_mass, nt_nfw, falloff_factor]
        double nfw_params[4]; // Parameters for NFW mass integrand function
        nfw_params[0] = current_profile_rc;
        nfw_params[1] = current_profile_halo_mass; // Target total halo mass
        nfw_params[2] = 1.0; // Initial guess for nt_nfw normalization scaler
        nfw_params[3] = current_profile_falloff_C; // Falloff transition factor C
        F_nfw_calc.params = nfw_params;

        // Calculate normalization for NFW profile
        int status_norm = gsl_integration_qag(&F_nfw_calc, 0.0, rmax, 1e-12, 1e-12, 1000,
                            GSL_INTEG_GAUSS51, w, &nfw_result, &nfw_error);
        normalization = nfw_result;

        if (g_doDebug) {
            log_message("DEBUG", "Initial normalization integral (int r^2 * rho_guess dr):");
            log_message("DEBUG", "  rmax_for_norm_integral = %.3e kpc (factor=%.1f * rc=%.3f)", rmax, current_profile_rmax_norm_factor, current_profile_rc);
            log_message("DEBUG", "  GSL QAG status for norm_integral: %s", gsl_strerror(status_norm));
            log_message("DEBUG", "  Raw norm_integral_result (nfw_result for norm) = %.6e", nfw_result);
            log_message("DEBUG", "  Raw norm_integral_error_est = %.6e", nfw_error);
            log_message("DEBUG", "  Final 'normalization' variable = %.6e", normalization);
        }

        if (normalization <= 1e-30) {
            fprintf(stderr, "NFW_PATH: Normalization is zero or negative (%.3e). Exiting.\n", normalization);
            CLEAN_EXIT(1);
        }

        // Update nt_nfw with proper normalization
        nfw_params[2] = current_profile_halo_mass / (4.0 * M_PI * normalization);

        if (g_doDebug) {
            log_message("DEBUG", "Calculated nt_nfw (density scale factor):");
            log_message("DEBUG", "  current_profile_halo_mass (target M_total) = %.3e Msun", current_profile_halo_mass);
            log_message("DEBUG", "  nt_nfw = %.3e / (4pi * %.3e) = %.6e", current_profile_halo_mass, normalization, nfw_params[2]);
            if (!isfinite(nfw_params[2]) || (fabs(nfw_params[2]) < 1e-100 && fabs(nfw_params[2]) > 0)) {
                log_message("WARNING", "nt_nfw is NaN, Inf, or extremely small/large: %.6e", nfw_params[2]);
            }
        }

        log_message("INFO", "NFW Profile: RC=%.3f kpc, Halo Mass=%.3e Msun, Rmax_norm_calc=%.3f kpc, nt_nfw_scaler=%.6e",
               current_profile_rc, current_profile_halo_mass, rmax, nfw_params[2]);

        /**
         * @brief Calculate mass profile M(r) for NFW.
         */
        mass = (double *)malloc(num_points * sizeof(double));
        radius = (double *)malloc(num_points * sizeof(double));
        radius_monotonic_grid_nfw = (double *)malloc(num_points * sizeof(double));
        if (!mass || !radius || !radius_monotonic_grid_nfw) {
            fprintf(stderr, "NFW_PATH: Failed to allocate mass/radius arrays\n");
            CLEAN_EXIT(1);
        }

        mass[0] = 0.0;
        radius[0] = 0.0;                                  // For the y-values of r(Psi) spline later
        radius_monotonic_grid_nfw[0] = 0.0;             // For x-axes of M(r), Psi(r), maxv2f(r)


        for (i_nfw = 1; i_nfw < num_points; i_nfw++) {
            double r_current = (double)i_nfw * rmax / (num_points - 1);
            if (i_nfw == num_points - 1) r_current = rmax; // Ensure exact endpoint

            gsl_integration_qag(&F_nfw_calc, 0.0, r_current, 1e-12, 1e-12,
                                1000, GSL_INTEG_GAUSS51, w, &nfw_result, &nfw_error);
            mass[i_nfw] = 4.0 * M_PI * nfw_result;
            radius[i_nfw] = r_current; // This 'radius' array will be sorted with nPsivalues
            radius_monotonic_grid_nfw[i_nfw] = r_current; // This 'radius_monotonic_grid_nfw' stays sorted by r
        }

        if (g_doDebug) {
            log_message("DEBUG", "M(r) spline data summary (num_points=%d):", num_points);
            log_message("DEBUG", "  Target M_total for sampling = %.3e Msun", current_profile_halo_mass);
            log_message("DEBUG", "  nt_nfw used for M(r) calcs = %.3e", nfw_params[2]);
            log_message("DEBUG", "  rmax for M(r) array = %.3e kpc", rmax);
            if (num_points > 0) {
                log_message("DEBUG", "  Final Mass at rmax (radius[num_points-1]=%.3e kpc): %.3e Msun", radius[num_points-1], mass[num_points-1]);
                if (fabs(mass[num_points-1] - current_profile_halo_mass) / current_profile_halo_mass > 0.1) {
                    log_message("WARNING", "Mass at rmax (%.3e) differs significantly from target halo mass (%.3e)!", mass[num_points-1], current_profile_halo_mass);
                }
            }
            log_message("DEBUG", "End of M(r) data summary.");
        }

        // Create mass spline
        enclosedmass = gsl_interp_accel_alloc();
        splinemass = gsl_spline_alloc(gsl_interp_cspline, num_points);
        if (!enclosedmass || !splinemass) {
            fprintf(stderr, "NFW_PATH: Failed to allocate mass spline\n");
            CLEAN_EXIT(1);
        }
        gsl_spline_init(splinemass, radius_monotonic_grid_nfw, mass, num_points);

        // Write generic mass profile for plotting script
        if (g_doDebug) {
            fp = fopen("data/massprofile.dat", "wb");
            if (fp) {
                for (double r_write = 0.0; r_write < radius_monotonic_grid_nfw[num_points-1]; r_write += rmax / 900.0) {
                    double mass_at_r = gsl_spline_eval(splinemass, r_write, enclosedmass);
                    fprintf(fp, "%e %e\n", r_write, mass_at_r);
                }
                fclose(fp);
            }
        }

        /**
         * @brief Calculate gravitational potential Psi(r) for NFW.
         */
        Psivalues = (double *)malloc(num_points * sizeof(double));
        nPsivalues = (double *)malloc(num_points * sizeof(double));
        if (!Psivalues || !nPsivalues) {
            fprintf(stderr, "NFW_PATH: Failed to allocate Psi arrays\n");
            CLEAN_EXIT(1);
        }

        // Prepare Psiintegrand parameters for NFW
        Psiintegrand_params psi_params_nfw;
        psi_params_nfw.massintegrand_func = &massintegrand_profile_nfwcutoff;
        psi_params_nfw.params_for_massintegrand = nfw_params;

        gsl_function F_for_psi_nfw;
        F_for_psi_nfw.function = &Psiintegrand;
        F_for_psi_nfw.params = &psi_params_nfw;

        for (i_nfw = 0; i_nfw < num_points; i_nfw++) {
            double r_current = radius[i_nfw];
            double r1_psi = fmax(r_current, current_profile_rc / 1000000.0);

            gsl_integration_qagiu(&F_for_psi_nfw, r1_psi, 1e-12, 1e-12,
                                  1000, w, &nfw_result, &nfw_error);
            double first_term_psi = G_CONST * gsl_spline_eval(splinemass, r1_psi, enclosedmass) / r1_psi;
            double second_term_psi = G_CONST * 4.0 * M_PI * nfw_result;
            if (r1_psi < 1e-9) first_term_psi = 0; // Avoid division by zero
            Psivalues[i_nfw] = (first_term_psi + second_term_psi);
            nPsivalues[i_nfw] = -Psivalues[i_nfw];
        }

        if (g_doDebug) {
            log_message("DEBUG", "Psi(r) and r(Psi) spline data summary (num_points=%d):", num_points);
            if (num_points > 1) {
                log_message("DEBUG", "  Psivalues[0] (Psimax candidate) = %.6e", Psivalues[0]);
                log_message("DEBUG", "  Psivalues[num_points-1] (Psimin candidate) = %.6e", Psivalues[num_points-1]);
                if (Psivalues[0] <= Psivalues[num_points-1]) {
                    log_message("WARNING", "Psi(r) may not be monotonic decreasing (Psivalues[0]=%.3e <= Psivalues[end]=%.3e)!", Psivalues[0], Psivalues[num_points-1]);
                }
            }
            log_message("DEBUG", "End of Psi(r) data summary.");
        }

        // Create Psi splines
        Psiinterp = gsl_interp_accel_alloc();
        splinePsi = gsl_spline_alloc(gsl_interp_cspline, num_points);
        if (!Psiinterp || !splinePsi) {
            fprintf(stderr, "NFW_PATH: Failed to allocate Psi spline\n");
            CLEAN_EXIT(1);
        }
        gsl_spline_init(splinePsi, radius_monotonic_grid_nfw, Psivalues, num_points);

        // Write generic Psi profile for plotting script
        if (g_doDebug) {
            fp = fopen("data/Psiprofile.dat", "wb");
            if (fp) {
                for (double r_write = 0.0; r_write < radius_monotonic_grid_nfw[num_points-1]; r_write += rmax / 900.0) {
                    double psi_at_r = gsl_spline_eval(splinePsi, r_write, Psiinterp);
                    fprintf(fp, "%e %e\n", r_write, psi_at_r);
                }
                fclose(fp);
            }
        }

        // Create r(Psi) spline
        rofPsiinterp = gsl_interp_accel_alloc();
        splinerofPsi = gsl_spline_alloc(gsl_interp_cspline, num_points);
        if (!rofPsiinterp || !splinerofPsi) {
            fprintf(stderr, "NFW_PATH: Failed to allocate r(Psi) spline\n");
            CLEAN_EXIT(1);
        }

        // Use temporary copies for r(Psi) spline to preserve the original radius grid
        double *nPsivalues_for_rPsi_spline = (double *)malloc(num_points * sizeof(double));
        double *radius_values_for_rPsi_spline = (double *)malloc(num_points * sizeof(double));

        if (!nPsivalues_for_rPsi_spline || !radius_values_for_rPsi_spline) {
            fprintf(stderr, "NFW_PATH: Failed to allocate temp arrays for r(Psi) spline data.\n");
            if (nPsivalues_for_rPsi_spline) free(nPsivalues_for_rPsi_spline);
            if (radius_values_for_rPsi_spline) free(radius_values_for_rPsi_spline);
            CLEAN_EXIT(1);
        }

        // Copy nPsivalues and the corresponding radius_monotonic_grid_nfw values
        memcpy(nPsivalues_for_rPsi_spline, nPsivalues, num_points * sizeof(double));
        memcpy(radius_values_for_rPsi_spline, radius_monotonic_grid_nfw, num_points * sizeof(double));

        // Sort nPsivalues_for_rPsi_spline and apply identical swaps to radius_values_for_rPsi_spline
        for (int k_sort = 0; k_sort < num_points - 1; k_sort++) {
            for (int j_sort = k_sort + 1; j_sort < num_points; j_sort++) {
                if (nPsivalues_for_rPsi_spline[k_sort] > nPsivalues_for_rPsi_spline[j_sort]) {
                    // Swap nPsivalues_for_rPsi_spline
                    double temp_npsi = nPsivalues_for_rPsi_spline[k_sort];
                    nPsivalues_for_rPsi_spline[k_sort] = nPsivalues_for_rPsi_spline[j_sort];
                    nPsivalues_for_rPsi_spline[j_sort] = temp_npsi;

                    // Swap corresponding radius_values_for_rPsi_spline
                    double temp_rad = radius_values_for_rPsi_spline[k_sort];
                    radius_values_for_rPsi_spline[k_sort] = radius_values_for_rPsi_spline[j_sort];
                    radius_values_for_rPsi_spline[j_sort] = temp_rad;
                }
            }
        }

        // Debug check for nPsivalues_for_rPsi_spline monotonicity
        if (g_doDebug) {
            int mono_violations_npsi = 0;
            for (int k_chk = 0; k_chk < num_points - 1; ++k_chk) {
                if (!(nPsivalues_for_rPsi_spline[k_chk+1] > nPsivalues_for_rPsi_spline[k_chk])) {
                    if (mono_violations_npsi < 5) log_message("DEBUG", "nPsivalues_for_rPsi_spline not strictly increasing at index %d", k_chk);
                    mono_violations_npsi++;
                }
            }
            if (mono_violations_npsi > 0) log_message("WARNING", "Total nPsivalues violations in r(Psi) spline: %d", mono_violations_npsi);
            else if (g_doDebug) log_message("DEBUG", "nPsivalues_for_rPsi_spline confirmed strictly monotonic for r(Psi) spline.");
        }

        // Initialize splinerofPsi with the sorted temporary arrays
        gsl_spline_init(splinerofPsi, nPsivalues_for_rPsi_spline, radius_values_for_rPsi_spline, num_points);

        // Free the temporary sorted copies
        free(nPsivalues_for_rPsi_spline);
        free(radius_values_for_rPsi_spline);

        /**
         * @brief Calculate f(E) distribution function for NFW using Eddington's formula.
         */
        Psimin = Psivalues[num_points - 1];
        Psimax = Psivalues[0];

        if (g_doDebug) {
            log_message("DEBUG", "Potential range for I(E) calculation: Psimin=%.6e, Psimax=%.6e", Psimin, Psimax);
            log_message("DEBUG", "Continuing with I(E) calculation.");
        }

        if (Psimax <= Psimin) {
            fprintf(stderr, "NFW_PATH: Potential not monotonic (Psimax=%.3e <= Psimin=%.3e)\n",
                    Psimax, Psimin);
            CLEAN_EXIT(1);
        }

        innerintegrandvalues = (double *)malloc((num_points + 1) * sizeof(double));
        Evalues = (double *)malloc((num_points + 1) * sizeof(double));
        if (!innerintegrandvalues || !Evalues) {
            fprintf(stderr, "NFW_PATH: Failed to allocate f(E) arrays\n");
            CLEAN_EXIT(1);
        }

        // NFW uses conservative tolerance for f(E) integral
        F_nfw_calc.function = &fEintegrand_nfw;

        innerintegrandvalues[0] = 0.0;
        Evalues[0] = Psimin;

        // GSL integration status tracking
        int status_fE_nfw_local;

        for (i_nfw = 1; i_nfw <= num_points; i_nfw++) {
            nfw_calE = Psimin + (Psimax - Psimin) * ((double)i_nfw) / ((double)num_points);

            // Enforce strict monotonicity for Evalues
            if (i_nfw > 0 && nfw_calE <= Evalues[i_nfw-1]) {
                // If current nfw_calE is not strictly greater than previous, add a tiny increment.
                // DBL_EPSILON for the scale of Evalues might be too small if Evalues are large.
                // A small fraction of the typical step, or a fixed small number relative to Evalues scale.
                double previous_E = Evalues[i_nfw-1];
                double ideal_step = (Psimax - Psimin) / (double)num_points;
                double increment = ideal_step * 1e-6; // Small fraction of an ideal step
                if (increment == 0.0) increment = DBL_MIN * fabs(previous_E) + DBL_MIN; // Absolute minimum if ideal step is zero
                if (increment == 0.0) increment = 1e-20; // Fallback if previous_E is zero

                nfw_calE = previous_E + increment;

                if (g_doDebug && (i_nfw <= 10 || i_nfw > num_points -10 || i_nfw % (num_points/50<1?1:num_points/50) == 0) ) { // Log adjustment sparsely
                     log_message("DEBUG","Adjusted Evalues[%d] from ideal %.17e to %.17e (prev E: %.17e)",
                            i_nfw, Psimin + (Psimax - Psimin) * ((double)i_nfw) / ((double)num_points),
                            nfw_calE, previous_E);
                }
            }

            // Create NFW-specific parameter structure

            fE_integrand_params_NFW_t params_for_fE_integrand_nfw = {
                nfw_calE,               // E_current_shell
                splinerofPsi,           // spline_r_of_Psi (this is r_of_nPsi from nPsivalues)
                rofPsiinterp,           // accel_r_of_Psi
                splinemass,             // spline_M_of_r (this is M(r) for NFW, from corrected NFMP.1)
                enclosedmass,           // accel_M_of_r
                G_CONST,                // const_G_universal
                current_profile_rc,     // profile_rc_const
                nfw_params[2],          // profile_nt_norm_const (this is the scaled nt_nfw)
                nfw_params[3],          // profile_falloff_C_const (falloff factor C)
                Psimin,                 // Psimin_global
                Psimax                  // Psimax_global
            };
            F_nfw_calc.params = &params_for_fE_integrand_nfw;

            double E_current_shell = nfw_calE; // E for which I(E) is being computed

            // Integration will be over t_prime = sqrt(E_shell - Psi_true)
            // As Psi_true goes from Psimin_global to E_shell, t_prime goes from sqrt(E_shell - Psimin_global) down to 0.
            // So, integrate t_prime from 0 to sqrt(E_shell - Psimin_global).

            double t_integration_upper_bound = sqrt(fmax(0.0, E_current_shell - Psimin));
            double t_integration_lower_bound;

            if (t_integration_upper_bound < 1e-9) { // If E_current_shell is very close to Psimin (or below)
                t_integration_lower_bound = 0.0;
                t_integration_upper_bound = 0.0;
            } else {
                // Set a very small, but strictly positive, lower bound relative to the upper bound,
                // or an absolute small number if t_upper_bound is itself very small.
                // This helps GSL avoid evaluating exactly at t=0 if there's a 1/t or 1/sqrt(t) type issue.
                // The term 1/sqrt(E-Psi) in d(rho)/d(Psi) / sqrt(E-Psi) becomes 1/t when Psi = E-t^2.
                // Our fEintegrand_nfw is 2 * d(rho)/d(Psi), so it does not have this explicit 1/t.
                // Using t_upper_bound / 1.0e4 scaling for numerical consistency.
                t_integration_lower_bound = t_integration_upper_bound / 1.0e4;
                // If t_integration_lower_bound becomes extremely small (e.g. < DBL_MIN), GSL might treat it as zero.
                // Ensure it's at least some representable small positive number if t_upper_bound is positive.
                if (t_integration_lower_bound == 0.0 && t_integration_upper_bound > 0.0) {
                    t_integration_lower_bound = DBL_EPSILON * t_integration_upper_bound; // Or just DBL_EPSILON if t_upper is very small
                    if (t_integration_lower_bound == 0.0) t_integration_lower_bound = 1e-20; // Absolute floor
                }
            }

            // Ensure lower bound is strictly less than upper bound for GSL
            if (t_integration_lower_bound >= t_integration_upper_bound - 1e-12) { // Adjusted epsilon
                t_integration_upper_bound = 0.0; // Force zero integration range
                t_integration_lower_bound = 0.0;
            }

            if (g_doDebug && (i_nfw <= 5 || i_nfw > num_points - 5 || i_nfw % (num_points/10 < 1 ? 1 : num_points/10) == 0) ) {
                log_message("DEBUG", "I(E) integral setup: E_shell=%.3e, Psimin=%.3e, integrating fEintegrand_nfw(t) from t_low=%.3e to t_high=%.3e",
                       E_current_shell, Psimin, t_integration_lower_bound, t_integration_upper_bound);
            }

            if (t_integration_upper_bound <= t_integration_lower_bound + 1e-10) { // If range is zero or too small
                nfw_result = 0.0;
                status_fE_nfw_local = GSL_SUCCESS;
                 if (g_doDebug && (i_nfw <= 5 || i_nfw > num_points - 5 || i_nfw % (num_points/10 < 1 ? 1 : num_points/10) == 0) ) {
                    log_message("DEBUG", "NFEFE_INTEGRAL_SETUP: Skipping t-integration, range invalid/tiny (t_high=%.3e, t_low=%.3e)", t_integration_upper_bound, t_integration_lower_bound);
                 }
            } else {
                status_fE_nfw_local = gsl_integration_qag(&F_nfw_calc, t_integration_lower_bound, t_integration_upper_bound,
                                    1e-8, 1e-8, 1000, GSL_INTEG_GAUSS61, // Using conservative GSL tolerances
                                    w, &nfw_result, &nfw_error);

                if (g_doDebug && (i_nfw <= 5 || i_nfw > num_points - 5 || i_nfw % (num_points/10 < 1 ? 1 : num_points/10) == 0) ) {
                    log_message("DEBUG", "NFEFE_INTEGRAL_RESULT: I(E=%.3e) = %.6e, error=%.3e, status=%s",
                           E_current_shell, nfw_result, nfw_error,
                           (status_fE_nfw_local == GSL_SUCCESS) ? "SUCCESS" : "ERROR");
                }
            }
            innerintegrandvalues[i_nfw] = nfw_result;
            Evalues[i_nfw] = nfw_calE;
        }


        // Create f(E) interpolation using cspline interpolation for NFW (smoother dI/dE)
        g_main_fofEinterp = gsl_interp_alloc(gsl_interp_cspline, num_points + 1);
        g_main_fofEacc = gsl_interp_accel_alloc();
        if (!g_main_fofEinterp || !g_main_fofEacc) {
            fprintf(stderr, "NFW_PATH: Failed to allocate f(E) interpolation\n");
            CLEAN_EXIT(1);
        }
        // ADD THIS BLOCK BEFORE gsl_interp_init:
        if (g_doDebug) {
            int monotonicity_violations = 0;
            log_message("DEBUG", "Checking Evalues for strict monotonicity (%d points) before I(E) spline init.", num_points + 1);
            // Evalues has num_points + 1 elements, indexed 0 to num_points.
            for (int chk_e = 0; chk_e < num_points; ++chk_e) { // Loop up to num_points-1 to check Evalues[chk_e+1] vs Evalues[chk_e]
                if (!(Evalues[chk_e+1] > Evalues[chk_e])) {
                    if (monotonicity_violations < 20) { // Print first few violations
                        fprintf(stderr, "  MONOTONICITY_VIOLATION_PRE_SPLINE: Evalues[%d]=%.17e, Evalues[%d]=%.17e (Diff: %.3e)\n",
                               chk_e, Evalues[chk_e], chk_e+1, Evalues[chk_e+1], Evalues[chk_e+1] - Evalues[chk_e]);
                    }
                    monotonicity_violations++;
                }
            }
            if (monotonicity_violations > 0) {
                fprintf(stderr, "  NFW_CRITICAL_SPLINE_INIT: Total Evalues monotonicity violations: %d. GSL interp_init will likely fail. Exiting.\n", monotonicity_violations);
                CLEAN_EXIT(1); // Add explicit exit if violations found.
            } else if (g_doDebug) { // Only log success if in debug mode
                log_message("DEBUG", "Evalues array confirmed strictly monotonic before I(E) spline init.");
            }
        }
        // END ADDED BLOCK

        gsl_interp_init(g_main_fofEinterp, Evalues, innerintegrandvalues, num_points + 1);

        gsl_integration_workspace_free(w);
        w = NULL;

        log_message("INFO", "NFW theoretical calculation for IC splines complete.");

        /**
         * @brief NFW Sample Generator - Generate particle positions and velocities.
         * @details Uses rejection sampling with the NFW density profile and f(E) distribution.
         */
        if (doReadInit) {
            // Read initial conditions from file
            printf("NFW_PATH: Reading initial conditions from %s...\n", readInitFilename);
            read_initial_conditions(particles, npts_initial, readInitFilename);
        } else if (!g_doRestart) {
            if (tidal_fraction > 0.0) {
                log_message("INFO", "NFW IC Gen: Initial particle count before stripping: %d", npts_initial);
            }
            log_message("INFO", "NFW IC Gen: Generating %d initial particle positions and velocities...", npts_initial);

            // Print overall Psimin, Psimax for NFW path once
            if (npts_initial > 0) { // Avoid printing if no particles
                if (Psimax <= Psimin) {
                }
            }

            /**
             * @brief Allocate memory for the particle data array.
             * @details 2D array particles[5][npts_initial] where:
             *          [0] = radius, [1] = velocity, [2] = angular momentum,
             *          [3] = particle ID, [4] = orientation (mu)
             */
            particles = (double **)malloc(5 * sizeof(double *));
            if (particles == NULL) {
                fprintf(stderr, "NFW_PATH: Memory allocation failed for particle array\n");
                CLEAN_EXIT(1);
            }
            for (i = 0; i < 5; i++) {
                particles[i] = (double *)malloc(npts_initial * sizeof(double));
                if (particles[i] == NULL) {
                    fprintf(stderr, "NFW_PATH: Memory allocation failed for particles[%d]\n", i);
                    CLEAN_EXIT(1);
                }
            }

            /**
             * @brief Calculate maximum velocity squared at each radius for rejection sampling.
             */
            double *maxv2f_nfw = (double *)malloc(num_maxv2f * sizeof(double));
            double *radius_maxv2f_nfw = (double *)malloc(num_maxv2f * sizeof(double));
            if (!maxv2f_nfw || !radius_maxv2f_nfw) {
                fprintf(stderr, "NFW_PATH: Failed to allocate maxv2f arrays\n");
                CLEAN_EXIT(1);
            }

            double nfw_vel, nfw_ratio, nfw_Psir, nfw_mu, nfw_maxv, nfw_maxvalue;

            // Create spline for r(M) - radius as function of enclosed mass
            gsl_interp_accel *rofMaccel_nfw = gsl_interp_accel_alloc();
            gsl_spline *splinerofM_nfw = gsl_spline_alloc(gsl_interp_cspline, num_points);
            if (!rofMaccel_nfw || !splinerofM_nfw) {
                fprintf(stderr, "NFW_PATH: Failed to allocate r(M) spline\n");
                CLEAN_EXIT(1);
            }
            // ADD THIS DIAGNOSTIC BLOCK:
        if (g_doDebug) {
            int monotonicity_violations_mass_spline = 0;
            log_message("DEBUG", "Checking mass array for strict monotonicity (size: %d)", num_points);
            // 'mass' array has num_points elements. Loop up to num_points-2 to check mass[chk+1] vs mass[chk].
            if (num_points >= 2) { // Need at least 2 points to check monotonicity
                for (int chk_m = 0; chk_m < num_points - 1; ++chk_m) {
                    if (!(mass[chk_m+1] > mass[chk_m])) {
                        if (monotonicity_violations_mass_spline < 20) {
                            fprintf(stderr, "  Mass array monotonicity violation: mass[%d]=%.17e >= mass[%d]=%.17e (Diff: %.3e)\n",
                                   chk_m, mass[chk_m], chk_m+1, mass[chk_m+1], mass[chk_m+1] - mass[chk_m]);
                        }
                        monotonicity_violations_mass_spline++;
                    }
                }
            }
            if (monotonicity_violations_mass_spline > 0) {
                fprintf(stderr, "  NFW_CRITICAL_MONO_CHECK_ROFM: Total 'mass' array monotonicity violations: %d. gsl_spline_init for splinerofM_nfw will fail.\n", monotonicity_violations_mass_spline);
                 if (monotonicity_violations_mass_spline > 20) fprintf(stderr, "  (Further violations suppressed)\n");
            } else {
                log_message("DEBUG", "Mass array confirmed strictly monotonic.");
            }
        }
        // END ADDED DIAGNOSTIC BLOCK

        gsl_spline_init(splinerofM_nfw, mass, radius, num_points);

            // Calculate max v^2 * f(E) at each radius
            radius_maxv2f_nfw[0] = 0.0;
            for (int i_r_nfw = 1; i_r_nfw < num_maxv2f; i_r_nfw++) {
                nfw_maxvalue = 0.0;
                double r_maxv2f = (double)i_r_nfw * rmax / (num_maxv2f - 1);
                radius_maxv2f_nfw[i_r_nfw] = r_maxv2f;
                nfw_Psir = evaluatespline(splinePsi, Psiinterp, r_maxv2f);
                nfw_maxv = sqrt(2.0 * (nfw_Psir - Psimin));

                // Find maximum of v^2 * dI/dE over velocity range
                for (int j_v_nfw = 1; j_v_nfw < num_maxv2f - 2; j_v_nfw++) {
                    nfw_vel = nfw_maxv * ((double)j_v_nfw) / ((double)num_maxv2f);
                    double E_test_nfw = nfw_Psir - 0.5 * nfw_vel * nfw_vel;
                    double currentvalue_nfw = 0.0;

                    if (E_test_nfw >= Psimin && E_test_nfw <= Psimax) {
                        currentvalue_nfw = nfw_vel * nfw_vel *
                            fabs(gsl_interp_eval_deriv(g_main_fofEinterp, Evalues,
                                                       innerintegrandvalues, E_test_nfw, g_main_fofEacc));
                    }
                    if (isfinite(currentvalue_nfw) && currentvalue_nfw > nfw_maxvalue) {
                        nfw_maxvalue = currentvalue_nfw;
                    }
                }
                maxv2f_nfw[i_r_nfw] = nfw_maxvalue;
            }

            // Extrapolate for r=0
            if (num_maxv2f >= 3) {
                maxv2f_nfw[0] = 2.0 * maxv2f_nfw[1] - maxv2f_nfw[2];
                if (maxv2f_nfw[0] < 0) maxv2f_nfw[0] = 0;
            } else {
                maxv2f_nfw[0] = maxv2f_nfw[1];
            }

            // Create spline for max v^2 * f(E)
            gsl_interp_accel *maxv2faccel_nfw = gsl_interp_accel_alloc();
            gsl_spline *splinemaxv2f_nfw = gsl_spline_alloc(gsl_interp_cspline, num_maxv2f);
            if (!maxv2faccel_nfw || !splinemaxv2f_nfw) {
                fprintf(stderr, "NFW_PATH: Failed to allocate maxv2f spline\n");
                CLEAN_EXIT(1);
            }
            // ADD THIS DIAGNOSTIC BLOCK:
        if (g_doDebug) {
            int monotonicity_violations_rad_spline = 0;
            log_message("DEBUG", "Checking radius array for strict monotonicity (size: %d)", num_points);
            // 'radius' array has num_points elements. Loop up to num_points-2.
            if (num_points >= 2) {
                for (int chk_r = 0; chk_r < num_points - 1; ++chk_r) {
                    if (!(radius[chk_r+1] > radius[chk_r])) {
                        if (monotonicity_violations_rad_spline < 20) {
                            fprintf(stderr, "  Radius array monotonicity violation: radius[%d]=%.17e >= radius[%d]=%.17e (Diff: %.3e)\n",
                                   chk_r, radius[chk_r], chk_r+1, radius[chk_r+1], radius[chk_r+1] - radius[chk_r]);
                        }
                        monotonicity_violations_rad_spline++;
                    }
                }
            }
            if (monotonicity_violations_rad_spline > 0) {
                fprintf(stderr, "  NFW_CRITICAL_MONO_CHECK_MAXV2F: Total 'radius' array monotonicity violations: %d. gsl_spline_init for splinemaxv2f_nfw will fail.\n", monotonicity_violations_rad_spline);
                if (monotonicity_violations_rad_spline > 20) fprintf(stderr, "  (Further violations suppressed)\n");
            } else {
                log_message("DEBUG", "Radius array confirmed strictly monotonic.");
            }
        }
        // END ADDED DIAGNOSTIC BLOCK

        gsl_spline_init(splinemaxv2f_nfw, radius_maxv2f_nfw, maxv2f_nfw, num_maxv2f);

            /**
             * @brief Generate particles using rejection sampling.
             */
            for (int k_nfw = 0; k_nfw < npts_initial; k_nfw++) {
                if (k_nfw < 5 || k_nfw % (npts_initial / 10 < 1 ? 1 : npts_initial/10) == 0) { // Log for first few & periodically
                    fflush(stdout);
                }

                // Sample radius from mass distribution
                double mass_frac_sample_nfw = gsl_rng_uniform(g_rng) * 0.999999;
                double mass_sample_nfw = mass_frac_sample_nfw * current_profile_halo_mass;
                particles[0][k_nfw] = evaluatespline(splinerofM_nfw, rofMaccel_nfw, mass_sample_nfw);

                if (k_nfw < 5 || k_nfw % (npts_initial / 10 < 1 ? 1 : npts_initial/10) == 0) {
                }

                nfw_maxvalue = evaluatespline(splinemaxv2f_nfw, maxv2faccel_nfw, particles[0][k_nfw]);
                nfw_Psir = evaluatespline(splinePsi, Psiinterp, particles[0][k_nfw]);

                // Check for problematic values
                if (!isfinite(nfw_Psir)) {
                    if (k_nfw < 5 || k_nfw % (npts_initial / 10 < 1 ? 1 : npts_initial/10) == 0) {
                    }
                    particles[1][k_nfw] = 0.0; // Assign zero velocity
                    nfw_mu = (2.0 * gsl_rng_uniform(g_rng) - 1.0);
                    particles[2][k_nfw] = 0.0; // L = 0 since v = 0
                    particles[4][k_nfw] = nfw_mu;
                    particles[3][k_nfw] = (double)k_nfw;
                    continue;
                }

                if (nfw_Psir <= Psimin + 1e-9 * fabs(Psimin)) { // Check if Psir is too close to Psimin
                    if (k_nfw < 5 || k_nfw % (npts_initial / 10 < 1 ? 1 : npts_initial/10) == 0) {
                    }
                    particles[1][k_nfw] = 0.0; // No kinetic energy possible
                } else {
                    // Sample velocity using rejection method
                    nfw_maxv = sqrt(fmax(0.0, 2.0 * (nfw_Psir - Psimin)));
                    if (!isfinite(nfw_maxv) || nfw_maxv < 1e-9) {
                        particles[1][k_nfw] = 0.0;
                    } else {

                        if (!isfinite(nfw_maxvalue) || nfw_maxvalue <= 1e-30) { // If envelope is effectively zero
                            particles[1][k_nfw] = 0.0;
                        } else {
                            // Velocity Rejection Sampling Loop
                            int vflag_nfw = 0;
                            int v_trials_nfw = 0;

                    while (vflag_nfw == 0 && v_trials_nfw < 20000) {
                        v_trials_nfw++;
                        nfw_vel = gsl_rng_uniform(g_rng) * nfw_maxv;
                        double E_test_nfw = nfw_Psir - 0.5 * nfw_vel * nfw_vel;
                        double target_func_val_nfw = 0.0;

                        double deriv_val_dIdE = 0.0;
                        if (E_test_nfw >= Psimin - 1e-9*fabs(Psimin) && E_test_nfw <= Psimax + 1e-9*fabs(Psimax)) { // Looser check for spline domain
                            deriv_val_dIdE = gsl_interp_eval_deriv(g_main_fofEinterp, Evalues,
                                                                   innerintegrandvalues, E_test_nfw, g_main_fofEacc);
                        }

                        // Add diagnostic for dI/dE values
                        if (g_doDebug && v_trials_nfw <= 2 && k_nfw < 5) { // Only for very first few trials of first few particles
                        }

                        target_func_val_nfw = nfw_vel * nfw_vel * fabs(deriv_val_dIdE);
                        if (!isfinite(target_func_val_nfw) || target_func_val_nfw < 0) target_func_val_nfw = 0.0; // Ensure non-negative

                        nfw_ratio = target_func_val_nfw / nfw_maxvalue; // maxvalue should be >0 here
                        if (nfw_ratio < 0) nfw_ratio = 0;
                        if (nfw_ratio > 1.001) { // If ratio is slightly > 1 due to numerics
                            nfw_ratio = 1.0;
                        }

                        if ((k_nfw < 2 && v_trials_nfw < 5) || (v_trials_nfw % 5000 == 0 && v_trials_nfw > 0) ) {
                        }

                        // Enhanced high trial count diagnostics
                        if (g_doDebug && (v_trials_nfw % 4000 == 0 && v_trials_nfw > 0)) {
                        }

                        // Diagnostic for zero dI/dE in valid energy range
                        if (g_doDebug && fabs(deriv_val_dIdE) < 1e-20 && (E_test_nfw > Psimin + 1e-6*fabs(Psimin) && E_test_nfw < Psimax - 1e-6*fabs(Psimax)) && (v_trials_nfw % 100 == 0) && v_trials_nfw > 0 && k_nfw < 100) {
                        }

                        if (gsl_rng_uniform(g_rng) < nfw_ratio) {
                            particles[1][k_nfw] = nfw_vel;
                            vflag_nfw = 1;
                        }
                    }
                    if (!vflag_nfw) {
                        particles[1][k_nfw] = 0.0; // Failed to find velocity
                    }
                        } // End else (maxvalue_envelope is finite and positive)
                    } // End else (maxv is finite and positive)
                } // End else (Psir > Psimin)

                // Sample angular momentum direction
                nfw_mu = 2.0 * gsl_rng_uniform(g_rng) - 1.0;
                // Ensure L is non-negative and well-defined even if particles[1][k_nfw] (velocity magnitude) is 0
                double L_val_nfw = 0.0;
                if (particles[1][k_nfw] > 1e-9) { // If velocity is non-zero
                    L_val_nfw = particles[1][k_nfw] * particles[0][k_nfw] * sqrt(fmax(0.0, 1.0 - nfw_mu * nfw_mu));
                }
                particles[2][k_nfw] = L_val_nfw;
                particles[3][k_nfw] = (double)k_nfw; // Particle ID
                particles[4][k_nfw] = nfw_mu;        // Orientation

                if (k_nfw < 5 || k_nfw % (npts_initial / 10 < 1 ? 1 : npts_initial/10) == 0) {
                }
            }

            // Clean up NFW sample generator allocations
            gsl_spline_free(splinerofM_nfw);
            gsl_interp_accel_free(rofMaccel_nfw);
            gsl_spline_free(splinemaxv2f_nfw);
            gsl_interp_accel_free(maxv2faccel_nfw);
            free(maxv2f_nfw);
            free(radius_maxv2f_nfw);

            log_message("INFO", "NFW IC Gen: Successfully generated %d particles.", npts_initial);
        } // End NFW sample generator


    } else { // Default: Use Cored Plummer-like Profile (Original Pathway)
        log_message("INFO", "Starting IC generation using Cored Plummer-like profile pathway.");
        log_message("INFO", "Generating Initial Conditions using Cored Plummer-like profile (original method)...");


        // ORIGINAL CORED PLUMMER-LIKE IC GENERATION PATHWAY
        // This is the entire block from nsphere.c.main.may20_1554.txt starting with
        // its "Theoretical Calculation Loop (Diagnostic)" down to the end of its
        // "SAMPLE GENERATOR" block.
        // It uses RC, HALO_MASS macros, its original massintegrand/drhodr,
        // and its original GSL settings.

        /** @brief Arrays defining integration and spline point counts for theoretical calculations. */
        int integration_points_array[2] = {1000, 10000};
        int spline_points_array[2] = {1000, 10000};

    /**
     * @brief Theoretical Calculation Loop (Eddington's Formula - Multiple Params).
     * @details Calculates theoretical profiles (mass, potential, f(E), density) based on the
     *          assumed initial density profile using Eddington's formula and GSL integration/splines.
     *          This loop iterates through different numbers of integration points (`Nintegration`)
     *          and spline points (`Nspline`) from the arrays above to generate reference files
     *          (e.g., massprofile_NiX_NsY.dat). These files are generated for comparison/validation
     *          but are *not* directly used in the main simulation timestepping or IC generation.
     *          They demonstrate the calculation process with varying numerical precision settings.
     */
    {
        for (int ii_ip = 0; ii_ip < 2; ii_ip++) // Loop over Nintegration values
        {
            for (int ii_sp = 0; ii_sp < 2; ii_sp++) // Loop over Nspline values
            {
                int Nintegration = integration_points_array[ii_ip];
                int Nspline = spline_points_array[ii_sp];

                double result, error, r;
                double calE;
                gsl_integration_workspace *w = gsl_integration_workspace_alloc(Nintegration);
                int i;

                gsl_function F;
                F.function = &massintegrand;
                F.params = NULL;

                double rmax = g_cored_profile_rmax_factor * g_cored_profile_rc;
                /** @note Calculate mass normalization factor for these params. */
                gsl_integration_qag(&F, 0.0, rmax, 0, 1.0e-12, Nintegration, 5, w, &result, &error);
                normalization = result;

                /** @note Calculate M(r) and create mass spline for these params. */
                int num_points = Nspline * 10; // Use more points for spline data generation than for integration
                double *mass = (double *)malloc(num_points * sizeof(double));
                double *radius = (double *)malloc(num_points * sizeof(double));

                for (int i = 0; i < num_points; i++)
                {
                    double r = (double)i * rmax / (num_points);
                    gsl_integration_qag(&F, 0.0, r, 0, 1.0e-12, Nintegration, 5, w, &result, &error);
                    mass[i] = result * g_cored_profile_halo_mass / normalization;
                    radius[i] = r;
                }

                gsl_interp_accel *enclosedmass = gsl_interp_accel_alloc();
                gsl_spline *splinemass = gsl_spline_alloc(gsl_interp_cspline, num_points);
                gsl_spline_init(splinemass, radius, mass, num_points);
                double rlow = radius[0];
                double rhigh = radius[num_points - 1];

                /** @note Write mass profile file for these params (e.g., data/massprofile_Ni1k_Ns1k.dat). */
                char fname[256];
                FILE *fp;
                char base_filename_massprofile[256];
                snprintf(base_filename_massprofile, sizeof(base_filename_massprofile), "data/massprofile_Ni%d_Ns%d.dat", Nintegration, Nspline);
                get_suffixed_filename(base_filename_massprofile, 1, fname, sizeof(fname));
                fp = fopen(fname, "wb"); // Binary mode for fprintf_bin
                for (r = 0.0; r < rhigh; r += rmax / 900.0)
                {
                    if (r < rlow || r > rhigh)
                    {
                        printf("r out of range\n");
                        CLEAN_EXIT(1);
                    }
                    fprintf_bin(fp, "%f %f\n", r, gsl_spline_eval(splinemass, r, enclosedmass));
                }
                fclose(fp);

                /** @note Calculate Psi(r) and create potential spline for these params. */
                double *Psivalues = (double *)malloc(num_points * sizeof(double));
                double *nPsivalues = (double *)malloc(num_points * sizeof(double)); // For inverse spline r(Psi)
                // Prepare Psiintegrand parameters for diagnostic loop
                Psiintegrand_params psi_params_diag;
                psi_params_diag.massintegrand_func = &massintegrand;
                psi_params_diag.params_for_massintegrand = NULL;

                gsl_function F_for_psi_diag;
                F_for_psi_diag.function = &Psiintegrand;
                F_for_psi_diag.params = &psi_params_diag;

                for (i = 0; i < num_points; i++)
                {
                    double r = (double)i * rmax / ((double)num_points);
                    double r1 = fmax(r, g_cored_profile_rc / 1000000.0);
                    gsl_integration_qagiu(&F_for_psi_diag, r1, 0, 1e-12, Nintegration, w, &result, &error);
                    double M_at_r1 = gsl_spline_eval(splinemass, r1, enclosedmass);
                    double first_term = G_CONST * M_at_r1 / r1;
                    double second_term = G_CONST * result * g_cored_profile_halo_mass / normalization;
                    Psivalues[i] = (first_term + second_term);
                    nPsivalues[i] = -Psivalues[i];
                }


                gsl_interp_accel *Psiinterp = gsl_interp_accel_alloc();
                gsl_spline *splinePsi = gsl_spline_alloc(gsl_interp_cspline, num_points);
                gsl_spline_init(splinePsi, radius, Psivalues, num_points);

                /** @note Write potential profile file for these params (e.g., data/Psiprofile_Ni1k_Ns1k.dat). */
                char base_filename_psiprofile[256];
                snprintf(base_filename_psiprofile, sizeof(base_filename_psiprofile), "data/Psiprofile_Ni%d_Ns%d.dat", Nintegration, Nspline);
                get_suffixed_filename(base_filename_psiprofile, 1, fname, sizeof(fname));
                fp = fopen(fname, "wb"); // Binary mode for fprintf_bin
                for (r = 0.0; r < ((double)num_points - 1.0) / ((double)num_points) * rmax; r += rmax / 900.0)
                {
                    if (r < rlow || r > rhigh)
                    {
                        printf("r out of range\n");
                        CLEAN_EXIT(1);
                    }
                    fprintf_bin(fp, "%f %f\n", r, evaluatespline(splinePsi, Psiinterp, r));
                }
                fclose(fp);

                /** @note Create inverse spline r(Psi) for these params. */
                gsl_interp_accel *rofPsiinterp = gsl_interp_accel_alloc();
                gsl_spline *splinerofPsi = gsl_spline_alloc(gsl_interp_cspline, num_points);
                gsl_spline_init(splinerofPsi, nPsivalues, radius, num_points);

                /** @note Calculate inner integral for f(E) and create f(E) spline for these params. */
                double *innerintegrandvalues = (double *)malloc((num_points + 1) * sizeof(double));
                double *Evalues = (double *)malloc((num_points + 1) * sizeof(double));
                double Psimin = Psivalues[num_points - 1];
                double Psimax = Psivalues[0];

                char base_filename_integrand[256];
                snprintf(base_filename_integrand, sizeof(base_filename_integrand), "data/integrand_Ni%d_Ns%d.dat", Nintegration, Nspline);
                get_suffixed_filename(base_filename_integrand, 1, fname, sizeof(fname));
                fp = fopen(fname, "wb"); // Binary mode for fprintf_bin
                calE = Psivalues[0];
                F.function = &fEintegrand;
                fEintegrand_params params = {calE, splinerofPsi, splinemass, rofPsiinterp, enclosedmass};
                F.params = &params;
                for (i = 0; i < num_points; i++)
                {
                    double t = sqrt(calE - Psimin) * ((double)i) / ((double)num_points);
                    fprintf_bin(fp, "%f %f\n", t, fEintegrand(t, &params));
                }
                fclose(fp);

                if (Psimax <= Psimin) { // Check for diagnostic loop
                    if (g_doDebug) log_message("DEBUG", "Diagnostic: Psimax (%.6e) <= Psimin (%.6e) in diagnostic loop", Psimax, Psimin);
                    // Continue with diagnostic but note the issue
                }

                innerintegrandvalues[0] = 0.0;
                Evalues[0] = Psimin;
                for (i = 1; i <= num_points; i++)
                {
                    calE = Psimin + (Psimax - Psimin) * ((double)i) / ((double)num_points);
                    if (i > 0 && calE <= Evalues[i-1]) { // Adjust if not strictly increasing
                        double prev_E_diag = Evalues[i-1];
                        double ideal_step_diag = (Psimax - Psimin) / (double)num_points;
                        double incr_diag = ideal_step_diag * 1e-6;
                        if(incr_diag == 0.0) incr_diag = DBL_MIN * fabs(prev_E_diag) + DBL_MIN;
                        if(incr_diag == 0.0) incr_diag = 1e-20;
                        calE = prev_E_diag + incr_diag;
                        if (g_doDebug && (i <= 3 || i > num_points - 3)) {
                            log_message("DEBUG", "Evalues[%d] adjusted to %.6e", i, calE);
                        }
                    }
                    Evalues[i] = calE;

                    fEintegrand_params params2 = {calE, splinerofPsi, splinemass, rofPsiinterp, enclosedmass};
                    F.params = &params2;
                    // Ensure sqrt argument is non-negative
                    double sqrt_arg_diag = calE - Psimin;
                    if (sqrt_arg_diag < 0) sqrt_arg_diag = 0.0;
                    gsl_integration_qag(&F, sqrt(sqrt_arg_diag) / 1.0e4, sqrt(sqrt_arg_diag), 1.0e-12, 1.0e-12, Nintegration, 6, w, &result, &error);
                    innerintegrandvalues[i] = result;
                }

                gsl_interp *fofEinterp = gsl_interp_alloc(gsl_interp_cspline, num_points + 1);
                gsl_interp_init(fofEinterp, Evalues, innerintegrandvalues, num_points + 1);
                gsl_interp_accel *fofEacc = gsl_interp_accel_alloc();

                /** @note Write theoretical density profile file for these params (e.g., data/density_profile_NiX_NsY.dat). */
                char base_filename[256];
                snprintf(base_filename, sizeof(base_filename), "data/density_profile_Ni%d_Ns%d.dat", Nintegration, Nspline);
                get_suffixed_filename(base_filename, 1, fname, sizeof(fname));
                fp = fopen(fname, "wb"); // Binary mode for fprintf_bin
                for (i = 0; i < num_points; i++)
                {
                    double rr = radius[i];
                    double rho_r = g_cored_profile_halo_mass / normalization * (1.0 / cube(1.0 + sqr(rr / g_cored_profile_rc)));
                    fprintf_bin(fp, "%f %f\n", rr, rho_r);
                }
                fclose(fp);

                /** @note Write dPsi/dr file (data/dpsi_dr<suffix>.dat) (overwrites previous if suffix same). */
                get_suffixed_filename("data/dpsi_dr.dat", 1, fname, sizeof(fname));
                fp = fopen(fname, "wb"); // Binary mode for fprintf_bin
                for (i = 0; i < num_points; i++)
                {
                    double rr = radius[i];
                    if (rr > 0.0)
                    {
                        double Menc = gsl_spline_eval(splinemass, rr, enclosedmass);
                        double dpsidr = -(G_CONST * Menc) / (rr * rr);
                        fprintf_bin(fp, "%f %f\n", rr, dpsidr);
                    }
                }
                fclose(fp);

                /** @note Write drho/dPsi file (data/drho_dpsi<suffix>.dat) (overwrites previous if suffix same). */
                get_suffixed_filename("data/drho_dpsi.dat", 1, fname, sizeof(fname));
                fp = fopen(fname, "wb"); // Binary mode for fprintf_bin
                for (i = 1; i < num_points - 1; i++)
                {
                    double rr = radius[i];
                    double rho_left = g_cored_profile_halo_mass / normalization * (1.0 / cube(1.0 + sqr(radius[i - 1] / g_cored_profile_rc)));
                    double rho_right = g_cored_profile_halo_mass / normalization * (1.0 / cube(1.0 + sqr(radius[i + 1] / g_cored_profile_rc)));
                    double drho_dr_num = (rho_right - rho_left) / (radius[i + 1] - radius[i - 1]);
                    double Menc = gsl_spline_eval(splinemass, rr, enclosedmass);
                    double dPsidr = -(G_CONST * Menc) / (rr * rr);

                    if (dPsidr != 0.0)
                    {
                        double Psi_val = evaluatespline(splinePsi, Psiinterp, rr);
                        fprintf_bin(fp, "%f %f\n", Psi_val, drho_dr_num / dPsidr);
                    }
                }
                fclose(fp);

                /** @note Write f(E) = dI/dE / const file for these params (e.g., data/f_of_E_NiX_NsY.dat). */
                char base_filename_fofe[256];
                snprintf(base_filename_fofe, sizeof(base_filename_fofe), "data/f_of_E_Ni%d_Ns%d.dat", Nintegration, Nspline);
                get_suffixed_filename(base_filename_fofe, 1, fname, sizeof(fname));
                fp = fopen(fname, "wb"); // Binary mode for fprintf_bin
                for (i = 0; i <= num_points; i++)
                {
                    double E = Evalues[i];
                    double deriv = 0.0;
                    if (i > 0 && i < num_points + 1)
                    {
                        if (i > 0 && i < num_points)
                        {
                            deriv = (innerintegrandvalues[i + 1] - innerintegrandvalues[i - 1]) / (Evalues[i + 1] - Evalues[i - 1]);
                        }
                        else if (i == 0)
                        {
                            deriv = (innerintegrandvalues[i + 1] - innerintegrandvalues[i]) / (Evalues[i + 1] - Evalues[i]);
                        }
                        else if (i == num_points)
                        {
                            deriv = (innerintegrandvalues[i] - innerintegrandvalues[i - 1]) / (Evalues[i] - Evalues[i - 1]);
                        }
                    }
                    double fE = fabs(deriv) / (sqrt(8.0) * PI * PI);
                    if (E == 0.0 || !isfinite(fE))
                        fE = 0.0;
                    fprintf_bin(fp, "%f %f\n", E, fE);
                }
                fclose(fp);

                get_suffixed_filename("data/df_fixed_radius.dat", 1, fname, sizeof(fname));
                fp = fopen(fname, "wb"); // Binary mode for fprintf_bin
                {
                    double r_fixed = 200.0;
                    double Psi_rf = evaluatespline(splinePsi, Psiinterp, r_fixed);
                    Psi_rf *= VEL_CONV_SQ;
                    int vsteps = 100;
                    for (int vv = 0; vv <= vsteps; vv++)
                    {
                        double vtest = (double)vv * (sqrt(2.0 * Psi_rf) / (vsteps));
                        double Etest = Psi_rf - 0.5 * vtest * vtest;
                        double fEval = 0.0;
                        if (Etest > Psimin && Etest < Psimax)
                        {
                            double dval = gsl_interp_eval_deriv(fofEinterp, Evalues, innerintegrandvalues, Etest, fofEacc);
                            fEval = dval / (sqrt(8.0) * PI * PI) * vtest * vtest * r_fixed * r_fixed;
                        }
                        fprintf_bin(fp, "%f %f\n", vtest, fEval);
                    }
                }
                fclose(fp);

                gsl_spline_free(splinePsi);
                gsl_spline_free(splinerofPsi);
                gsl_interp_accel_free(Psiinterp);
                gsl_interp_accel_free(rofPsiinterp);
                gsl_spline_free(splinemass);
                gsl_interp_accel_free(enclosedmass);
                gsl_interp_free(fofEinterp);
                gsl_interp_accel_free(fofEacc);
                free(mass);
                free(radius);
                if (radius_monotonic_grid_nfw != NULL) {
                    free(radius_monotonic_grid_nfw);
                    radius_monotonic_grid_nfw = NULL;
                }
                free(Psivalues);
                free(nPsivalues);
                free(innerintegrandvalues);
                free(Evalues);
                gsl_integration_workspace_free(w);
                w = NULL;
            }
        }
    }

    /**
     * @brief Main Theoretical Calculation (using default Nintegration=1000, Nspline=10000).
     * @details Repeats the theoretical profile calculations using fixed, default parameters
     *          (Nintegration=1000, Nspline=10000). The results from *this* block (splines:
     *          `splinemass`, `splinePsi`, `splinerofPsi`, `fofEinterp` and accelerators)
     *          are the ones used for generating the initial particle distribution and potentially
     *          for comparison during the simulation (e.g., debug energy calculation).
     *          Generates primary output files like `massprofile<suffix>.dat`, `Psiprofile<suffix>.dat`, `f_of_E<suffix>.dat`.
     */
    double r;

    /** @brief Allocate workspace for GSL integration operations. */
    w = gsl_integration_workspace_alloc(1000);

    gsl_function F;
    F.function = &massintegrand;
    F.params = NULL;

    rmax = g_cored_profile_rmax_factor * g_cored_profile_rc;
    gsl_integration_qag(&F, 0.0, rmax, 0, 1.0e-12, 1000, 5, w, &result, &error);
    normalization = result;

    num_points = 10000;
    /** @brief Allocate arrays for mass profile calculation. */
    mass = (double *)malloc(num_points * sizeof(double));
    radius = (double *)malloc(num_points * sizeof(double));

    for (int i = 0; i < num_points; i++)
    {
        double r = (double)i * rmax / (num_points);
        gsl_integration_qag(&F, 0.0, r, 0, 1.0e-12, 1000, 5, w, &result, &error);
        mass[i] = result * g_cored_profile_halo_mass / normalization;
        radius[i] = r;
    }

    /** @brief Create mass interpolation spline for M(r). */
    enclosedmass = gsl_interp_accel_alloc();
    splinemass = gsl_spline_alloc(gsl_interp_cspline, num_points);
    if (!check_strict_monotonicity(radius, num_points, "radius (main splinemass)")) {
        fprintf(stderr, "CRITICAL: radius array not monotonic for main splinemass\n");
        fflush(stderr);
        CLEAN_EXIT(1);
    }
    gsl_spline_init(splinemass, radius, mass, num_points);
    double rlow = radius[0];
    double rhigh = radius[num_points - 1];

    /** @brief Create 'data' directory again (harmless if exists). */
    {
        struct stat st = {0};
        if (stat("data", &st) == -1)
        {
            mkdir("data", 0755);
        }
    }

    /** @brief Write main mass profile file (data/massprofile<suffix>.dat). */
    get_suffixed_filename("data/massprofile.dat", 1, fname, sizeof(fname));
    fp = fopen(fname, "wb"); // Binary mode for fprintf_bin
    for (r = 0.0; r < rhigh; r += rmax / 900.0)
    {
        if (r < rlow || r > rhigh)
        {
            printf("r out of range\n");
            CLEAN_EXIT(1);
        }
        fprintf_bin(fp, "%f %f\n", r, gsl_spline_eval(splinemass, r, enclosedmass));
    }
    fclose(fp);

    /** @brief Calculate Psi(r) array (num_points=10000) and create primary `splinePsi`. */
    Psivalues = (double *)malloc(num_points * sizeof(double));
    nPsivalues = (double *)malloc(num_points * sizeof(double));
    // Prepare Psiintegrand parameters for Cored profile
    Psiintegrand_params psi_params_cored;
    psi_params_cored.massintegrand_func = &massintegrand;
    psi_params_cored.params_for_massintegrand = NULL;

    gsl_function F_for_psi_cored;
    F_for_psi_cored.function = &Psiintegrand;
    F_for_psi_cored.params = &psi_params_cored;

    for (i = 0; i < num_points; i++)
    {
        double r = (double)i * rmax / ((double)num_points);
        double r1 = fmax(r, g_cored_profile_rc / 1000000.0);
        gsl_integration_qagiu(&F_for_psi_cored, r1, 0, 1e-12, 1000, w, &result, &error);
        double first_term = G_CONST * gsl_spline_eval(splinemass, r1, enclosedmass) / r1;
        double second_term = G_CONST * result * g_cored_profile_halo_mass / normalization;
        Psivalues[i] = (first_term + second_term);
        nPsivalues[i] = -Psivalues[i];
    }

    /** @brief Create potential interpolation spline for Psi(r). */
    Psiinterp = gsl_interp_accel_alloc();
    splinePsi = gsl_spline_alloc(gsl_interp_cspline, num_points);
    if (!check_strict_monotonicity(radius, num_points, "radius (main splinePsi)")) {
        fprintf(stderr, "CRITICAL: radius array not monotonic for main splinePsi\n");
        fflush(stderr);
        CLEAN_EXIT(1);
    }
    gsl_spline_init(splinePsi, radius, Psivalues, num_points);


    /** @brief Write main potential profile file (data/Psiprofile<suffix>.dat). */
    get_suffixed_filename("data/Psiprofile.dat", 1, fname, sizeof(fname));
    fp = fopen(fname, "wb"); // Binary mode for fprintf_bin
    for (r = 0.0; r < ((double)num_points - 1.0) / ((double)num_points) * rmax; r += rmax / 900.0)
    {
        if (r < rlow || r > rhigh)
        {
            printf("r out of range\n");
            CLEAN_EXIT(1);
        }
        fprintf_bin(fp, "%f %f\n", r, evaluatespline(splinePsi, Psiinterp, r), evaluatespline(splinemass, enclosedmass, r));
    }
    fclose(fp);

    /** @brief Create inverse spline r(Psi) for radius lookup from potential. */
    rofPsiinterp = gsl_interp_accel_alloc();
    splinerofPsi = gsl_spline_alloc(gsl_interp_cspline, num_points);

    gsl_spline_init(splinerofPsi, nPsivalues, radius, num_points);

    /** @brief Calculate distribution function f(E) using Eddington's formula. */
    innerintegrandvalues = (double *)malloc((num_points + 1) * sizeof(double));
    Evalues = (double *)malloc((num_points + 1) * sizeof(double));
    Psimin = Psivalues[num_points - 1];
    Psimax = Psivalues[0];
    get_suffixed_filename("data/integrand.dat", 1, fname, sizeof(fname));
    fp = fopen(fname, "wb"); // Binary mode for fprintf_bin
    calE = Psivalues[0];
    F.function = &fEintegrand;
    fEintegrand_params params = {calE, splinerofPsi, splinemass, rofPsiinterp, enclosedmass};
    F.params = &params;
    for (i = 0; i < num_points; i++)
    {
        double t = sqrt(calE - Psimin) * ((double)i) / ((double)num_points);
        fprintf_bin(fp, "%f %f\n", t, fEintegrand(t, &params));
    }
    fclose(fp);


    innerintegrandvalues[0] = 0.0;
    Evalues[0] = Psimin;

    for (i = 1; i <= num_points; i++)
    {
        calE = Psimin + (Psimax - Psimin) * ((double)i) / ((double)num_points);
        fEintegrand_params params2 = {calE, splinerofPsi, splinemass, rofPsiinterp, enclosedmass};
        F.params = &params2;
        gsl_integration_qag(&F, sqrt(calE - Psimin) / 1.0e4, sqrt(calE - Psimin), 1.0e-12, 1.0e-12, 1000, 6, w, &result, &error);
        innerintegrandvalues[i] = result;
        Evalues[i] = calE;
    }

    /** @brief Create f(E) interpolation for particle generation. */

    g_main_fofEinterp = gsl_interp_alloc(gsl_interp_cspline, num_points + 1);
    gsl_interp_init(g_main_fofEinterp, Evalues, innerintegrandvalues, num_points + 1);
    g_main_fofEacc = gsl_interp_accel_alloc();

    /**
     * @brief Initialize particle arrays with default values (0.0, ID=index).
     * @details Uses the already allocated `particles` array from outer scope.
     *          Component indices: 0=radius, 1=velocity magnitude (initially), 2=ang. mom.,
     *          3=ID (initial index 0..npts_initial-1), 4=orientation(mu).
     */
    for (i = 0; i < npts_initial; i++)
    {
        particles[0][i] = 0.0;
        particles[1][i] = 0.0;
        particles[2][i] = 0.0;
        particles[3][i] = (double)i; // Initial ID is the index
        particles[4][i] = 0.0;
    }

    /**
     * Seed determination logic:
     * 1. If a specific seed (`--initial-cond-seed` or `--sidm-seed`) is provided, use it.
     * 2. Else if `--master-seed` is provided, derive specific seeds from it.
     * 3. Else if `--load-seeds` is specified (or by default if files exist and seeds not given), try to load from last_X_seed_{suffix}.dat.
     * 4. Else (no seeds provided, no load requested/possible), generate new seeds from time/pid.
     * Finally, save the seeds actually used to last_X_seed_{suffix}.dat and link last_X_seed.dat.
     */

    unsigned long int current_time_pid_seed = (unsigned long int)time(NULL) ^ (unsigned long int)getpid();
    char seed_filepath[512];
    FILE *fp_seed;

    // Determine Initial Conditions Seed
    if (!g_initial_cond_seed_provided) {
        if (g_master_seed_provided) {
            g_initial_cond_seed = g_master_seed + 1; // Deterministic offset
        } else if (g_attempt_load_seeds) {
            get_suffixed_filename(g_initial_cond_seed_filename_base, 1, seed_filepath, sizeof(seed_filepath));
            fp_seed = fopen(seed_filepath, "r");
            if (fp_seed) {
                if (fscanf(fp_seed, "%lu", &g_initial_cond_seed) == 1) {
                    log_message("INFO", "Loaded initial conditions seed %lu from %s", g_initial_cond_seed, seed_filepath);
                } else {
                    g_initial_cond_seed = current_time_pid_seed + 100; // Fallback if read fails
                    log_message("WARNING", "Failed to read IC seed from %s, generating new: %lu", seed_filepath, g_initial_cond_seed);
                }
                fclose(fp_seed);
            } else {
                g_initial_cond_seed = current_time_pid_seed + 100; // File not found, generate
                log_message("INFO", "No IC seed file found, generating new: %lu", g_initial_cond_seed);
            }
        } else {
            g_initial_cond_seed = current_time_pid_seed + 100; // Default generation
            log_message("INFO", "Generating new IC seed: %lu", g_initial_cond_seed);
        }
    } else {
        log_message("INFO", "Using user-provided IC seed: %lu", g_initial_cond_seed);
    }

    // Determine SIDM Seed
    if (!g_sidm_seed_provided) {
        if (g_master_seed_provided) {
            g_sidm_seed = g_master_seed + 2; // Deterministic offset, different from IC seed
        } else if (g_attempt_load_seeds) {
            get_suffixed_filename(g_sidm_seed_filename_base, 1, seed_filepath, sizeof(seed_filepath));
            fp_seed = fopen(seed_filepath, "r");
            if (fp_seed) {
                if (fscanf(fp_seed, "%lu", &g_sidm_seed) == 1) {
                    log_message("INFO", "Loaded SIDM seed %lu from %s", g_sidm_seed, seed_filepath);
                } else {
                    g_sidm_seed = current_time_pid_seed + 200; // Fallback
                    log_message("WARNING", "Failed to read SIDM seed from %s, generating new: %lu", seed_filepath, g_sidm_seed);
                }
                fclose(fp_seed);
            } else {
                g_sidm_seed = current_time_pid_seed + 200; // File not found, generate
                log_message("INFO", "No SIDM seed file found, generating new: %lu", g_sidm_seed);
            }
        } else {
            g_sidm_seed = current_time_pid_seed + 200; // Default generation
            log_message("INFO", "Generating new SIDM seed: %lu", g_sidm_seed);
        }
    } else {
        log_message("INFO", "Using user-provided SIDM seed: %lu", g_sidm_seed);
    }

    // If --readinit is used, we should generally use a specified/loaded SIDM seed
    // or a newly generated one, NOT one derived from IC seed, as ICs are fixed.
    if (doReadInit && !g_sidm_seed_provided && !g_master_seed_provided && !g_attempt_load_seeds) {
        // If reading ICs and no SIDM/master seed is given and not told to load, ensure SIDM seed is fresh.
        // This case might have already generated g_sidm_seed from current_time_pid_seed, which is fine.
        // If it was derived from g_initial_cond_seed (which itself might have been from time), it's also fine.
        // The logic above should correctly make g_sidm_seed independent of g_initial_cond_seed
        // if g_master_seed_provided is false.
    }

    // Save the seeds that will actually be used
    // Save Initial Conditions Seed
    get_suffixed_filename(g_initial_cond_seed_filename_base, 1, seed_filepath, sizeof(seed_filepath));
    fp_seed = fopen(seed_filepath, "w");
    if (fp_seed) {
        fprintf(fp_seed, "%lu\n", g_initial_cond_seed);
        fclose(fp_seed);
        log_message("INFO", "Saved initial conditions seed %lu to %s", g_initial_cond_seed, seed_filepath);
        // Create link/copy
        char linkname_ic[512];
        snprintf(linkname_ic, sizeof(linkname_ic), "%s.dat", g_initial_cond_seed_filename_base); // e.g. data/last_initial_seed.dat

        // Platform-dependent link/copy code (similar to lastparams.dat)
        #if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
            // Windows: use copy
            char copy_cmd[1024];
            sprintf(copy_cmd, "copy \"%s\" \"%s\"", seed_filepath, linkname_ic);
            if (system(copy_cmd) != 0) {
                log_message("WARNING", "Failed to copy IC seed file from %s to %s", seed_filepath, linkname_ic);
            }
        #else
            // Unix/Linux/macOS: use symbolic link
            unlink(linkname_ic); // Remove existing link if present
            // Extract basename for relative symlink within data directory
            const char *basename_ic = strrchr(seed_filepath, '/');
            basename_ic = basename_ic ? basename_ic + 1 : seed_filepath; // Skip the '/' or use full name if no '/'
            if (symlink(basename_ic, linkname_ic) != 0) {
                log_message("WARNING", "Failed to create symbolic link from %s to %s", basename_ic, linkname_ic);
            }
        #endif
    } else {
        log_message("ERROR", "Failed to save IC seed to %s", seed_filepath);
    }

    // Save SIDM Seed
    get_suffixed_filename(g_sidm_seed_filename_base, 1, seed_filepath, sizeof(seed_filepath));
    fp_seed = fopen(seed_filepath, "w");
    if (fp_seed) {
        fprintf(fp_seed, "%lu\n", g_sidm_seed);
        fclose(fp_seed);
        log_message("INFO", "Saved SIDM seed %lu to %s", g_sidm_seed, seed_filepath);
        // Create link/copy
        char linkname_sidm[512];
        snprintf(linkname_sidm, sizeof(linkname_sidm), "%s.dat", g_sidm_seed_filename_base); // e.g. data/last_sidm_seed.dat

        // Platform-dependent link/copy code
        #if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
            // Windows: use copy
            char copy_cmd[1024];
            sprintf(copy_cmd, "copy \"%s\" \"%s\"", seed_filepath, linkname_sidm);
            if (system(copy_cmd) != 0) {
                log_message("WARNING", "Failed to copy SIDM seed file from %s to %s", seed_filepath, linkname_sidm);
            }
        #else
            // Unix/Linux/macOS: use symbolic link
            unlink(linkname_sidm); // Remove existing link if present
            // Extract basename for relative symlink within data directory
            const char *basename_sidm = strrchr(seed_filepath, '/');
            basename_sidm = basename_sidm ? basename_sidm + 1 : seed_filepath; // Skip the '/' or use full name if no '/'
            if (symlink(basename_sidm, linkname_sidm) != 0) {
                log_message("WARNING", "Failed to create symbolic link from %s to %s", basename_sidm, linkname_sidm);
            }
        #endif
    } else {
        log_message("ERROR", "Failed to save SIDM seed to %s", seed_filepath);
    }

    /**
     * @brief INITIAL CONDITION HANDLING block.
     * @details Determines whether to generate new initial conditions or load existing ones,
     *          or skip if restarting.
     *          - If `doReadInit` is true: Loads from `readInitFilename` via `read_initial_conditions`.
     *          - If `g_doRestart` is true: Skips generation (assumes simulation data exists or will be checked).
     *          - Otherwise: Generates new particles using the Sample Generator.
     */

    if (doReadInit)
    {
        /** @brief Load existing initial conditions from specified file. */
        printf("Reading initial conditions from %s...\n", readInitFilename);
        read_initial_conditions(particles, npts_initial, readInitFilename);
    }
    else if (!g_doRestart) // Only generate if NOT reading init file and NOT restarting
    {
        /**
         * @brief SAMPLE GENERATOR block.
         * @details Generates the initial particle distribution (radius, velocity magnitude, orientation)
         *          for `npts_initial` particles, based on the theoretical equilibrium distribution
         *          function `f(E)` derived via Eddington's formula (using primary splines).
         *          Uses inverse transform sampling for radius (via M(r) spline) and rejection
         *          sampling for velocity magnitude (using `f(E)` derivative spline).
         */
        if (tidal_fraction > 0.0) log_message("INFO", "Cored IC Gen: Initial particle count before stripping: %d", npts_initial);
        log_message("INFO", "Cored IC Gen: Generating %d initial particle positions and velocities...", npts_initial);

        double vel, ratio, Psir, mu, maxv, maxvalue;

        /** @brief Set up GSL spline for radius as a function of enclosed mass: r(M).
         *         Used for inverse transform sampling of radius. */
        gsl_interp_accel *rofMaccel = gsl_interp_accel_alloc();
        gsl_spline *splinerofM = gsl_spline_alloc(gsl_interp_cspline, num_points);
        // Add checks for allocation failure
        if (!rofMaccel || !splinerofM) { /* Handle error */ CLEAN_EXIT(1); }
        gsl_spline_init(splinerofM, mass, radius, num_points);

        /** @brief Set up GSL spline for the maximum of `v^2 * f(E)` envelope at each radius `r`.
         *         Used for rejection sampling efficiency. f(E) proportional to dI/dE. */
        double *maxv2f = (double *)malloc(num_points * sizeof(double));
        if (!maxv2f) { /* Handle error */ CLEAN_EXIT(1); }

        for (int i_r = 1; i_r < num_points; i_r++)
        {
            maxvalue = 0.0;
            Psir = evaluatespline(splinePsi, Psiinterp, radius[i_r]);
            maxv = sqrt(2.0 * (Psir - Psimin));
            for (int j_v = 1; j_v < num_points - 2; j_v++)
            {
                vel = maxv * ((double)j_v) / ((double)num_points);
                double currentvalue = vel * vel * gsl_interp_eval_deriv(g_main_fofEinterp, Evalues, innerintegrandvalues, Psir - (0.5) * vel * vel, g_main_fofEacc);
                if (currentvalue > maxvalue) maxvalue = currentvalue;
            }
            maxv2f[i_r] = maxvalue;
        }
        // Extrapolate for r=0 assuming linear behavior near origin based on points 1 and 2
        if (num_points >= 3) {
            maxv2f[0] = 2.0 * maxv2f[1] - maxv2f[2];
            if (maxv2f[0] < 0.0) maxv2f[0] = 0.0; // Ensure non-negative
        } else if (num_points == 2) {
             maxv2f[0] = maxv2f[1]; // Simple fallback
        } else {
             maxv2f[0] = 0.0; // Fallback for very few points
        }

        gsl_interp_accel *maxv2faccel = gsl_interp_accel_alloc();
        gsl_spline *splinemaxv2f = gsl_spline_alloc(gsl_interp_cspline, num_points);
        // Add checks for allocation failure
        if (!maxv2faccel || !splinemaxv2f) { /* Handle error */ free(maxv2f); CLEAN_EXIT(1); }
        gsl_spline_init(splinemaxv2f, radius, maxv2f, num_points);

        /**
         * @brief Generate `npts_initial` particle samples using sampling methods.
         * @note Uses GSL random number generation. For thread-safety in parallel code,
         *       thread-specific RNG states should be used instead of the global state.
         */
        for (i = 0; i < npts_initial; i++) // Loop over particles to generate
        {
            /** @note 1. Choose radius `r` using inverse transform sampling on M(r). */
            double mass_frac_sample = gsl_rng_uniform(g_rng) * 0.999999; // Avoid sampling exactly M_total
            double mass_sample = mass_frac_sample * g_cored_profile_halo_mass;
            particles[0][i] = evaluatespline(splinerofM, rofMaccel, mass_sample);

            /** @note 2. Find velocity magnitude `v` using rejection sampling against `max(v^2*f(E))`. */
            maxvalue = evaluatespline(splinemaxv2f, maxv2faccel, particles[0][i]);
            Psir = evaluatespline(splinePsi, Psiinterp, particles[0][i]);
            maxv = sqrt(2.0 * (Psir - Psimin));
            int vflag = 0;
            while (vflag == 0)
            {
                vel = gsl_rng_uniform(g_rng) * maxv;
                // Evaluate target function (proportional to v^2 * f(E))
                double target_func_val = vel * vel * gsl_interp_eval_deriv(g_main_fofEinterp, Evalues, innerintegrandvalues, Psir - (0.5) * vel * vel, g_main_fofEacc);
                ratio = (maxvalue > 1e-15) ? (target_func_val / maxvalue) : 0.0; // Avoid division by zero
                if (gsl_rng_uniform(g_rng) < ratio)
                {
                    particles[1][i] = vel;
                    vflag = 1; // Accept
                }
            } // End rejection loop

            /** @note 3. Generate random velocity orientation `mu = cos(theta)`. */
            mu = (2.0 * gsl_rng_uniform(g_rng) - 1.0); // Uniform distribution in [-1, 1] for isotropy

            /** @note 4. Calculate angular momentum `L = r * v_tangential = r * v * sqrt(1-mu^2)`. */
            particles[2][i] = particles[0][i] * particles[1][i] * sqrt(1.0 - mu * mu);

            /** @note 5. Store orientation parameter `mu` (needed later for v_radial). */
            particles[4][i] = mu;

            /** @note 6. Store initial index as particle ID. */
            particles[3][i] = (double)i;
        } // End particle generation loop

        /** @brief Cleanup Sample Generator resources (splines, accelerators, temp arrays). */
        gsl_spline_free(splinerofM);
        gsl_interp_accel_free(rofMaccel);
        gsl_spline_free(splinemaxv2f);
        gsl_interp_accel_free(maxv2faccel);
        free(maxv2f);

    } // End Sample Generator block (if !doReadInit && !g_doRestart)

    } // End if/else for profile selection for IC generation

    /**
     * @brief Save generated initial conditions to file if requested (`--writeinit`).
     * @details Saves the initial state if `doWriteInit` is true and not in restart mode.
     *          The state saved is *before* tidal stripping and unit conversion, using `write_initial_conditions`.
     * @see write_initial_conditions
     */
    if (doWriteInit && !g_doRestart)
    {
        printf("Saving initial conditions to %s...\n", writeInitFilename);
        write_initial_conditions(particles, npts_initial, writeInitFilename);
    }

    /**
     * @brief TIDAL STRIPPING IMPLEMENTATION block.
     * @details If `tidal_fraction` > 0 and not in restart mode, simulates tidal stripping by removing the outermost
     *          fraction of particles based on radius. It first sorts the `npts_initial`
     *          particles by radius, then keeps only the innermost `npts` particles.
     *          It reallocates the `particles` array to the final size `npts` and remaps
     *          the original indices stored in `particles[3]` to ranks [0, npts-1] using
     *          `reassign_orig_ids_with_rank`.
     * @see reassign_orig_ids_with_rank
     * @see sort_particles_with_alg
     */
    if (!g_doRestart) // Skip stripping if restarting
    {
        /** @note Only show stripping message if `--ftidal` was used. */
        if (tidal_fraction > 0.0) printf("Tidal stripping: sorting and retaining inner %.1f%% of particles...\n", (1.0 - tidal_fraction) * 100.0);

        /** @note Sort all `npts_initial` particles by radius using basic quadsort. */
        sort_particles_with_alg(particles, npts_initial, "quadsort"); // Sorts by particles[0]

        /** @note Allocate new smaller arrays (`final_particles`) for the `npts` particles to keep. */
        double **final_particles = (double **)malloc(5 * sizeof(double *));
        if (final_particles == NULL)
        {
            fprintf(stderr, "Memory allocation failed for final_particles\n");
            CLEAN_EXIT(1);
        }

        /** @brief Copy innermost `npts` particles to final arrays and replace `particles` pointers. */
        for (int i = 0; i < 5; i++) // Loop over components
        {
            final_particles[i] = (double *)malloc(npts * sizeof(double));
            if (final_particles[i] == NULL)
            {
                fprintf(stderr, "Memory allocation failed for final_particles[%d]\n", i);
                CLEAN_EXIT(1);
            }
            /** @note Copy only the first `npts` elements (innermost after sort). */
            memcpy(final_particles[i], particles[i], npts * sizeof(double));

            /** @note Free original oversized array and update `particles[i]` pointer. */
            free(particles[i]);                // Free the original oversized array
            particles[i] = final_particles[i]; // particles[i] now points to the smaller array
        }
        free(final_particles); // Free the temporary ** structure, not the data arrays

        /** @note Only show completion message if `--ftidal` was used. */
        if (tidal_fraction > 0.0)
        {
            printf("Tidal stripping complete: %d particles retained.\n\n", npts);
        }

        /**
         * @brief Remap original IDs (now in `particles[3]` for the kept particles) to ranks [0, npts-1].
         * @details Ensures `particles[3][i]` holds the final rank ID (0 to npts-1)
         *          for the particle currently at index `i` after stripping and sorting.
         */
        reassign_orig_ids_with_rank(particles[3], npts);
    } // End tidal stripping block (!g_doRestart)

    /**
     * @brief VELOCITY UNIT CONVERSION and ORIENTATION block.
     * @details Converts particle velocity magnitude (`particles[1]`) and angular momentum
     *          (`particles[2]`) from simulation generation units (implicitly km/s from `f(E)`)
     *          to physical units used in timestepping (kpc/Myr). It also applies the
     *          orientation parameter `mu = v_radial / v_total` (stored in `particles[4]`)
     *          to `particles[1]` to get the actual radial velocity component for integration.
     *          The original velocity magnitude in `particles[1]` is overwritten.
     *          Skipped in restart mode (`g_doRestart`).
     * @see kmsec_to_kpcmyr
     */
    if (!g_doRestart)
    {
        for (i = 0; i < npts; i++) // Loop over final npts particles
        {
            // particles[1] holds velocity magnitude 'v' from Sample Generator
            // particles[4] holds orientation 'mu' from Sample Generator
            particles[1][i] *= particles[4][i]; ///< Apply orientation: v_rad = v * mu
            particles[1][i] *= kmsec_to_kpcmyr; ///< Convert v_rad [km/s] to [kpc/Myr]
            // particles[2] holds angular momentum L = r*v*sqrt(1-mu^2)
            particles[2][i] *= kmsec_to_kpcmyr; ///< Convert L [kpc*km/s] to [kpc^2/Myr]
            // particles[4] (mu) is no longer needed after this step.
        }
    }

    /**
     * @brief PARTICLE DATA EXPORT block (Initial State for Simulation).
     * @details Writes the initial state of all `npts` particles (after potential stripping,
     *          unit conversion, orientation application, and ID remapping) to a file
     *          named `data/particles<suffix>.dat`. This file represents the state at t=0
     *          entering the simulation loop.
     *          Format (binary via fprintf_bin):
     *          radius(kpc, float) v_radial(kpc/Myr, float) ang_mom(kpc^2/Myr, float) final_rank_id(float? check fprintf_bin)
     *          Skipped if `skip_file_writes` is true (restart mode).
     * @see fprintf_bin
     */
    if (!skip_file_writes)
    {
        char filename[256];
        get_suffixed_filename("data/particles.dat", 1, filename, sizeof(filename));
        FILE *fpp = fopen(filename, "wb"); // Binary mode for fprintf_bin output
        if (fpp == NULL)
        {
            fprintf(stderr, "Error opening file %s for writing initial particles\n", filename);
            exit(1); // Use CLEAN_EXIT?
        }

        /** @brief Write each particle's state (using fprintf_bin). */
        for (int i = 0; i < npts; i++) // i is current index (0..npts-1)
        {
            fprintf_bin(fpp, "%f %f %f %f\n", // Format string likely ignored by fprintf_bin beyond types
                        particles[0][i], // Radius (kpc)
                        particles[1][i], // Radial velocity (kpc/Myr)
                        particles[2][i], // Angular momentum (kpc²/Myr)
                        particles[3][i]);// Particle ID (final rank, written as float)
        }
        fclose(fpp);
        printf("Wrote initial particle state to %s\n", filename);
    } // End skip_file_writes block for particles.dat

    /**
     * @brief SIMULATION TIMESTEP CALCULATION block.
     * @details Calculates the characteristic dynamical time (`tdyn`) based on core radius (`RC`)
     *          and total mass (`HALO_MASS`). Uses this to determine the total simulation
     *          duration (`totaltime = tfinal_factor * tdyn`) and the individual timestep
     *          size (`dt = totaltime / Ntimes`) used in the integration loop.
     * @see tdyn
     * @see totaltime
     * @see dt
     */
    double characteristic_radius_for_tdyn;
    if (g_use_nfw_profile) {
        characteristic_radius_for_tdyn = g_nfw_profile_rc; // Set from g_scale_radius_param
    } else {
        characteristic_radius_for_tdyn = g_cored_profile_rc; // Set from g_scale_radius_param
    }
    double tdyn = 1.0 / sqrt((VEL_CONV_SQ * G_CONST) * g_active_halo_mass / cube(characteristic_radius_for_tdyn));
    double totaltime = (double)tfinal_factor * tdyn; ///< Total simulation time (Myr)
    double dt = totaltime / ((double)Ntimes);        ///< Individual timestep size (Myr)
    printf("Dynamical time tdyn = %.4f Myr\n", tdyn);
    printf("Total simulation time = %.4f Myr (%.1f tdyn)\n", totaltime, (double)tfinal_factor);
    printf("Timestep dt = %.6f Myr\n\n", dt);

    /** @brief Initialize simulation time tracking and progress reporting. */
    double time = 0.0;                   ///< Current simulation time (Myr)
    double start_time = omp_get_wtime(); ///< Wall-clock start time for timing
    /** @brief Setup progress reporting steps (array `print_steps` holding step numbers for 0%, 5%, ..., 100%). */
    int print_steps[21];
    for (int k = 0; k <= 20; k++) print_steps[k] = (int)floor(k * 0.05 * Ntimes); // Calculate steps for progress output

    /** @brief Flag to determine if simulation phase can be skipped. */
    int skip_simulation = 0;

    /** @brief Set up simulation tracking variables. */

    /**
     * @brief TRAJECTORY TRACKING SETUP block.
     * @details Allocates memory arrays (`trajectories`, `energies`, `velocities_arr`, etc.)
     *          for tracking the evolution of a small number (`num_traj_particles`, max 10)
     *          of selected particles over time. These particles are identified by their
     *          *final rank ID* (0 to `upper_npts_num_traj - 1`). Uses `ext_Ntimes` for array
     *          size to accommodate potential loop overruns or post-loop access.
     * @note The specific particles tracked are those with final rank IDs 0, 1, ..., 9 (or fewer if npts < 10).
     */
    /** @brief Determine number of particles to track (max 10 or npts). */
    int num_traj_particles = 10; ///< Max number of low-ID particles to track
    int upper_npts_num_traj = (num_traj_particles < npts) ? num_traj_particles : npts; ///< Actual number tracked

    /** @brief Allocate trajectory arrays [tracked_particle_index][time_step]. */
    // Index `p` corresponds to the particle with final_rank_id `p`.
    double **trajectories = (double **)malloc(upper_npts_num_traj * sizeof(double *)); // Radius (kpc)
    // Add checks for allocation failure
    for (i = 0; i < upper_npts_num_traj; i++) trajectories[i] = (double *)malloc(ext_Ntimes * sizeof(double)); // Add checks for allocation failure
    double **energies = (double **)malloc(upper_npts_num_traj * sizeof(double *)); // Relative Energy E_rel (per unit mass)
    // Add checks for allocation failure
    for (i = 0; i < upper_npts_num_traj; i++) energies[i] = (double *)malloc(ext_Ntimes * sizeof(double)); // Add checks for allocation failure

    /** @brief Allocate arrays for additional tracked properties [tracked_particle_index][time_step]. */
    double **velocities_arr = (double **)malloc(upper_npts_num_traj * sizeof(double *)); // Radial velocity (kpc/Myr)
    // Add checks for allocation failure
    double **mu_arr = (double **)malloc(upper_npts_num_traj * sizeof(double *));         // Radial direction cosine (v_rad / v_tot)
    // Add checks for allocation failure
    // Note: E_arr seems redundant with 'energies' array tracking E_rel. Verify necessity.
    double **E_arr = (double **)malloc(upper_npts_num_traj * sizeof(double *));          // Total relative energy E_rel
    // Add checks for allocation failure
    double **L_arr = (double **)malloc(upper_npts_num_traj * sizeof(double *));          // Angular momentum (kpc^2/Myr)
    // Add checks for allocation failure

    // Allocate memory for each particle's timestep history
    for (i = 0; i < upper_npts_num_traj; i++)
    {
        velocities_arr[i] = (double *)malloc(ext_Ntimes * sizeof(double)); // Add checks for allocation failure
        mu_arr[i] = (double *)malloc(ext_Ntimes * sizeof(double)); // Add checks for allocation failure
        E_arr[i] = (double *)malloc(ext_Ntimes * sizeof(double)); // If kept, needs allocation check
        L_arr[i] = (double *)malloc(ext_Ntimes * sizeof(double)); // Add checks for allocation failure
    }

    /**
     * @brief ENERGY AND ANGULAR MOMENTUM INITIALIZATION block.
     * @details Calculates the initial relative energy `E_rel = Psi - KE` and angular
     *          momentum `L` for *all* particles based on their initial state (after
     *          stripping/conversion/remapping). Stores these initial values in `E_i_arr` and `L_i_arr`,
     *          indexed by the particle's final rank ID (`particles[3]`). Uses the
     *          theoretical potential spline `splinePsi` for the potential energy term.
     *          These arrays `E_i_arr`, `L_i_arr` store the *initial* values for later comparison.
     * @see E_i_arr
     * @see L_i_arr
     */
    /** @brief Allocate arrays for initial E and L, indexed by final rank ID. */
    double *E_i_arr = (double *)malloc(npts * sizeof(double)); // Stores initial E_rel[final_rank_id]
    // Add checks for allocation failure
    double *L_i_arr = (double *)malloc(npts * sizeof(double)); // Stores initial L[final_rank_id]
    // Add checks for allocation failure

    /** @brief Calculate initial E_rel and L for each particle and store by final rank ID. */
    for (i = 0; i < npts; i++) // Loop through particles 0..npts-1 (current index 'i')
    {
        double rr = particles[0][i];            // Radius (at current index i)
        double vrad = particles[1][i];          // Radial velocity (at current index i)
        double ell = particles[2][i];           // Angular momentum (at current index i)
        int final_rank_id = (int)particles[3][i]; // Final rank ID (stored at current index i)

        // Ensure final_rank_id is within bounds [0, npts-1]
        if (final_rank_id < 0 || final_rank_id >= npts)
        {
            log_message("ERROR", "Invalid remapped_id %d encountered at index %d during initial E/L calculation.", final_rank_id, i);
            // Handle error appropriately, maybe skip or exit
            continue;
        }

        /** @brief Calculate theoretical potential Psi(r) using initial theoretical spline. */
        // Uses Psiinterp accelerator associated with theoretical splinePsi
        double Psi_val = evaluatespline(splinePsi, Psiinterp, rr);
        Psi_val *= VEL_CONV_SQ; // Convert to physical units (kpc/Myr)^2

        /** @brief Calculate relative energy E = Psi - (1/2)(v_r^2 + L^2/r^2). */
        double E_rel = Psi_val - 0.5 * (vrad * vrad + (ell * ell) / (rr * rr));

        /** @brief Store initial E_rel and L using the final rank ID as index. */
        E_i_arr[final_rank_id] = E_rel;
        L_i_arr[final_rank_id] = ell;
    } // End initial E/L calculation loop

    /**
     * @brief LOW ANGULAR MOMENTUM PARTICLE SELECTION block.
     * @details Identifies the `nlowest` particles with either the lowest initial absolute angular
     *          momentum (`use_closest_to_Lcompare` = 0) or initial angular momentum closest
     *          to a reference value `Lcompare` (`use_closest_to_Lcompare` = 1), based on
     *          the values stored in `L_i_arr` (indexed by final_rank_id).
     *          Stores the final rank IDs of these selected particles in the `chosen` array.
     *          Allocates tracking arrays (`lowestL_r`, `lowestL_E`, `lowestL_L`) for these particles.
     * @see LAndIndex
     * @see cmp_LAI
     * @see chosen
     * @see lowestL_r
     * @see lowestL_E
     * @see lowestL_L
     */
    {
        /** @brief Create temporary LAndIndex array (size npts) to facilitate sorting by L. */
        LAndIndex *LAI = (LAndIndex *)malloc(npts * sizeof(LAndIndex));
        if (!LAI) { fprintf(stderr, "Error: Failed to allocate LAI array\n"); CLEAN_EXIT(1); }
        for (int i = 0; i < npts; i++) // 'i' here is the final_rank_id
        {
            if (use_closest_to_Lcompare) // Mode 1: Closest to Lcompare
            {
                double L_initial = L_i_arr[i];
                LAI[i].L = (L_initial - Lcompare) * (L_initial - Lcompare); // Store squared difference
                // Store sign of L_initial, used later for reconstruction or potentially unused
                LAI[i].sign = (L_initial >= 0.0) ? 1 : -1;
            }
            else // Mode 0: Lowest absolute L (using signed L for now based on original code)
            {
                LAI[i].L = L_i_arr[i]; // Store initial L value (signed)
                LAI[i].sign = 0;       // Sign field not used for sorting in this mode
            }
            LAI[i].idx = i; // Store the final_rank_id associated with this L value
        }

        /** @brief Sort LAI array by L member (ascending L or ascending L_diff^2). */
        qsort(LAI, npts, sizeof(LAndIndex), cmp_LAI);

        /** @brief If using Mode 1 (closest to Lcompare), reconstruct actual L values in LAI[].L.
         *         Note: This is primarily for the debug log output below; the primary goal is selecting indices. */
        if (use_closest_to_Lcompare) {
            for (int i = 0; i < npts; i++) {
                // Reconstruct L = Lcompare +/- sqrt(diff^2). The stored sign determines +/-.
                // Assumes LAI[i].sign stored the sign of (L_initial - Lcompare) or similar intention.
                // If sign was sign(L_initial), this reconstruction isn't quite right but matches original apparent logic.
                LAI[i].L = Lcompare + LAI[i].sign * sqrt(LAI[i].L);
            }
        }

        /** @brief Log the lowest L particles (L value and final_rank_id) if not restarting. */
        if (!g_doRestart) {
             log_message("DEBUG", "Selecting %d lowest L particles (mode=%d):", nlowest, use_closest_to_Lcompare);
             for (int i = 0; i < nlowest; i++) {
                 log_message("DEBUG", "  Rank %d: L=%.6f ID=%d", i, LAI[i].L, LAI[i].idx);
             }
        }

        /** @brief Allocate `chosen` array to store the final_rank_ids of the selected particles. */
        chosen = (int *)malloc(nlowest * sizeof(int));
        if (!chosen) { fprintf(stderr, "Error: Failed to allocate chosen array\n"); CLEAN_EXIT(1); }

        /** @brief Store the final rank IDs of the nlowest L particles into `chosen`. */
        for (int i = 0; i < nlowest; i++) chosen[i] = LAI[i].idx;

        free(LAI);

        /** @brief Allocate tracking arrays for chosen low-L particles [chosen_index][time_step]. */
        // Index 'p' corresponds to the p-th particle in the 'chosen' array.
        lowestL_r = (double **)malloc(nlowest * sizeof(double *)); // Radius
        // Add checks for allocation failure
        lowestL_E = (double **)malloc(nlowest * sizeof(double *)); // Energy
        // Add checks for allocation failure
        lowestL_L = (double **)malloc(nlowest * sizeof(double *)); // Angular Momentum
        // Add checks for allocation failure

        // Allocate memory for each selected particle's time history
        for (int p = 0; p < nlowest; p++) {
            lowestL_r[p] = (double *)malloc(ext_Ntimes * sizeof(double)); // Add checks for allocation failure
            lowestL_E[p] = (double *)malloc(ext_Ntimes * sizeof(double)); // Add checks for allocation failure
            lowestL_L[p] = (double *)malloc(ext_Ntimes * sizeof(double)); // Add checks for allocation failure
        }
    } // End Low-L selection block

    /**
     * @brief SIMULATION DATA TRACKING SETUP block.
     * @details Allocates memory for tracking particle data during the simulation.
     *          Includes the `inverse_map` array (mapping final_rank_id to current array index after sorting),
     *          calculates `deltaM` (mass per particle), and allocates block storage arrays
     *          (`L_block`, `Rank_block`, `R_block`, `Vrad_block`) if `g_doAllParticleData` is enabled.
     *          These blocks store data for `block_size` timesteps before being written to disk.
     * @see inverse_map
     * @see deltaM
     * @see L_block
     * @see Rank_block
     * @see R_block
     * @see Vrad_block
     */
    /** @brief Allocate index map: `inverse_map[final_rank_id]` will store the current index `i` after sorting. */
    int *inverse_map = (int *)malloc(npts * sizeof(int));
    if (!inverse_map) { fprintf(stderr, "Error: Failed to allocate inverse_map\n"); CLEAN_EXIT(1); }

    /** @brief Allocate particle scatter state array for AB3 history management after SIDM scattering. */
    g_particle_scatter_state = (int *)calloc(npts, sizeof(int)); // Use calloc to initialize all to 0
    if (!g_particle_scatter_state) {
        fprintf(stderr, "Error: Failed to allocate g_particle_scatter_state array\n");
        CLEAN_EXIT(1);
    }

    /** @brief Calculate mass per particle based on *initial* particle count before stripping. Used for M(rank). */
    double deltaM = g_active_halo_mass / (double)npts_initial;

    /** @brief Define block storage size (number of timesteps per block written to all_particle_data.dat). */
    int block_size = 100;



    /** @brief Create filename for the main all-particle data output file. */
    char apd_filename[256]; // Filename for data/all_particle_data<suffix>.dat
    get_suffixed_filename("data/all_particle_data.dat", 1, apd_filename, sizeof(apd_filename));

    /** @brief Allocate primary block storage arrays (float/int for memory efficiency).
     *         These store data for `block_size` steps, indexed [step_in_block * npts + final_rank_id].
     *         Freed later via `cleanup_all_particle_data()`. */
    L_block = (float *)malloc((size_t)npts * block_size * sizeof(float));    // Angular momentum
    Rank_block = (int *)malloc((size_t)npts * block_size * sizeof(int));     // Particle rank (sorted index) at that step
    R_block = (float *)malloc((size_t)npts * block_size * sizeof(float));    // Radius
    Vrad_block = (float *)malloc((size_t)npts * block_size * sizeof(float)); // Radial velocity
    // Check allocation results
    if (!L_block || !Rank_block || !R_block || !Vrad_block) {
        fprintf(stderr, "Error: Failed to allocate block storage arrays.\n");
        free(L_block); free(Rank_block); free(R_block); free(Vrad_block); // Free any that were allocated
        L_block = NULL; Rank_block = NULL; R_block = NULL; Vrad_block = NULL; // Prevent double free in cleanup
        CLEAN_EXIT(1);
    }

    /**
     * @brief RESTART MODE DATA VERIFICATION block (`all_particle_data.dat` check).
     * @details If running in restart mode (`g_doRestart`), this block checks if the main particle
     *          data output file (`all_particle_data<suffix>.dat`) already exists and is non-empty.
     *          If it exists and contains data, the flag `skip_simulation` is set to 1. This flag
     *          will cause the main timestepping loop to be bypassed entirely, allowing the program
     *          to proceed directly to the post-simulation snapshot analysis phase.
     * @note This check happens *before* the main loop. A more detailed check using
     *       `find_last_processed_snapshot` occurs *later* if snapshot analysis needs restarting.
     * @see skip_simulation
     * @see find_last_processed_snapshot
     */
    if (g_doRestart)
    {
        /** @brief Check for existence and size of the particle data file. */
        FILE *check_file = fopen(apd_filename, "rb");
        if (check_file)
        {
            /** @brief Determine file size by seeking to end and getting position. */
            fseek(check_file, 0, SEEK_END);
            long file_size = ftell(check_file);
            fclose(check_file);

            if (file_size > 0)
            {
                /** @brief Valid data file found - enable simulation phase bypass. */
                char human_size[32];
                format_file_size(file_size, human_size, sizeof(human_size));
                printf("Restart mode: Found existing all_particle_data file '%s' (%s).\n",
                       apd_filename, human_size);
                skip_simulation = 1;

                /** @brief Display appropriate message based on particle count comparison. */
                if (npts == npts_initial) {
                    printf("Skipping simulation phase and proceeding to post-processing. "
                           "Initial condition generation already performed, skipping these steps.\n\n");
                } else {
                    printf("Skipping simulation phase and proceeding to post-processing. "
                           "Initial condition generation and tidal stripping already performed, "
                           "skipping these steps.\n\n");
                }
            }
            else {
                /** @brief File exists but contains no data - simulation required. */
                printf("Restart mode: Found empty all_particle_data file '%s'. Will run simulation.\n", apd_filename);
            }
        }
        else {
            /** @brief File not found - simulation required. */
            printf("Restart mode: all_particle_data file '%s' not found. Will run simulation.\n", apd_filename);
        }
    } // End restart check block

    /**
     * @brief OUTPUT FILE INITIALIZATION block (`all_particle_data.dat`).
     * @details Creates (or overwrites) an empty binary file `all_particle_data<suffix>.dat`
     *          if saving all particle data (`g_doAllParticleData` is true) AND the simulation
     *          is *not* being skipped (`skip_simulation` is false). This file will be appended to
     *          incrementally during the simulation timestepping loop via block writes.
     * @see append_all_particle_data_chunk_to_file
     * @see apd_filename
     * @see g_doAllParticleData
     * @see skip_simulation
     */
    if (g_doAllParticleData && !skip_simulation)
    {
        // Create/truncate the output file in binary write mode.
        FILE *fapd = fopen(apd_filename, "wb");
        if (!fapd) {
            fprintf(stderr, "Error: cannot create all_particle_data output file %s\n", apd_filename);
            CLEAN_EXIT(1);
        }
        fclose(fapd); // Close immediately, file is now ready for appending.
        printf("Initialized empty file for all particle data: %s\n", apd_filename);

        // Calculate and display expected file size
        long long expected_size = (long long)total_writes * (long long)npts * 16LL; // 16 bytes per particle record
        double size_gb = expected_size / (1024.0 * 1024.0 * 1024.0);
        double size_mb = expected_size / (1024.0 * 1024.0);
        double size_kb = expected_size / 1024.0;

        if (size_gb >= 1.0) {
            printf("All particle data file requires: %.1f GB (%lld bytes)\n", size_gb, expected_size);
        } else if (size_mb >= 1.0) {
            printf("All particle data file requires: %.1f MB (%lld bytes)\n", size_mb, expected_size);
        } else if (size_kb >= 1.0) {
            printf("All particle data file requires: %.1f KB (%lld bytes)\n", size_kb, expected_size);
        } else {
            printf("All particle data file requires: %lld bytes\n", expected_size);
        }

        // Calculate and display expected snapshot file sizes
        // Each snapshot has 2 files: unsorted (28 bytes/particle) and sorted (32 bytes/particle)
        long long snapshot_size = (long long)npts * (28LL + 32LL); // Total per snapshot pair
        long long total_snapshot_size = snapshot_size * (long long)noutsnaps;
        double snap_size_gb = total_snapshot_size / (1024.0 * 1024.0 * 1024.0);
        double snap_size_mb = total_snapshot_size / (1024.0 * 1024.0);
        double snap_size_kb = total_snapshot_size / 1024.0;

        printf("%d time snapshot files will require: ", noutsnaps);
        if (snap_size_gb >= 1.0) {
            printf("%.1f GB (%lld bytes)\n", snap_size_gb, total_snapshot_size);
        } else if (snap_size_mb >= 1.0) {
            printf("%.1f MB (%lld bytes)\n", snap_size_mb, total_snapshot_size);
        } else if (snap_size_kb >= 1.0) {
            printf("%.1f KB (%lld bytes)\n", snap_size_kb, total_snapshot_size);
        } else {
            printf("%lld bytes\n", total_snapshot_size);
        }

        // Calculate and display total disk space
        long long total_disk_space = expected_size + total_snapshot_size;
        double total_gb = total_disk_space / (1024.0 * 1024.0 * 1024.0);
        double total_mb = total_disk_space / (1024.0 * 1024.0);
        double total_kb = total_disk_space / 1024.0;

        printf("Total disk space required: ");
        if (total_gb >= 1.0) {
            printf("%.1f GB (%lld bytes)\n", total_gb, total_disk_space);
        } else if (total_mb >= 1.0) {
            printf("%.1f MB (%lld bytes)\n", total_mb, total_disk_space);
        } else if (total_kb >= 1.0) {
            printf("%.1f KB (%lld bytes)\n", total_kb, total_disk_space);
        } else {
            printf("%lld bytes\n", total_disk_space);
        }

        // Check available disk space
        long long available_space = get_available_disk_space("data/");
        if (available_space > 0) {
            double avail_gb = available_space / (1024.0 * 1024.0 * 1024.0);
            double avail_mb = available_space / (1024.0 * 1024.0);
            double avail_kb = available_space / 1024.0;

            printf("Available disk space: ");
            if (avail_gb >= 1.0) {
                printf("%.1f GB (%lld bytes)\n", avail_gb, available_space);
            } else if (avail_mb >= 1.0) {
                printf("%.1f MB (%lld bytes)\n", avail_mb, available_space);
            } else if (avail_kb >= 1.0) {
                printf("%.1f KB (%lld bytes)\n", avail_kb, available_space);
            } else {
                printf("%lld bytes\n", available_space);
            }

            // Check if we're within 5% of total available or insufficient
            double usage_after = (double)(available_space - total_disk_space) / (double)available_space;

            if (available_space < total_disk_space) {
                // Insufficient space
                fprintf(stderr, "\nError: Insufficient disk space!\n");
                fprintf(stderr, "Required: %.1f GB\n", total_gb);
                fprintf(stderr, "Available: %.1f GB\n", avail_gb);
                fprintf(stderr, "Shortfall: %.1f GB\n", total_gb - avail_gb);
                CLEAN_EXIT(1);
            } else if (usage_after < 0.05) {
                // Within 5% of capacity after simulation
                printf("\nWarning: Simulation will use %.1f%% of available disk space!\n",
                       (100.0 * total_disk_space / available_space));
                printf("After simulation: %.1f GB free (%.1f%% remaining)\n",
                       (available_space - total_disk_space) / (1024.0 * 1024.0 * 1024.0),
                       usage_after * 100.0);

                if (!prompt_yes_no("Continue")) {
                    printf("Aborting simulation.\n");
                    CLEAN_EXIT(0);
                }
            }
        } else {
            fprintf(stderr, "Warning: Could not determine available disk space.\n");
            if (!prompt_yes_no("Continue without disk space check")) {
                printf("Aborting simulation.\n");
                CLEAN_EXIT(0);
            }
        }

        printf("\n");

        /** @brief Display initial simulation progress. */
        printf("0%% complete, timestep 0/%d, time=0.0000 Myr, elapsed=0.00 s\n", Ntimes);
    }

    /**
     * @brief BLOCK STORAGE INITIALIZATION block (Initial L).
     * @details Copies initial angular momentum values (`L_i_arr`, indexed by final_rank_id)
     *          into the first timestep slot (index 0 implicitly) of the `L_block` storage array.
     *          Converts from double precision (`L_i_arr`) to single precision (`L_block`).
     *          This prepares the block storage for the first chunk of simulation data.
     * @note Assumes L_block is indexed [step_in_block * npts + final_rank_id].
     */
    for (i = 0; i < npts; i++) // Loop over final_rank_id 'i'
    {
        double l_val = L_i_arr[i];
        float lf = (float)l_val;
        // Store initial L in the slot for step 0 for this final_rank_id
        L_block[i] = lf; // Index 'i' corresponds to final_rank_id for the first block slot (step 0)
    }

    /**
     * @brief Store initial (t=0) approximate energy for the debug particle.
     * @details If `g_doDebug` is enabled, this block calculates the theoretical
     *          approximate energy `E = Psi - KE` for the particle with final rank ID
     *          `DEBUG_PARTICLE_ID` based on its *initial* state (before any
     *          timesteps) and stores it using `store_debug_approxE` at snapshot index 0.
     * @note Assumes `DEBUG_PARTICLE_ID` refers to the *final rank ID* after potential
     *       stripping and remapping. Retrieves the initial state using this ID as the index
     *       into the `particles` array *before* the first sort in the main loop.
     */
    if (g_doDebug)
    {
        int debug_id = DEBUG_PARTICLE_ID;
        if (debug_id >= 0 && debug_id < npts) {
            // Retrieve initial state using final_rank_id as index into initial particles array
            double r0 = particles[0][debug_id];
            double v0 = particles[1][debug_id]; // This is v_rad
            double l0 = particles[2][debug_id];

            /** @brief Evaluate theoretical Psi(r) at initial radius using initial spline. */
            double psi_val = evaluatespline(splinePsi, Psiinterp, r0) * VEL_CONV_SQ;

            /** @brief Calculate initial approximate energy E = Psi - KE. */
            double E_approx0 = psi_val - 0.5 * (v0 * v0 + (l0 * l0) / (r0 * r0));
            double time0 = 0.0;

            /** @brief Store initial energy in debug arrays at snapshot index 0. */
            store_debug_approxE(0, E_approx0, time0);
        } else {
             fprintf(stderr, "Warning: Invalid DEBUG_PARTICLE_ID %d (must be 0 <= ID < %d)\n",
                     debug_id, npts);
        }
    }

    static int nwrite_total = 0; ///< Counter for number of dtwrite-interval writes performed.

    /**
     * @brief Perform initial sort of particles by radius before time integration.
     * @details Uses `quadsort` for this initial sort. Executed only by the master thread
     *          via `#pragma omp single`. After this sort, the index `i` of `particles[:][i]`
     *          corresponds to its rank based on radius.
     */
    #pragma omp single
    {
        sort_particles_with_alg(particles, npts, "quadsort");
    }

    /**
     * @brief Main Simulation Timestepping Loop.
     * @details This loop iterates from j = 0 to Ntimes + dtwrite - 1.
     *          In each iteration `j` (representing timestep from \f$t_j\f$ to \f$t_{j+1}\f$):
     *          1. The appropriate N-body integration method (selected by `method_select`)
     *             is called to advance all particle positions and velocities over `dt` due
     *             to gravitational forces. This typically involves sorting particles by radius
     *             for rank-based force calculation and updating an inverse map from original
     *             particle ID to current sorted rank.
     *          2. If SIDM is enabled (`g_enable_sidm_scattering`), the `handle_sidm_step()`
     *             function is called to apply stochastic SIDM scattering to particle pairs,
     *             further modifying their velocities and angular momenta.
     *          3. Simulation time `time` is incremented by `dt`.
     *          4. Trajectory and energy data for selected particles (low-ID and low-L) are recorded.
     *          5. If it's a `dtwrite` interval, a block of full particle data (rank, R, Vrad, L)
     *             is stored and potentially appended to `all_particle_data.dat`. Debug energy
     *             for `DEBUG_PARTICLE_ID` is also calculated and stored if it's a snapshot step.
     *          6. Progress is printed to the console at 5% intervals.
     *          This entire loop is skipped if `skip_simulation` is true (e.g., in restart mode
     *          where only post-processing of existing `all_particle_data.dat` is required).
     */
    if (!skip_simulation)
    {
        for (int j = 0; j < Ntimes + dtwrite; j++) // Main time loop
        {
            int current_step;

            // Method_select = 1; Flag now.
            if (method_select == 0)
            {
/****************************/
// EULER METHOD
/****************************/
#pragma omp single
                {
                    sort_particles(particles, npts);
                }

#pragma omp barrier

#pragma omp parallel for default(shared) schedule(static)
                {
                    for (int idx = 0; idx < npts; idx++)
                    {
                        int orig_id = (int)particles[3][idx];
                        inverse_map[orig_id] = idx;
                    }
                }

#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    double r = particles[0][i];
                    double vrad = particles[1][i];
                    double ell = particles[2][i];
                    double drdt = vrad;

                    double force = gravitational_force(r, i, npts, G_CONST, g_active_halo_mass);
                    double dvdt = force + effective_angular_force(r, ell);

                    particles[0][i] += drdt * dt;
                    particles[1][i] += dvdt * dt;
                }

                // SIDM scattering: profile-aware scale radius selection and execution mode handling
                double current_active_rc_for_sidm = g_use_nfw_profile ? g_nfw_profile_rc : g_cored_profile_rc;
                handle_sidm_step(particles, npts, dt, time, current_active_rc_for_sidm, display_method, 0);

#pragma omp single
                {
                    current_step = j + 1;
                    time += dt;
                }
#pragma omp parallel for if (upper_npts_num_traj > 1000) schedule(static)
                for (int p = 0; p < upper_npts_num_traj; p++)
                {
                    int idx = inverse_map[p];
                    double rr = particles[0][idx];
                    double vrad = particles[1][idx];
                    double ell = particles[2][idx];
                    double Psi_val = evaluatespline(splinePsi, Psiinterp, rr);
                    Psi_val *= VEL_CONV_SQ;
                    double vtot = sqrt(vrad * vrad + (ell * ell) / (rr * rr));
                    double mu_val = vrad / vtot;
                    double E_rel = Psi_val - 0.5 * (vrad * vrad + (ell * ell) / (rr * rr));
                    double l_current = ell;
                    trajectories[p][j] = rr;
                    energies[p][j] = E_rel;
                    velocities_arr[p][j] = vrad;
                    mu_arr[p][j] = mu_val;
                    E_arr[p][j] = E_rel;
                    L_arr[p][j] = l_current;
                }
            }
            else if (method_select == 1)
            {
                /****************************/
                // LEAPFROG METHOD (POSITION HALF STEP)
                /****************************/

#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    particles[0][i] += particles[1][i] * (dt / 2.0);
                }

#pragma omp single
                {
                    sort_particles(particles, npts);
                }

#pragma omp barrier
#pragma omp parallel for default(shared) schedule(static)
                for (int idx = 0; idx < npts; idx++)
                {
                    int orig_id = (int)particles[3][idx];
                    inverse_map[orig_id] = idx;
                }
#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    double r = particles[0][i];
                    double ell = particles[2][i];

                    double force = gravitational_force(r, i, npts, G_CONST, g_active_halo_mass);
                    double dvdt = force + effective_angular_force(r, ell);

                    particles[1][i] += dvdt * dt;
                    particles[0][i] += particles[1][i] * (dt / 2.0);
                }

                // SIDM scattering after leapfrog drift completion
                double current_active_rc_for_sidm = g_use_nfw_profile ? g_nfw_profile_rc : g_cored_profile_rc;
                handle_sidm_step(particles, npts, dt, time, current_active_rc_for_sidm, display_method, 0);

#pragma omp single
                {
                    current_step = j + 1;
                    time += dt;
                }
#pragma omp parallel for if (upper_npts_num_traj > 1000) schedule(static)

                for (int p = 0; p < upper_npts_num_traj; p++)
                {
                    int idx = inverse_map[p];
                    double rr = particles[0][idx];
                    double vrad = particles[1][idx];
                    double ell = particles[2][idx];
                    double Psi_val = evaluatespline(splinePsi, Psiinterp, rr);
                    Psi_val *= VEL_CONV_SQ;
                    double vtot = sqrt(vrad * vrad + (ell * ell) / (rr * rr));
                    double mu_val = vrad / vtot;
                    double E_rel = Psi_val - 0.5 * (vrad * vrad + (ell * ell) / (rr * rr));
                    double l_current = ell;
                    trajectories[p][j] = rr;
                    energies[p][j] = E_rel;
                    velocities_arr[p][j] = vrad;
                    mu_arr[p][j] = mu_val;
                    E_arr[p][j] = E_rel;
                    L_arr[p][j] = l_current;
                }
            }
            else if (method_select == 2)
            {

#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    double r = particles[0][i];
                    double vrad = particles[1][i];
                    double ell = particles[2][i];

                    double force = gravitational_force(r, i, npts, G_CONST, g_active_halo_mass);
                    double dvdt = force + effective_angular_force(r, ell);

                    particles[1][i] = vrad + 0.5 * dvdt * dt;
                }

#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    double vrad = particles[1][i];
                    particles[0][i] += vrad * dt;
                }

#pragma omp single
                {
                    sort_particles(particles, npts);
                }
#pragma omp barrier

#pragma omp parallel for default(shared) schedule(static)
                for (int idx = 0; idx < npts; idx++)
                {
                    int orig_id = (int)particles[3][idx];
                    inverse_map[orig_id] = idx;
                }

#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    double r = particles[0][i];
                    double vrad = particles[1][i];
                    double ell = particles[2][i];

                    double force = gravitational_force(r, i, npts, G_CONST, g_active_halo_mass);
                    double dvdt = force + effective_angular_force(r, ell);

                    particles[1][i] = vrad + 0.5 * dvdt * dt;
                }

                // SIDM scattering after velocity half-step completion
                double current_active_rc_for_sidm = g_use_nfw_profile ? g_nfw_profile_rc : g_cored_profile_rc;
                handle_sidm_step(particles, npts, dt, time, current_active_rc_for_sidm, display_method, 0);

#pragma omp single
                {
                    current_step = j + 1;
                    time += dt;
                }
#pragma omp parallel for if (upper_npts_num_traj > 1000) schedule(static)

                for (int p = 0; p < upper_npts_num_traj; p++)
                {
                    int idx = inverse_map[p];
                    double rr = particles[0][idx];
                    double vrad = particles[1][idx];
                    double ell = particles[2][idx];
                    double Psi_val = evaluatespline(splinePsi, Psiinterp, rr);
                    Psi_val *= VEL_CONV_SQ;
                    double vtot = sqrt(vrad * vrad + (ell * ell) / (rr * rr));
                    double mu_val = vrad / vtot;
                    double E_rel = Psi_val - 0.5 * (vrad * vrad + (ell * ell) / (rr * rr));
                    double l_current = ell;
                    trajectories[p][j] = rr;
                    energies[p][j] = E_rel;
                    velocities_arr[p][j] = vrad;
                    mu_arr[p][j] = mu_val;
                    E_arr[p][j] = E_rel;
                    L_arr[p][j] = l_current;
                }
            }
            else if (method_select == 3)
            {
                /**
                 * @brief Full-step adaptive leapfrog integration.
                 * @details Performs single step integration from (r_n, v_n) to (r_{n+1}, v_{n+1})
                 *          using adaptive timestep control and the doAdaptiveFullLeap function.
                 */

                double velocity_tol = 1.0e-5;
                double radius_tol = 1.0e-5;
                // Int max_subdiv = 8383608;
                int max_subdiv = 1;
                int out_type = 0;

                /**
                 * @brief Radial sorting phase - organize particles by radius.
                 * @details Sorting improves cache locality and force calculation efficiency.
                 */
#pragma omp single
                {
                    sort_particles(particles, npts);
                }
                /**
                 * @brief Main integration loop - adaptive leapfrog update for each particle.
                 * @details Each particle is advanced independently using adaptive timestepping.
                 */
#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    double r = particles[0][i];
                    double v = particles[1][i];
                    double ell = particles[2][i];

                    // One full step => h = dt.
                    double r_new, v_new;
                    doAdaptiveFullLeap(
                        i, npts,
                        r, v,
                        ell,
                        dt,
                        radius_tol, velocity_tol,
                        max_subdiv,
                        G_CONST,
                        out_type,
                        &r_new, &v_new);

                    particles[0][i] = r_new;
                    particles[1][i] = v_new;
                }

                /**
                 * Particle sorting after integration (commented out)
                 *
                 * Optional re-sorting could be performed after each full step
                 * Currently disabled for performance reasons
                 */
                // #pragma omp single.
                // Sort_particles(particles, npts);.

#pragma omp parallel for default(shared) schedule(static)
                for (int idx = 0; idx < npts; idx++)
                {
                    int orig_id = (int)particles[3][idx];
                    inverse_map[orig_id] = idx;
                }
                // }

                // SIDM scattering after adaptive leapfrog completion
                double current_active_rc_for_sidm = g_use_nfw_profile ? g_nfw_profile_rc : g_cored_profile_rc;
                handle_sidm_step(particles, npts, dt, time, current_active_rc_for_sidm, display_method, 0);

                /**
                 * @brief Timestep wrap-up phase - update time and record particle states.
                 * @details Updates global time counter and records trajectory information
                 *          for selected particles.
                 */
#pragma omp single
                {
                    current_step = j + 1;
                    time += dt;
                }
#pragma omp parallel for if (upper_npts_num_traj > 1000) schedule(static)

                for (int p = 0; p < upper_npts_num_traj; p++)
                {
                    int idx = inverse_map[p];
                    double rr = particles[0][idx];
                    double vrad = particles[1][idx];
                    double ell = particles[2][idx];

                    double Psi_val = evaluatespline(splinePsi, Psiinterp, rr);
                    Psi_val *= VEL_CONV_SQ;

                    double vtot = sqrt(vrad * vrad + (ell * ell) / (rr * rr));
                    double mu_val = vrad / vtot;

                    double E_rel = Psi_val - 0.5 * (vrad * vrad + (ell * ell) / (rr * rr));

                    energies[p][j] = E_rel;
                    velocities_arr[p][j] = vrad;
                    mu_arr[p][j] = mu_val;
                    E_arr[p][j] = E_rel;
                    L_arr[p][j] = ell;
                }
            }
            else if (method_select == 4)
            {
                /**
                 * @brief Hybrid integration with adaptive method selection.
                 * @details Uses Levi-Civita regularization for close encounters (r < r_crit)
                 *          and standard leapfrog otherwise. Radius threshold r_crit is
                 *          dynamically calculated for each particle.
                 */

                double velocity_tol = 1.0e-8;
                double radius_tol = 1.0e-8;
                int max_subdiv = 4096 * 4096;
                int out_type = 2;
                int N_taumin = 1000;

                double alpha_param = 0.05;

#pragma omp single
                {
                    sort_particles(particles, npts);
                }
                /**
                 * @brief Main integration loop - method selection based on orbital parameters.
                 * @details Dynamically selects between standard leapfrog and Levi-Civita
                 *          regularization based on particle's radius and angular momentum.
                 */
#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    double r = particles[0][i];
                    double v = particles[1][i];
                    double ell = particles[2][i];

                    // Define critical radius r_crit for switching integration method:
                    // r_crit = alpha_param * (ell^2) / (G * M(r))
                    double r_crit = 0.0;
                    if (ell != 0.0)
                    {
                        double M_enc = ((double)i / (double)npts) * g_active_halo_mass;
                        double gravPart = (VEL_CONV_SQ * G_CONST) * M_enc;
                        r_crit = (ell * ell) * (alpha_param) / (gravPart);
                    }
                    else
                    {
                        r_crit = 0.0;
                    }

                    double r_new, v_new;

                    if ((r > 1.0e-30) && (r < r_crit))
                    {
                        doLeviCivitaLeapfrog(
                            i, npts,
                            r, v,
                            ell,
                            dt,
                            N_taumin,
                            G_CONST,
                            &r_new, &v_new);
                    }
                    else
                    {
                        doAdaptiveFullLeap(
                            i, npts,
                            r, v,
                            ell,
                            dt,
                            radius_tol, velocity_tol,
                            max_subdiv,
                            G_CONST,
                            out_type,
                            &r_new, &v_new);
                    }

                    particles[0][i] = r_new;
                    particles[1][i] = v_new;
                }

#pragma omp parallel for default(shared) schedule(static)
                for (int idx = 0; idx < npts; idx++)
                {
                    int orig_id = (int)particles[3][idx];
                    inverse_map[orig_id] = idx;
                }
                // }

                // SIDM scattering after hybrid integrator completion (Levi-Civita/adaptive leapfrog)
                double current_active_rc_for_sidm = g_use_nfw_profile ? g_nfw_profile_rc : g_cored_profile_rc;
                handle_sidm_step(particles, npts, dt, time, current_active_rc_for_sidm, display_method, 0);

#pragma omp single
                {
                    current_step = j + 1;
                    time += dt;
                }
#pragma omp parallel for if (upper_npts_num_traj > 1000) schedule(static)

                for (int p = 0; p < upper_npts_num_traj; p++)
                {
                    int idx = inverse_map[p];
                    double rr = particles[0][idx];
                    double vrad = particles[1][idx];
                    double ell = particles[2][idx];

                    double Psi_val = evaluatespline(splinePsi, Psiinterp, rr);
                    Psi_val *= VEL_CONV_SQ;

                    double vtot = sqrt(vrad * vrad + (ell * ell) / (rr * rr));
                    double mu_val = vrad / vtot;

                    double E_rel = Psi_val - 0.5 * (vrad * vrad + (ell * ell) / (rr * rr));

                    energies[p][j] = E_rel;
                    velocities_arr[p][j] = vrad;
                    mu_arr[p][j] = mu_val;
                    E_arr[p][j] = E_rel;
                    L_arr[p][j] = ell;
                }
            }
            else if (method_select == 5)
            {
                double velocity_tol = 1.0e-7;
                double radius_tol = 1.0e-7;
                int max_subdiv = 4096 * 4096 * 16;
                int out_type = 2;
                int N_taumin = 10;
                double alpha_param = 0.05;

#pragma omp single
                {
                    sort_particles(particles, npts);
                }

#pragma omp parallel for default(shared) schedule(static)
                for (int idx = 0; idx < npts; idx++)
                {
                    int orig_id = (int)particles[3][idx];
                    inverse_map[orig_id] = idx;
                }

                // SIDM scattering before adaptive orbital integration with Levi-Civita regularization
                double current_active_rc_for_sidm = g_use_nfw_profile ? g_nfw_profile_rc : g_cored_profile_rc;
                handle_sidm_step(particles, npts, dt, time, current_active_rc_for_sidm, display_method, 0);

#pragma omp parallel for default(shared) schedule(static)
                for (int i = 0; i < npts; i++)
                {
                    double r = particles[0][i];
                    double v = particles[1][i];
                    double ell = particles[2][i];

                    // Calculate r_crit, using special handling for particle i=0 (M_enc=0) to avoid division by zero.
                    double r_crit = 0.0;
                    if (fabs(ell) > 1.0e-30)
                    {
                        double M_enc;
                        if (i == 0)
                        {
                            // Special handling for i=0 to avoid zero mass in denominator
                            M_enc = 0.1 * (1.0 / (double)npts) * g_active_halo_mass;
                        }
                        else
                        {
                            M_enc = ((double)i / (double)npts) * g_active_halo_mass;
                        }
                        double gravPart = (VEL_CONV_SQ * G_CONST) * M_enc;
                        r_crit = (ell * ell) * alpha_param / gravPart;
                    }

                    double r_new, v_new;
                    if (r > 1.0e-30 && r < r_crit) // Switch based on critical radius
                    {
                        doAdaptiveFullLeviCivita(
                            i, npts, r, v, ell, dt, N_taumin,
                            radius_tol, velocity_tol, max_subdiv,
                            G_CONST, out_type, &r_new, &v_new);
                    }
                    else
                    {
                        doAdaptiveFullLeap(
                            i, npts, r, v, ell, dt,
                            radius_tol, velocity_tol, max_subdiv,
                            G_CONST, out_type, &r_new, &v_new);
                    }
                    particles[0][i] = r_new;
                    particles[1][i] = v_new;
                }

#pragma omp single
                {
                    current_step = j + 1;
                    time += dt;
                }
#pragma omp parallel for if (upper_npts_num_traj > 1000) schedule(static)
                for (int p = 0; p < upper_npts_num_traj; p++)
                {
                    int idx = inverse_map[p];
                    double rr = particles[0][idx];
                    double vr = particles[1][idx];
                    double ell = particles[2][idx];

                    double Psi_val = evaluatespline(splinePsi, Psiinterp, rr);
                    Psi_val *= VEL_CONV_SQ;

                    double vtot = sqrt(vr * vr + (ell * ell) / (rr * rr));
                    double mu_val = vr / vtot;
                    double E_rel = Psi_val - 0.5 * (vr * vr + (ell * ell) / (rr * rr));

                    trajectories[p][j] = rr;
                    energies[p][j] = E_rel;
                    velocities_arr[p][j] = vr;
                    mu_arr[p][j] = mu_val;
                    E_arr[p][j] = E_rel;
                    L_arr[p][j] = ell;
                }
            }
            else if (method_select == 6)
            {

                // Coefficients for 4th-order Forest-Ruth-Yoshida integrator (c1=c3).
                // Derived from: c1 = 1 / (2 - 2^(1/3)), c2 = 1 - 2*c1
                double c1 = 0.6756035959798289;
                double c2 = -0.3512071919596578; // = 1.0 - 2.0 * c1
                double c3 = c1;

                /**
                 * @brief STEP 1: Kick by (c1 * dt/2).
                 * @details Velocity update using the old position.
                 */
#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    double r = particles[0][i];
                    double vrad = particles[1][i];
                    double ell = particles[2][i];

                    double force = gravitational_force(r, i, npts, G_CONST, g_active_halo_mass);
                    double dvdt = force + effective_angular_force(r, ell);
                    particles[1][i] = vrad + 0.5 * c1 * dt * dvdt;
                }

                /**
                 * @brief STEP 2: Drift by (c1 * dt).
                 * @details Position update using the intermediate velocity v^*.
                 */
#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    double vrad = particles[1][i];
                    particles[0][i] += vrad * (c1 * dt);
                }

#pragma omp single
                {
                    sort_particles(particles, npts);
                }

#pragma omp barrier
#pragma omp parallel for default(shared) schedule(static)
                for (int idx = 0; idx < npts; idx++)
                {
                    int orig_id = (int)particles[3][idx];
                    inverse_map[orig_id] = idx;
                }
                // }

                /**
                 * @brief STEP 3: Kick by ((c1 + c2) * dt/2).
                 * @details Velocity update using the new position after the first drift.
                 */
#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    double r = particles[0][i];
                    double vrad = particles[1][i];
                    double ell = particles[2][i];

                    double force = gravitational_force(r, i, npts, G_CONST, g_active_halo_mass);
                    double dvdt = force + effective_angular_force(r, ell);
                    double coeff = 0.5 * (c1 + c2);
                    particles[1][i] = vrad + coeff * dt * dvdt;
                }

                /**
                 * @brief STEP 4: Drift by (c2 * dt).
                 * @details Position update using the intermediate velocity.
                 */
#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    double vrad = particles[1][i];
                    particles[0][i] += vrad * (c2 * dt);
                }

                /**
                 * @brief STEP 5: Kick by ((c2 + c3) * dt/2).
                 * @details Velocity update using the new position after the second drift.
                 */
#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    double r = particles[0][i];
                    double vrad = particles[1][i];
                    double ell = particles[2][i];

                    // Recompute acceleration at new position.
                    double force = gravitational_force(r, i, npts, G_CONST, g_active_halo_mass);
                    double dvdt = force + effective_angular_force(r, ell);

                    // Combination of the remaining half of c2 and half of c3.
                    double coeff = 0.5 * (c2 + c3);
                    particles[1][i] = vrad + coeff * dt * dvdt;
                }

                /**
                 * @brief STEP 6: Drift by (c3 * dt).
                 * @details Final position update in this integration step.
                 */
#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    double vrad = particles[1][i];
                    particles[0][i] += vrad * (c3 * dt);
                }

                /**
                 * @brief STEP 7: Kick by (c3 * dt/2).
                 * @details Final velocity update (half-kick) to complete the integration step.
                 */
#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    double r = particles[0][i];
                    double vrad = particles[1][i];
                    double ell = particles[2][i];

                    double force = gravitational_force(r, i, npts, G_CONST, g_active_halo_mass);
                    double dvdt = force + effective_angular_force(r, ell);
                    particles[1][i] = vrad + 0.5 * c3 * dt * dvdt;
                }

                // SIDM scattering after 4th-order Yoshida symplectic integration
                double current_active_rc_for_sidm = g_use_nfw_profile ? g_nfw_profile_rc : g_cored_profile_rc;
                handle_sidm_step(particles, npts, dt, time, current_active_rc_for_sidm, display_method, 0);

#pragma omp single
                {
                    current_step = j + 1;
                    time += dt;
                }
#pragma omp parallel for if (upper_npts_num_traj > 1000) schedule(static)

                for (int p = 0; p < upper_npts_num_traj; p++)
                {
                    int idx = inverse_map[p];
                    double rr = particles[0][idx];
                    double vrad = particles[1][idx];
                    double ell = particles[2][idx];
                    double Psi_val = evaluatespline(splinePsi, Psiinterp, rr);
                    Psi_val *= VEL_CONV_SQ;
                    double vtot = sqrt(vrad * vrad + (ell * ell) / (rr * rr));
                    double mu_val = vrad / vtot;
                    double E_rel = Psi_val - 0.5 * (vrad * vrad + (ell * ell) / (rr * rr));
                    double l_current = ell;
                    trajectories[p][j] = rr;
                    energies[p][j] = E_rel;
                    velocities_arr[p][j] = vrad;
                    mu_arr[p][j] = mu_val;
                    E_arr[p][j] = E_rel;
                    L_arr[p][j] = l_current;
                }
            }
            else if (method_select == 7)
            {
                /****************************/
                // RK4 METHOD
                // In similar style as
                // the Euler code above
                /****************************/

                double *r_orig_by_id = (double *)malloc(npts * sizeof(double));
                double *v_orig_by_id = (double *)malloc(npts * sizeof(double));
                double *k1r_by_id = (double *)malloc(npts * sizeof(double));
                double *k1v_by_id = (double *)malloc(npts * sizeof(double));
                double *k2r_by_id = (double *)malloc(npts * sizeof(double));
                double *k2v_by_id = (double *)malloc(npts * sizeof(double));
                double *k3r_by_id = (double *)malloc(npts * sizeof(double));
                double *k3v_by_id = (double *)malloc(npts * sizeof(double));
                double *k4r_by_id = (double *)malloc(npts * sizeof(double));
                double *k4v_by_id = (double *)malloc(npts * sizeof(double));
                double h = dt;

// Store original state by orig_id.
#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    int orig_id = (int)particles[3][i];
                    r_orig_by_id[orig_id] = particles[0][i];
                    v_orig_by_id[orig_id] = particles[1][i];
                }

// K1 calculation.
#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    int orig_id = (int)particles[3][i];
                    double r = particles[0][i];
                    double vrad = particles[1][i];
                    double ell = particles[2][i];

                    double drdt = vrad;
                    // Use the gravitational_force and effective_angular_force functions.
                    double force = gravitational_force(r, i, npts, G_CONST, g_active_halo_mass);
                    double dvdt = force + effective_angular_force(r, ell);

                    k1r_by_id[orig_id] = drdt;
                    k1v_by_id[orig_id] = dvdt;
                }

#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    int orig_id = (int)particles[3][i];
                    double r_mid = particles[0][i];
                    double v_mid = particles[1][i];
                    double ell = particles[2][i];

                    double drdt = v_mid;
                    // Use the gravitational_force and effective_angular_force functions.
                    double force = gravitational_force(r_mid, i, npts, G_CONST, g_active_halo_mass);
                    double dvdt = force + effective_angular_force(r_mid, ell);

                    k2r_by_id[orig_id] = drdt;
                    k2v_by_id[orig_id] = dvdt;
                }

#pragma omp single
                {
                    sort_particles(particles, npts);
                }

#pragma omp barrier
#pragma omp parallel for default(shared) schedule(static)
                for (int idx = 0; idx < npts; idx++)
                {
                    int orig_id = (int)particles[3][idx];
                    inverse_map[orig_id] = idx;
                }
                // }

#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    int orig_id = (int)particles[3][i];
                    double r_mid = particles[0][i];
                    double v_mid = particles[1][i];
                    double ell = particles[2][i];

                    double drdt = v_mid;
                    // Use the gravitational_force and effective_angular_force functions.
                    double force = gravitational_force(r_mid, i, npts, G_CONST, g_active_halo_mass);
                    double dvdt = force + effective_angular_force(r_mid, ell);

                    k3r_by_id[orig_id] = drdt;
                    k3v_by_id[orig_id] = dvdt;
                }

// K4 calculation.
#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    int orig_id = (int)particles[3][i];
                    double r_end = particles[0][i];
                    double v_end = particles[1][i];
                    double ell = particles[2][i];

                    double drdt = v_end;
                    // Use the gravitational_force and effective_angular_force functions.
                    double force = gravitational_force(r_end, i, npts, G_CONST, g_active_halo_mass);
                    double dvdt = force + effective_angular_force(r_end, ell);

                    k4r_by_id[orig_id] = drdt;
                    k4v_by_id[orig_id] = dvdt;
                }

// Final RK4 combination.
#pragma omp parallel for default(shared) schedule(static)
                for (i = 0; i < npts; i++)
                {
                    int orig_id = (int)particles[3][i];

                    double r_new = r_orig_by_id[orig_id] + (h / 6.0) * (k1r_by_id[orig_id] + 2.0 * k2r_by_id[orig_id] + 2.0 * k3r_by_id[orig_id] + k4r_by_id[orig_id]);
                    double v_new = v_orig_by_id[orig_id] + (h / 6.0) * (k1v_by_id[orig_id] + 2.0 * k2v_by_id[orig_id] + 2.0 * k3v_by_id[orig_id] + k4v_by_id[orig_id]);

                    particles[0][i] = r_new;
                    particles[1][i] = v_new;
                }

                // SIDM scattering after RK4 state update completion
                double current_active_rc_for_sidm = g_use_nfw_profile ? g_nfw_profile_rc : g_cored_profile_rc;
                handle_sidm_step(particles, npts, dt, time, current_active_rc_for_sidm, display_method, 0);

#pragma omp single
                {
                    current_step = j + 1;
                    time += h;
                }
#pragma omp parallel for if (upper_npts_num_traj > 1000) schedule(static)

                for (int p = 0; p < upper_npts_num_traj; p++)
                {
                    int idx = inverse_map[p];
                    double rr = particles[0][idx];
                    double vrad = particles[1][idx];
                    double ell = particles[2][idx];
                    double Psi_val = evaluatespline(splinePsi, Psiinterp, rr);
                    Psi_val *= VEL_CONV_SQ;
                    double vtot = sqrt(vrad * vrad + (ell * ell) / (rr * rr));
                    double mu_val = vrad / vtot;
                    double E_rel = Psi_val - 0.5 * (vrad * vrad + (ell * ell) / (rr * rr));
                    double l_current = ell;
                    trajectories[p][j] = rr;
                    energies[p][j] = E_rel;
                    velocities_arr[p][j] = vrad;
                    mu_arr[p][j] = mu_val;
                    E_arr[p][j] = E_rel;
                    L_arr[p][j] = l_current;
                }

                // Free RK4 arrays.
                free(r_orig_by_id);
                free(v_orig_by_id);
                free(k1r_by_id);
                free(k1v_by_id);
                free(k2r_by_id);
                free(k2v_by_id);
                free(k3r_by_id);
                free(k3v_by_id);
                free(k4r_by_id);
                free(k4v_by_id);
            }
            else if (method_select == 8)
            {
                // Static variables for Adams-Bashforth 3rd Order (AB3) method
                static int ab3_bootstrap_done = 0;           ///< Flag indicating if the AB3 bootstrap phase has been completed (0=no, 1=yes).
                static double **f_ab3_r = NULL;               ///< History array for \f$dr/dt\f$ derivatives. Indexed by `[history_slot (0..2)][original_particle_id]`. Slot 2 is most recent (\f$f_n\f$).
                static double **f_ab3_v = NULL;               ///< History array for \f$dv_{rad}/dt\f$ derivatives. Indexed by `[history_slot (0..2)][original_particle_id]`. Slot 2 is most recent (\f$f_n\f$).
                static double h_ab3_bootstrap_step;         ///< Timestep size (`dt`) used during the Euler steps of the bootstrap phase.

                // Adams-Bashforth coefficients for different orders
                static int ab3_num[3] = {23, -16, 5};        ///< Numerator coefficients for the AB3 formula: \f$y_{n+1} = y_n + (h/12) \sum (\text{ab3_num}_i \cdot f_{n-i})\f$.
                static int ab2_num[2] = {18, -6};           ///< Numerator coefficients for the AB2 formula (for comparison or fallback).
                static int ab3_den = 12;                     ///< Common denominator for the AB3 formula coefficients.

                static int bootstrap_euler_steps_needed = 2; ///< Number of full Euler steps (each of size \f$h_{bootstrap}\f$) required to generate the initial 3 derivative history points (\f$f_0, f_1, f_2\f$).

                // Allocate AB3 history arrays once.
                if (f_ab3_r == NULL)
                {
                    f_ab3_r = (double **)malloc(3 * sizeof(double *));
                    f_ab3_v = (double **)malloc(3 * sizeof(double *));
                    for (int hh = 0; hh < 3; hh++)
                    {
                        f_ab3_r[hh] = (double *)malloc(npts * sizeof(double));
                        f_ab3_v[hh] = (double *)malloc(npts * sizeof(double));
                    }
                    h_ab3_bootstrap_step = dt; // Step size equals dt.
                }

                // If bootstrap not done yet, do it ONCE to fill the 3-step derivative history.
                if (!ab3_bootstrap_done)
                {

                    // Bootstrap for AB3: Need to compute f0, f1, f2.
                    // This requires 2 Euler steps to get states y1, y2.
                    // sub_step = 0: calc f0 (from y0), store f_ab3_x[0]. Euler y0->y1. state is y1.
                    // sub_step = 1: calc f1 (from y1), store f_ab3_x[1]. Euler y1->y2. state is y2.
                    // sub_step = 2: calc f2 (from y2), store f_ab3_x[2]. NO Euler update. state is y2.
#pragma omp single
                    {
                        for (int sub_step = 0; sub_step <= bootstrap_euler_steps_needed; sub_step++)
                        {
                            sort_particles(particles, npts);
#pragma omp parallel for default(shared) schedule(static)
                            for (int idx = 0; idx < npts; idx++)
                            {
                                int orig_id = (int)particles[3][idx];
                                inverse_map[orig_id] = idx;
                            }
                            // Compute derivatives => store in f_ab3_r[sub_step], f_ab3_v[sub_step].
#pragma omp parallel for default(shared) schedule(static)
                            for (int i_eval = 0; i_eval < npts; i_eval++)
                            {
                                int orig_id = (int)particles[3][i_eval];
                                double rr = particles[0][i_eval];
                                double vrad = particles[1][i_eval];
                                double ell = particles[2][i_eval];

                                double drdt = vrad;
                                // Use the gravitational_force and effective_angular_force functions.
                                double force = gravitational_force(rr, i_eval, npts, G_CONST, g_active_halo_mass);
                                double dvdt = force + effective_angular_force(rr, ell);

                                f_ab3_r[sub_step][orig_id] = drdt;
                                f_ab3_v[sub_step][orig_id] = dvdt;
                            }

                            // If sub_step < bootstrap_euler_steps_needed (i.e., for sub_step 0 and 1),
                            // perform mini-substeps to advance particles to the next state with higher accuracy.
                            if (sub_step < bootstrap_euler_steps_needed)
                            {
                                double dt_mini = h_ab3_bootstrap_step / (double)NUM_MINI_SUBSTEPS_BOOTSTRAP;

#pragma omp parallel for default(shared) schedule(static)
                                for (int i_part = 0; i_part < npts; i_part++)
                                {
                                    // Each particle is evolved independently over NUM_MINI_SUBSTEPS_BOOTSTRAP
                                    double current_r_mini = particles[0][i_part];
                                    double current_vrad_mini = particles[1][i_part];
                                    double current_ell_mini = particles[2][i_part]; // Angular momentum (constant)

                                    // Perform NUM_MINI_SUBSTEPS_BOOTSTRAP mini-steps
                                    for (int m = 0; m < NUM_MINI_SUBSTEPS_BOOTSTRAP; m++)
                                    {
                                        // Calculate derivatives based on current mini-step state
                                        double drdt_m = current_vrad_mini;
                                        double force_m = gravitational_force(current_r_mini, i_part, npts, G_CONST, g_active_halo_mass);
                                        double dvdt_m = force_m + effective_angular_force(current_r_mini, current_ell_mini);

                                        // Euler update for this mini-step
                                        current_r_mini += dt_mini * drdt_m;
                                        current_vrad_mini += dt_mini * dvdt_m;
                                    }

                                    // After all mini-steps, update the main particles array
                                    particles[0][i_part] = current_r_mini;
                                    particles[1][i_part] = current_vrad_mini;
                                }
                            }
                        }
                    }

#pragma omp single
                    {
                        // Mark bootstrap done.
                        ab3_bootstrap_done = 1;
                    }
                }
                else
                {
                    /**
                     * Normal AB3 step each iteration
                     */

#pragma omp single
                    {
                        sort_particles(particles, npts);
                    }

#pragma omp barrier
// #pragma omp single.
// {
#pragma omp parallel for default(shared) schedule(static)
                    for (int idx = 0; idx < npts; idx++)
                    {
                        int orig_id = (int)particles[3][idx];
                        inverse_map[orig_id] = idx;
                    }
                    // }

#pragma omp parallel for default(shared) schedule(static)
                    for (int i = 0; i < npts; i++)
                    {
                        // Adams-Bashforth 8th order integration step.
                        int orig_id = (int)particles[3][i];
                        double rr = particles[0][i];
                        double vrad = particles[1][i];

                        double sum_r = 0.0;
                        double sum_v = 0.0;
                        int particle_state = g_particle_scatter_state[orig_id];

                        // AB3 coefficients: b0=23/12, b1=-16/12, b2=5/12. Denom ab3_den=12.
                        // History: f_ab3_[r/v][2] is f_n (latest), [1] is f_{n-1}, [0] is f_{n-2}

                        if (particle_state == 1) { // Just scattered: Use AB1 (Euler-like)
                            // sum = 12 * f_n
                            sum_r = 12.0 * f_ab3_r[2][orig_id];
                            sum_v = 12.0 * f_ab3_v[2][orig_id];
                            if (g_doDebug && i < 5) { // Extremely sparse debug
                                 log_message("DEBUG", "AB3_RESET: Particle %d (orig_id) using AB1 step (state 1)", orig_id);
                            }
                        } else if (particle_state == 2) { // One step after scatter: Use AB2
                            // sum = ab2_num[0] * f_n + ab2_num[1] * f_{n-1}
                            sum_r = ab2_num[0] * f_ab3_r[2][orig_id] + ab2_num[1] * f_ab3_r[1][orig_id];
                            sum_v = ab2_num[0] * f_ab3_v[2][orig_id] + ab2_num[1] * f_ab3_v[1][orig_id];
                            if (g_doDebug && i < 5) {
                                 log_message("DEBUG", "AB3_RESET: Particle %d (orig_id) using AB2 step (state 2)", orig_id);
                            }
                        } else { // Normal AB3 step
                            sum_r = ab3_num[0] * f_ab3_r[2][orig_id] + ab3_num[1] * f_ab3_r[1][orig_id] + ab3_num[2] * f_ab3_r[0][orig_id];
                            sum_v = ab3_num[0] * f_ab3_v[2][orig_id] + ab3_num[1] * f_ab3_v[1][orig_id] + ab3_num[2] * f_ab3_v[0][orig_id];
                        }

                        double r_next = rr + (dt / (double)ab3_den) * sum_r;
                        double v_next = vrad + (dt / (double)ab3_den) * sum_v;

                        particles[0][i] = r_next;
                        particles[1][i] = v_next;
                    }

                    // SIDM scattering after Adams-Bashforth update (skip during bootstrap)
                    double current_active_rc_for_sidm = g_use_nfw_profile ? g_nfw_profile_rc : g_cored_profile_rc;
                    handle_sidm_step(particles, npts, dt, time, current_active_rc_for_sidm, display_method, !ab3_bootstrap_done);

                    // We re-sort & compute new derivatives to shift the AB3 history.
#pragma omp single
                    {
                        sort_particles(particles, npts);
                    }
#pragma omp barrier
                    // #pragma omp single.
                    // {

#pragma omp parallel for default(shared) schedule(static)
                    for (int idx = 0; idx < npts; idx++)
                    {
                        int orig_id = (int)particles[3][idx];
                        inverse_map[orig_id] = idx;
                    }

                    // Recompute the derivatives for the new time => goes into f_ab3_r[2], f_ab3_v[2].
                    double **f_new_r = (double **)malloc(sizeof(double *));
                    double **f_new_v = (double **)malloc(sizeof(double *));
                    f_new_r[0] = (double *)malloc(npts * sizeof(double));
                    f_new_v[0] = (double *)malloc(npts * sizeof(double));

#pragma omp parallel for default(shared) schedule(static)
                    for (int i_dbg = 0; i_dbg < npts; i_dbg++)
                    {
                        int orig_id = (int)particles[3][i_dbg];
                        double rr = particles[0][i_dbg];
                        double vrad = particles[1][i_dbg];
                        double ell = particles[2][i_dbg];

                        double drdt = vrad;
                        // Use the gravitational_force and effective_angular_force functions.
                        double force = gravitational_force(rr, i_dbg, npts, G_CONST, g_active_halo_mass);
                        double dvdt = force + effective_angular_force(rr, ell);

                        f_new_r[0][orig_id] = drdt;
                        f_new_v[0][orig_id] = dvdt;
                    }

#pragma omp single
                    {
                        // SHIFT AB3 HISTORY: f0 <- f1, f1 <- f2
                        for (int i_s = 0; i_s < npts; i_s++)
                        {
                            f_ab3_r[0][i_s] = f_ab3_r[1][i_s]; // f_{n-2} becomes old f_{n-1}
                            f_ab3_v[0][i_s] = f_ab3_v[1][i_s];

                            f_ab3_r[1][i_s] = f_ab3_r[2][i_s]; // f_{n-1} becomes old f_n
                            f_ab3_v[1][i_s] = f_ab3_v[2][i_s];
                        }
                        // Put the new derivative (f_n for the just-completed step) in slot #2
                        for (int i_s = 0; i_s < npts; i_s++)
                        {
                            f_ab3_r[2][i_s] = f_new_r[0][i_s]; // f_n (latest)
                            f_ab3_v[2][i_s] = f_new_v[0][i_s];
                        }

                        free(f_new_r[0]);
                        free(f_new_v[0]);
                        free(f_new_r);
                        free(f_new_v);

                        current_step = j + 1;
                        time += dt;

                        // Advance particle scatter states for next AB step
                        if (ab3_bootstrap_done) { // Only advance state if AB is active and past bootstrap
                            for (int k_pstate = 0; k_pstate < npts; k_pstate++) {
                                // k_pstate here is the original_id since g_particle_scatter_state is indexed by orig_id
                                if (g_particle_scatter_state[k_pstate] == 2) {
                                    g_particle_scatter_state[k_pstate] = 0; // Transition from AB2 to full AB3
                                } else if (g_particle_scatter_state[k_pstate] == 1) {
                                    g_particle_scatter_state[k_pstate] = 2; // Transition from AB1 to AB2
                                }
                                // If state is 0, it remains 0 unless SIDM sets it to 1 in the next call to handle_sidm_step
                            }
                        }
                    }
                } // End of the "else" block for normal AB3.
            }
            // Record trajectory data for selected low-ID particles
#pragma omp parallel for if (upper_npts_num_traj > 1000) schedule(static)

            for (int p = 0; p < upper_npts_num_traj; p++)
            {
                int idx = inverse_map[p];
                double rr = particles[0][idx];
                double vrad = particles[1][idx];
                double ell = particles[2][idx];
                double Psi_val = evaluatespline(splinePsi, Psiinterp, rr) * VEL_CONV_SQ;
                double vtot = sqrt(vrad * vrad + (ell * ell) / (rr * rr));
                double mu_val = vrad / vtot;

                double E_rel = Psi_val - 0.5 * (vrad * vrad + (ell * ell) / (rr * rr));
                trajectories[p][j] = rr;
                energies[p][j] = E_rel;
                velocities_arr[p][j] = vrad;
                mu_arr[p][j] = mu_val;
                E_arr[p][j] = E_rel;
                L_arr[p][j] = ell;
            }

            // Record trajectory data for selected low-L particles
            int max_threads = omp_get_max_threads();
            gsl_interp_accel **thread_accel = (gsl_interp_accel **)malloc(max_threads * sizeof(gsl_interp_accel *));
            for (int i = 0; i < max_threads; i++)
            {
                thread_accel[i] = gsl_interp_accel_alloc();
            }

#pragma omp parallel for if (nlowest > 1000) schedule(static)
            for (int p = 0; p < nlowest; p++)
            {
                int thread_id = omp_get_thread_num();
                gsl_interp_accel *thread_safe_accel = thread_accel[thread_id];

                int idx_lowest = inverse_map[chosen[p]];

                double rr = particles[0][idx_lowest];
                double vrad = particles[1][idx_lowest];
                double ell = particles[2][idx_lowest];

                if (rr < 0.0)
                {
                    rr = 0.0000000001;
                }
                if (rr > rmax)
                {
                    rr = rmax;
                }

                // Use thread-safe accelerator instead of shared one.
                double Psi_val = evaluatespline(splinePsi, thread_safe_accel, rr) * VEL_CONV_SQ;
                double E_rel = Psi_val - 0.5 * (vrad * vrad + (ell * ell) / (rr * rr));

                lowestL_r[p][j] = rr;    // Store R.
                lowestL_E[p][j] = E_rel; // Store E.
                lowestL_L[p][j] = ell;   // Store L.
            }

            // Clean up thread-local accelerators.
            for (int i = 0; i < max_threads; i++)
            {
                gsl_interp_accel_free(thread_accel[i]);
            }
            free(thread_accel);
            // }

#pragma omp single
            {
                if ((current_step % dtwrite) == 0)
                {
                    double elapsed = omp_get_wtime() - start_time;
                    printf("Write data at timestep %d after %.2f s.\n",
                           current_step, elapsed);

                    // Calculate the write index (0-based) corresponding to this timestep
                    int nwrite = current_step / dtwrite - 1;
                    if (g_doAllParticleData)
                    {
                        // Calculate index within the current block (0 to block_size-1)
                        int block_index_apd = nwrite % block_size;
                        for (int pi = 0; pi < npts; pi++)
                        {
                            int orig_id = pi;
                            int rank = inverse_map[orig_id];
                            double par_r = particles[0][rank];
                            double par_vrad = particles[1][rank];
                            double par_ell = particles[2][rank];

                            Rank_block[block_index_apd * npts + orig_id] = rank;
                            R_block[block_index_apd * npts + orig_id] = (float)par_r;
                            Vrad_block[block_index_apd * npts + orig_id] = (float)par_vrad;
                            L_block[block_index_apd * npts + orig_id] = (float)par_ell;

                            if (g_doDebug)
                            {
                                // Check if this write corresponds to a desired snapshot output time
                                if (nwrite % stepBetweenSnaps == 0)
                                {
                                    int snapIndex = nwrite / stepBetweenSnaps;
                                    if (snapIndex < noutsnaps) // Ensure snapIndex is valid
                                    {
                                        int debug_id = DEBUG_PARTICLE_ID;

                                        // Retrieve current state for debug particle from block arrays
                                        float r_valF = R_block[block_index_apd * npts + debug_id];
                                        float v_valF = Vrad_block[block_index_apd * npts + debug_id];
                                        float l_valF = L_block[block_index_apd * npts + debug_id];
                                        double r_val = (double)r_valF;
                                        double v_val = (double)v_valF;
                                        double l_val = (double)l_valF;

                                        // Evaluate theoretical potential using original spline
                                        double psi_val = 0.0;
                                        if (r_val >= 0.0 && r_val <= rmax)
                                        {
                                            psi_val = evaluatespline(splinePsi, Psiinterp, r_val) * VEL_CONV_SQ;
                                        }

                                        // Calculate approximate energy E = Psi - KE
                                        double E_approx = psi_val - 0.5 * (v_val * v_val + (l_val * l_val) / (r_val * r_val));
                                        double sim_time = time;

                                        // Store the approximate energy for comparison
                                        store_debug_approxE(snapIndex, E_approx, sim_time);
                                    }
                                }
                            }
                        }
                    }
                    if (g_doAllParticleData)
                    {
                        // Append block to file if block is full
                        if (((nwrite + 1) % block_size) == 0 && nwrite > 0)
                        {
                            append_all_particle_data_chunk_to_file(apd_filename,
                                                                   npts,
                                                                   block_size,
                                                                   L_block,
                                                                   Rank_block,
                                                                   R_block,
                                                                   Vrad_block);
                            printf("Appended block ending write index %d to %s\n", nwrite, apd_filename);
                        }
                    }
                    // Increment count of dtwrite-based writes
                    nwrite_total++;
                }
            }

            for (int k = 0; k <= 20; k++)
            {
                if (current_step == print_steps[k])
                {
                    double elapsed = omp_get_wtime() - start_time;
                    int percent = k * 5;
                    printf("%d%% complete, timestep %d/%d, time=%f Myr, elapsed=%.2f s\n",
                           percent, current_step, Ntimes, time, elapsed);
                    break;
                }
            }
        }
        // After the main simulation loop, flush any remaining data in the last partial block.
        {
            if (g_doAllParticleData)
            {
                int leftover = nwrite_total % block_size;
                if (leftover > 0)
                {
                    printf("Flushing leftover %d steps from block storage...\n", leftover);
                    append_all_particle_data_chunk_to_file(apd_filename,
                                                           npts,
                                                           leftover, // Number of steps in this partial block
                                                           L_block,
                                                           Rank_block,
                                                           R_block,
                                                           Vrad_block);
                }
            }
        }
    } // Close the skip_simulation if block.

    // Write final particle state if simulation was run
    if (!skip_file_writes)
    {
        char suffixed_filename[256];
        get_suffixed_filename("data/particlesfinal.dat", 1, suffixed_filename, sizeof(suffixed_filename));
        fp = fopen(suffixed_filename, "wb"); // Binary mode for fprintf_bin
        if (fp == NULL)
        {
            fprintf(stderr, "Error opening file %s for writing\n", suffixed_filename);
            exit(1);
        }

        for (i = 0; i < npts; i++)
        {
            fprintf_bin(fp, "%f %f %f  %f\n", particles[0][i], particles[1][i], particles[2][i], particles[3][i]);
        }
        fclose(fp);
    }

    // Write theoretical profiles (profile-specific formulas)
    char suffixed_filename[256];

    if (g_use_nfw_profile) {
        /**
         * @brief Write final theoretical NFW profile characteristics to .dat files.
         * @details This block outputs several files (massprofile, Psiprofile, density_profile,
         *          dpsi_dr, drho_dpsi, f_of_E, df_fixed_radius) using the splines
         *          (e.g., splinemass, splinePsi, g_main_fofEinterp) and parameters
         *          (e.g., num_points, radius, normalization, g_nfw_profile_rc, etc.)
         *          that were established during the main NFW initial condition generation phase.
         *          Analytical formulas for NFW density and its derivatives are used where appropriate.
         */
        log_message("INFO", "Writing NFW theoretical profiles to final .dat files...");

        // Write NFW theoretical mass profile
        get_suffixed_filename("data/massprofile.dat", 1, suffixed_filename, sizeof(suffixed_filename));
        fp = fopen(suffixed_filename, "wb");
        if (fp) {
            for (double r_plot = 0.0; r_plot < radius[num_points - 1]; r_plot += (radius[num_points - 1] / 900.0)) {
                if (r_plot >= radius[0]) {
                    fprintf_bin(fp, "%f %f\n", r_plot, gsl_spline_eval(splinemass, r_plot, enclosedmass));
                }
            }
            if (num_points > 0) {
                 fprintf_bin(fp, "%f %f\n", radius[num_points-1], gsl_spline_eval(splinemass, radius[num_points-1], enclosedmass));
            }
            fclose(fp);
        } else {
            log_message("ERROR", "Failed to open %s for final NFW mass profile", suffixed_filename);
        }

        // Write NFW theoretical potential profile
        get_suffixed_filename("data/Psiprofile.dat", 1, suffixed_filename, sizeof(suffixed_filename));
        fp = fopen(suffixed_filename, "wb");
        if (fp) {
            for (double r_plot = 0.0; r_plot < radius[num_points - 1]; r_plot += (radius[num_points - 1] / 900.0)) {
                if (r_plot >= radius[0]) {
                     fprintf_bin(fp, "%f %f\n", r_plot, evaluatespline(splinePsi, Psiinterp, r_plot));
                }
            }
             if (num_points > 0) {
                 fprintf_bin(fp, "%f %f\n", radius[num_points-1], evaluatespline(splinePsi, Psiinterp, radius[num_points-1]));
             }
            fclose(fp);
        } else {
            log_message("ERROR", "Failed to open %s for final NFW Psi profile", suffixed_filename);
        }

        // Write NFW theoretical density profile
        get_suffixed_filename("data/density_profile.dat", 1, suffixed_filename, sizeof(suffixed_filename));
        fp = fopen(suffixed_filename, "wb");
        if (fp) {
            double nt_nfw_scaler_final = g_nfw_profile_halo_mass / (4.0 * M_PI * normalization);
            for (i = 0; i < num_points; i++) {
                double rr = radius[i];
                double rs_k = rr / g_nfw_profile_rc;
                double term_s_k = rs_k + 0.01;
                if (term_s_k <= 1e-9) term_s_k = 1e-9;
                double term_n_k = (1.0 + rs_k) * (1.0 + rs_k);
                double term_c_base_k = rs_k / g_nfw_profile_falloff_factor;
                double term_c_k = 1.0 + pow(term_c_base_k, 10.0);
                double rho_shape_k = (term_s_k < 1e-9 || term_n_k < 1e-9 || term_c_k < 1e-9) ? 0.0 : (1.0 / (term_s_k * term_n_k * term_c_k));
                if (rr < 1e-6 && term_s_k < 1e-3 && rho_shape_k == 0.0) {
                   rho_shape_k = 1.0 / (term_s_k * term_n_k * term_c_k);
                }
                double rho_r_k = nt_nfw_scaler_final * rho_shape_k;
                fprintf_bin(fp, "%f %f\n", rr, rho_r_k);
            }
            fclose(fp);
        } else {
            log_message("ERROR", "Failed to open %s for final NFW density profile", suffixed_filename);
        }

        // Write NFW theoretical dPsi/dr profile
        get_suffixed_filename("data/dpsi_dr.dat", 1, suffixed_filename, sizeof(suffixed_filename));
        fp = fopen(suffixed_filename, "wb");
        if (fp) {
            for (i = 0; i < num_points; i++) {
                double rr = radius[i];
                if (rr > 0.0) {
                    double Menc = gsl_spline_eval(splinemass, rr, enclosedmass);
                    double dpsidr = -(G_CONST * Menc) / (rr * rr);
                    fprintf_bin(fp, "%f %f\n", rr, dpsidr);
                }
            }
            fclose(fp);
        } else {
            log_message("ERROR", "Failed to open %s for final NFW dpsi/dr profile", suffixed_filename);
        }

        // Write NFW theoretical drho/dPsi profile
        get_suffixed_filename("data/drho_dpsi.dat", 1, suffixed_filename, sizeof(suffixed_filename));
        fp = fopen(suffixed_filename, "wb");
        if (fp) {
            double nt_nfw_scaler_final = g_nfw_profile_halo_mass / (4.0 * M_PI * normalization);
            for (i = 1; i < num_points - 1; i++) {
                double rr = radius[i];
                if (rr <= 1e-9) continue;

                // Use NFW derivative function
                double drhodr_val_k = drhodr_profile_nfwcutoff(rr, g_nfw_profile_rc, nt_nfw_scaler_final, g_nfw_profile_falloff_factor);

                double Menc_k = gsl_spline_eval(splinemass, rr, enclosedmass);
                double dPsidr_mag_k = (G_CONST * Menc_k) / (rr * rr);

                if (fabs(dPsidr_mag_k) > 1e-30) {
                    double Psi_val_k = evaluatespline(splinePsi, Psiinterp, rr);
                    double drho_dPsi_val_k = drhodr_val_k / dPsidr_mag_k;
                    fprintf_bin(fp, "%f %f\n", Psi_val_k, drho_dPsi_val_k);
                }
            }
            fclose(fp);
        } else {
            log_message("ERROR", "Failed to open %s for final NFW drho/dpsi profile", suffixed_filename);
        }

        // Write NFW theoretical f(E) profile
        get_suffixed_filename("data/f_of_E.dat", 1, suffixed_filename, sizeof(suffixed_filename));
        fp = fopen(suffixed_filename, "wb");
        if (fp) {
            for (i = 0; i <= num_points; i++) {
                double E = Evalues[i];
                double deriv = 0.0;
                if (i > 0 && i < num_points + 1) {
                    if (i > 0 && i < num_points) {
                        deriv = (innerintegrandvalues[i + 1] - innerintegrandvalues[i - 1]) / (Evalues[i + 1] - Evalues[i - 1]);
                    }
                    else if (i == 0) {
                        deriv = (innerintegrandvalues[i + 1] - innerintegrandvalues[i]) / (Evalues[i + 1] - Evalues[i]);
                    }
                    else if (i == num_points) {
                        deriv = (innerintegrandvalues[i] - innerintegrandvalues[i - 1]) / (Evalues[i] - Evalues[i - 1]);
                    }
                }
                double fE = fabs(deriv) / (sqrt(8.0) * PI * PI);
                if (E == 0.0 || !isfinite(fE))
                    fE = 0.0;
                fprintf_bin(fp, "%f %f\n", E, fE);
            }
            fclose(fp);
        } else {
            log_message("ERROR", "Failed to open %s for final NFW f(E) profile", suffixed_filename);
        }

        // Write NFW distribution function at a fixed radius if simulation was run
        if (!skip_file_writes) {
            get_suffixed_filename("data/df_fixed_radius.dat", 1, suffixed_filename, sizeof(suffixed_filename));
            fp = fopen(suffixed_filename, "wb");
            if (fp) {
                double r_F = 2.0 * g_nfw_profile_rc;  // r_F = 2 × NFW scale radius
                double Psi_rf = evaluatespline(splinePsi, Psiinterp, r_F);
                Psi_rf *= VEL_CONV_SQ;
                double Psimin_test = VEL_CONV_SQ * Psimin; // Convert to (km/s)² for velocity calculation

                int vsteps = 10000;
                int reduce_vsteps = 300;
                for (int vv = 0; vv <= vsteps - reduce_vsteps; vv++) {
                    double sqrt_arg_v = Psi_rf - Psimin_test;
                    if (sqrt_arg_v < 0) sqrt_arg_v = 0;
                    double vtest = (double)vv * (sqrt(2.0 * sqrt_arg_v) / (vsteps));
                    double Etest = Psi_rf - 0.5 * vtest * vtest;
                    Etest = Etest / VEL_CONV_SQ; // Convert back to code units for bounds check
                    double fEval = 0.0;
                    if (Etest >= Psimin && Etest <= Psimax) { // Bounds check in code units
                        double derivative;
                        int status = gsl_interp_eval_deriv_e(g_main_fofEinterp, Evalues, innerintegrandvalues, Etest, g_main_fofEacc, &derivative);
                        if (status == GSL_SUCCESS) {
                            fEval = derivative / (sqrt(8.0) * PI * PI) * vtest * vtest * r_F * r_F;
                        }
                    }
                    if (!isfinite(fEval)) fEval = 0.0;
                    fprintf_bin(fp, "%f %f\n", vtest, fEval);
                }
                fclose(fp);
            } else {
                log_message("ERROR", "Failed to open %s for final NFW df_fixed_radius", suffixed_filename);
            }
        }

    } else {
        /**
         * @brief Write final theoretical Cored Plummer-like profile characteristics to .dat files.
         * @details This block outputs several files (massprofile, Psiprofile, density_profile,
         *          dpsi_dr, drho_dpsi, f_of_E, df_fixed_radius) using the splines
         *          (e.g., splinemass, splinePsi, g_main_fofEinterp) and parameters
         *          (e.g., num_points, radius, normalization, g_cored_profile_rc, etc.)
         *          that were established during the main Cored Plummer initial condition generation phase.
         *          Analytical formulas for the Cored density and its derivatives are used where appropriate.
         */
        log_message("INFO", "Writing Cored theoretical profiles to final .dat files...");

        // Write Cored theoretical mass profile
        get_suffixed_filename("data/massprofile.dat", 1, suffixed_filename, sizeof(suffixed_filename));
        fp = fopen(suffixed_filename, "wb");
        if (fp) {
            for (double r_plot = 0.0; r_plot < radius[num_points - 1]; r_plot += (radius[num_points - 1] / 900.0)) {
                if (r_plot >= radius[0]) {
                    fprintf_bin(fp, "%f %f\n", r_plot, gsl_spline_eval(splinemass, r_plot, enclosedmass));
                }
            }
            if (num_points > 0) {
                 fprintf_bin(fp, "%f %f\n", radius[num_points-1], gsl_spline_eval(splinemass, radius[num_points-1], enclosedmass));
            }
            fclose(fp);
        } else {
            log_message("ERROR", "Failed to open %s for final cored mass profile", suffixed_filename);
        }

        // Write Cored theoretical potential profile
        get_suffixed_filename("data/Psiprofile.dat", 1, suffixed_filename, sizeof(suffixed_filename));
        fp = fopen(suffixed_filename, "wb");
        if (fp) {
            for (double r_plot = 0.0; r_plot < radius[num_points - 1]; r_plot += (radius[num_points - 1] / 900.0)) {
                if (r_plot >= radius[0]) {
                     fprintf_bin(fp, "%f %f\n", r_plot, evaluatespline(splinePsi, Psiinterp, r_plot));
                }
            }
             if (num_points > 0) {
                 fprintf_bin(fp, "%f %f\n", radius[num_points-1], evaluatespline(splinePsi, Psiinterp, radius[num_points-1]));
             }
            fclose(fp);
        } else {
            log_message("ERROR", "Failed to open %s for final cored Psi profile", suffixed_filename);
        }

        // Write Cored theoretical density profile
        get_suffixed_filename("data/density_profile.dat", 1, suffixed_filename, sizeof(suffixed_filename));
        fp = fopen(suffixed_filename, "wb");
        if (fp) {
            for (i = 0; i < num_points; i++) {
                double rr = radius[i];
                double rho_r = g_cored_profile_halo_mass / normalization * (1.0 / cube(1.0 + sqr(rr / g_cored_profile_rc)));
                fprintf_bin(fp, "%f %f\n", rr, rho_r);
            }
            fclose(fp);
        } else {
            log_message("ERROR", "Failed to open %s for final cored density profile", suffixed_filename);
        }

        // Write Cored theoretical dPsi/dr profile
        get_suffixed_filename("data/dpsi_dr.dat", 1, suffixed_filename, sizeof(suffixed_filename));
        fp = fopen(suffixed_filename, "wb");
        if (fp) {
            for (i = 0; i < num_points; i++) {
                double rr = radius[i];
                if (rr > 0.0) {
                    double Menc = gsl_spline_eval(splinemass, rr, enclosedmass);
                    double dpsidr = -(G_CONST * Menc) / (rr * rr);
                    fprintf_bin(fp, "%f %f\n", rr, dpsidr);
                }
            }
            fclose(fp);
        } else {
            log_message("ERROR", "Failed to open %s for final cored dpsi/dr profile", suffixed_filename);
        }

        // Write Cored theoretical drho/dPsi profile
        get_suffixed_filename("data/drho_dpsi.dat", 1, suffixed_filename, sizeof(suffixed_filename));
        fp = fopen(suffixed_filename, "wb");
        if (fp) {
            for (i = 1; i < num_points - 1; i++) {
                double rr = radius[i];
                double rho_left = g_cored_profile_halo_mass / normalization * (1.0 / cube(1.0 + sqr(radius[i - 1] / g_cored_profile_rc)));
                double rho_right = g_cored_profile_halo_mass / normalization * (1.0 / cube(1.0 + sqr(radius[i + 1] / g_cored_profile_rc)));
                double drho_dr_num = (rho_right - rho_left) / (radius[i + 1] - radius[i - 1]);
                double Menc = gsl_spline_eval(splinemass, rr, enclosedmass);
                double dPsidr = -(G_CONST * Menc) / (rr * rr);
                if (dPsidr != 0.0) {
                    double Psi_val = evaluatespline(splinePsi, Psiinterp, rr);
                    fprintf_bin(fp, "%f %f\n", Psi_val, drho_dr_num / dPsidr);
                }
            }
            fclose(fp);
        } else {
            log_message("ERROR", "Failed to open %s for final cored drho/dpsi profile", suffixed_filename);
        }

        // Write theoretical f(E) profile
        get_suffixed_filename("data/f_of_E.dat", 1, suffixed_filename, sizeof(suffixed_filename));
        fp = fopen(suffixed_filename, "wb");
        if (fp) {
            for (i = 0; i <= num_points; i++) {
                double E = Evalues[i];
                double deriv = 0.0;
                if (i > 0 && i < num_points + 1) {
                    if (i > 0 && i < num_points) {
                        deriv = (innerintegrandvalues[i + 1] - innerintegrandvalues[i - 1]) / (Evalues[i + 1] - Evalues[i - 1]);
                    }
                    else if (i == 0) {
                        deriv = (innerintegrandvalues[i + 1] - innerintegrandvalues[i]) / (Evalues[i + 1] - Evalues[i]);
                    }
                    else if (i == num_points) {
                        deriv = (innerintegrandvalues[i] - innerintegrandvalues[i - 1]) / (Evalues[i] - Evalues[i - 1]);
                    }
                }
                double fE = fabs(deriv) / (sqrt(8.0) * PI * PI);
                if (E == 0.0 || !isfinite(fE))
                    fE = 0.0;
                fprintf_bin(fp, "%f %f\n", E, fE);
            }
            fclose(fp);
        } else {
            log_message("ERROR", "Failed to open %s for final f(E) profile (%s)", suffixed_filename, g_use_nfw_profile ? "NFW" : "Cored");
        }

        // Write distribution function at a fixed radius if simulation was run
        if (!skip_file_writes) {
            get_suffixed_filename("data/df_fixed_radius.dat", 1, suffixed_filename, sizeof(suffixed_filename));
            fp = fopen(suffixed_filename, "wb");
            if (fp) {
                double r_F = 2.0 * g_cored_profile_rc;  // r_F = 2 × Cored scale radius
                double Psi_rf = evaluatespline(splinePsi, Psiinterp, r_F);
                Psi_rf *= VEL_CONV_SQ;
                double Psimin_test = VEL_CONV_SQ * Psimin; // Convert to (km/s)² for velocity calculation

                int vsteps = 10000;
                int reduce_vsteps = 300;
                for (int vv = 0; vv <= vsteps - reduce_vsteps; vv++) {
                    double sqrt_arg_v = Psi_rf - Psimin_test;
                    if (sqrt_arg_v < 0) sqrt_arg_v = 0;
                    double vtest = (double)vv * (sqrt(2.0 * sqrt_arg_v) / (vsteps));
                    double Etest = Psi_rf - 0.5 * vtest * vtest;
                    Etest = Etest / VEL_CONV_SQ; // Convert back to code units for bounds check
                    double fEval = 0.0;
                    if (Etest >= Psimin && Etest <= Psimax) { // Bounds check in code units
                        double derivative;
                        int status = gsl_interp_eval_deriv_e(g_main_fofEinterp, Evalues, innerintegrandvalues, Etest, g_main_fofEacc, &derivative);
                        if (status == GSL_SUCCESS) {
                            fEval = derivative / (sqrt(8.0) * PI * PI) * vtest * vtest * r_F * r_F;
                        }
                    }
                    if (!isfinite(fEval)) fEval = 0.0;
                    fprintf_bin(fp, "%f %f\n", vtest, fEval);
                }
                fclose(fp);
            } else {
                log_message("ERROR", "Failed to open %s for final df_fixed_radius (%s)", suffixed_filename, g_use_nfw_profile ? "NFW" : "Cored");
            }
        }
    }

    char filename[256];
    get_suffixed_filename("data/particles.dat", 1, filename, sizeof(filename));
    FILE *finit = fopen(filename, "rb"); // Binary mode for fscanf_bin
    if (!finit)
    {
        printf("Error: can't open %s\n", filename);
        CLEAN_EXIT(1);
    }
    double *r_initial = (double *)malloc(npts * sizeof(double));
    double *rv_initial = (double *)malloc(npts * sizeof(double));
    double *v_initial = (double *)malloc(npts * sizeof(double));
    double *l_initial = (double *)malloc(npts * sizeof(double));
    double rank_id;
    for (i = 0; i < npts; i++)
    {
        fscanf_bin(finit, "%f %f %f %f\n", &r_initial[i], &rv_initial[i], &l_initial[i], &rank_id);
    }
    fclose(finit);

    for (i = 0; i < npts; i++)
    {
        v_initial[i] = sqrt(rv_initial[i] * rv_initial[i] + l_initial[i] * l_initial[i] / (r_initial[i] * r_initial[i]));
    }

    double *r_final = (double *)malloc(npts * sizeof(double));
    double *v_final = (double *)malloc(npts * sizeof(double));
    for (i = 0; i < npts; i++)
    {
        r_final[i] = particles[0][i];
        v_final[i] = sqrt(particles[1][i] * particles[1][i] + particles[2][i] * particles[2][i] / (particles[0][i] * particles[0][i]));
    }

    /**
     * @brief Calculate percentile-based ranges for dynamic histogram binning.
     * @details Combines initial and final particle distributions to determine
     * appropriate histogram ranges using the 99th percentile with 1.2x padding.
     * This ensures virtually all particles are captured while avoiding outliers.
     */
    double *r_all_sorted = (double *)malloc(2 * npts * sizeof(double));
    double *v_all_sorted = (double *)malloc(2 * npts * sizeof(double));

    // Combine initial and final data for percentile calculation
    for (i = 0; i < npts; i++) {
        r_all_sorted[i] = r_initial[i];
        r_all_sorted[npts + i] = r_final[i];
        v_all_sorted[i] = fabs(v_initial[i]) * (1.0 / kmsec_to_kpcmyr); // Convert to km/s
        v_all_sorted[npts + i] = fabs(v_final[i]) * (1.0 / kmsec_to_kpcmyr); // Convert to km/s
    }

    // Sort arrays to find percentiles
    qsort(r_all_sorted, 2 * npts, sizeof(double), double_cmp);
    qsort(v_all_sorted, 2 * npts, sizeof(double), double_cmp);

    // Calculate 99th percentile with 1.2x multiplier for padding
    int p99_index = (int)(0.99 * (2 * npts - 1));
    double max_r_all = r_all_sorted[p99_index] * 1.2;
    double max_v_all = v_all_sorted[p99_index] * 1.2;

    // Ensure minimum ranges for very concentrated distributions
    if (max_r_all < 50.0) max_r_all = 50.0;   // Minimum 50 kpc
    if (max_v_all < 50.0) max_v_all = 50.0;   // Minimum 50 km/s

    // Log the calculated ranges
    if (g_enable_logging) {
        log_message("INFO", "2D Histogram dynamic ranges: r=[0, %.1f] kpc, v=[0, %.1f] km/s",
                    max_r_all, max_v_all);
    }

    free(r_all_sorted);
    free(v_all_sorted);

    /**
     * @brief Define histogram parameters with 400x400 bins for higher resolution.
     * @details Bin widths are calculated dynamically based on the percentile ranges
     * to ensure optimal coverage of the particle distribution.
     */
    #define HIST_NBINS 400

    double rbin_width = max_r_all / HIST_NBINS;
    double vbin_width = max_v_all / HIST_NBINS;

    double bin_width = max_r_all / HIST_NBINS;  // For 1D histogram

    /**
     * @brief Generate and write 2D phase-space histograms (radius vs. velocity magnitude).
     * @details Calculates 2D histograms of particle counts in (r, |v|) bins for both the
     *          initial (t=0, after IC generation and any stripping/conversion) and final
     *          (end of simulation) particle distributions. The bin ranges are dynamically
     *          determined based on the 99th percentile of the combined initial and final
     *          radial positions and velocity magnitudes, with a 1.2x padding factor.
     *          Outputs `2d_hist_initial.dat` and `2d_hist_final.dat`.
     */
    int hist_initial[HIST_NBINS][HIST_NBINS];
    memset(hist_initial, 0, sizeof(hist_initial));
    int hist_final[HIST_NBINS][HIST_NBINS];
    memset(hist_final, 0, sizeof(hist_final));

    for (i = 0; i < npts; i++)
    {
        int rbini = (int)(r_initial[i] / rbin_width);
        int vbini = (int)(fabs(v_initial[i]) * (1.0 / kmsec_to_kpcmyr) / vbin_width);
        if (rbini < HIST_NBINS && vbini < HIST_NBINS && rbini >= 0 && vbini >= 0)
            hist_initial[rbini][vbini]++;

        rbini = (int)(r_final[i] / rbin_width);
        vbini = (int)(fabs(v_final[i]) * (1.0 / kmsec_to_kpcmyr) / vbin_width);
        if (rbini < HIST_NBINS && vbini < HIST_NBINS && rbini >= 0 && vbini >= 0)
            hist_final[rbini][vbini]++;
    }

    // Write initial 2D histogram if simulation was run
    if (!skip_file_writes)
    {
        char suffixed_filename[256];
        get_suffixed_filename("data/2d_hist_initial.dat", 1, suffixed_filename, sizeof(suffixed_filename));
        fp = fopen(suffixed_filename, "wb"); // Binary mode for fprintf_bin
        for (int rr = 0; rr < HIST_NBINS; rr++)
        {
            for (int vv = 0; vv < HIST_NBINS; vv++)
            {
                fprintf_bin(fp, "%f %f %d\n", (rr + 0.5) * rbin_width, (vv + 0.5) * vbin_width, hist_initial[rr][vv]);
            }
            fprintf_bin(fp, "\n");
        }
        fclose(fp);
    } // End skip_file_writes block.

    // Write final 2D histogram if simulation was run
    if (!skip_file_writes)
    {
        get_suffixed_filename("data/2d_hist_final.dat", 1, suffixed_filename, sizeof(suffixed_filename));
        fp = fopen(suffixed_filename, "wb"); // Binary mode for fprintf_bin
        for (int rr = 0; rr < HIST_NBINS; rr++)
        {
            for (int vv = 0; vv < HIST_NBINS; vv++)
            {
                fprintf_bin(fp, "%f %f %d\n", (rr + 0.5) * rbin_width, (vv + 0.5) * vbin_width, hist_final[rr][vv]);
            }
            fprintf_bin(fp, "\n");
        }
        fclose(fp);
    } // End skip_file_writes block.

    /**
     * @brief Generate and write 1D radial distribution histograms.
     * @details Calculates 1D histograms of particle counts in radial bins for both the
     *          initial and final particle distributions. Uses the same dynamically determined
     *          radial binning as the 2D histograms.
     *          Outputs `combined_histogram.dat` with columns: r_bin_center, count_initial, count_final.
     */
    int hist_i[HIST_NBINS];
    memset(hist_i, 0, sizeof(hist_i));
    int hist_f[HIST_NBINS];
    memset(hist_f, 0, sizeof(hist_f));

    for (i = 0; i < npts; i++)
    {
        int b = (int)(r_initial[i] / bin_width);
        if (b < HIST_NBINS && b >= 0)
            hist_i[b]++;
        b = (int)(r_final[i] / bin_width);
        if (b < HIST_NBINS && b >= 0)
            hist_f[b]++;
    }

    // Write combined 1D radius histogram if simulation was run
    if (!skip_file_writes)
    {
        char filename[256];
        get_suffixed_filename("data/combined_histogram.dat", 1, filename, sizeof(filename));
        fp = fopen(filename, "wb"); // Binary mode for fprintf_bin
        for (i = 0; i < HIST_NBINS; i++)
        {
            double bin_center = (i + 0.5) * bin_width;
            fprintf_bin(fp, "%f %d %d\n", bin_center, hist_i[i], hist_f[i]);
        }
        fclose(fp);
    } // End skip_file_writes block.

    /**
     * @brief Write trajectory data for a selection of low-original-ID particles.
     * @details Outputs the time evolution of radius, radial velocity, and the radial
     *          direction cosine (mu = v_rad / v_total) for the first `num_traj_particles`
     *          (typically 10, or fewer if npts < 10) particles, identified by their
     *          final rank ID after any stripping and remapping.
     *          The data is written to `trajectories.dat`.
     */
    if (!skip_file_writes)
    {
        get_suffixed_filename("data/trajectories.dat", 1, suffixed_filename, sizeof(suffixed_filename));
        fp = fopen(suffixed_filename, "wb"); // Binary mode for fprintf_bin
        for (int step = 0; step < Ntimes; step++)
        {
            fprintf_bin(fp, "%f", step * dt);
            for (int p = 0; p < num_traj_particles; p++)
            {
                fprintf_bin(fp, " %f %f %f", trajectories[p][step], velocities_arr[p][step], mu_arr[p][step]);
            }
            fprintf_bin(fp, "\n");
        }
        fclose(fp);
    } // End skip_file_writes block.

    /**
     * @brief Write trajectory data for the particle with final rank ID 0.
     * @details Outputs the time evolution of radius, radial velocity, and mu for the
     *          particle that ended up with rank ID 0. This is a subset of the data
     *          in `trajectories.dat`. Written to `single_trajectory.dat`.
     */
    if (!skip_file_writes)
    {
        get_suffixed_filename("data/single_trajectory.dat", 1, suffixed_filename, sizeof(suffixed_filename));
        fp = fopen(suffixed_filename, "wb"); // Binary mode for fprintf_bin
        for (int step = 0; step < Ntimes; step++)
        {
            fprintf_bin(fp, "%f %f %f %f\n", step * dt, trajectories[0][step], velocities_arr[0][step], mu_arr[0][step]);
        }
        fclose(fp);
    } // End skip_file_writes block.

    /**
     * @brief Write energy and angular momentum evolution for low-original-ID particles.
     * @details Outputs the time evolution of the current relative energy (E_cur), initial
     *          relative energy (E_i), current angular momentum (L_cur), and initial angular
     *          momentum (L_i) for the same set of `num_traj_particles` tracked for `trajectories.dat`.
     *          Written to `energy_and_angular_momentum_vs_time.dat`.
     */
    // Write energy/angular momentum evolution for low-ID particles if simulation was run
    if (!skip_file_writes)
    {
        get_suffixed_filename("data/energy_and_angular_momentum_vs_time.dat", 1, suffixed_filename, sizeof(suffixed_filename));
        fp = fopen(suffixed_filename, "wb"); // Binary mode for fprintf_bin
        for (int step = 0; step < Ntimes; step++)
        {
            fprintf_bin(fp, "%f", step * dt);
            for (int p = 0; p < num_traj_particles; p++)
            {
                double E_i = E_i_arr[p];
                double l_i = L_i_arr[p];
                double Ecur = E_arr[p][step];
                double lcur = L_arr[p][step];
                fprintf_bin(fp, " %f %f %f %f", Ecur, E_i, lcur, l_i);
            }
            fprintf_bin(fp, "\n");
        }
        fclose(fp);
    } // End skip_file_writes block.

    {
        /**
         * @brief Write trajectory data (radius, energy, angular momentum) for selected low-L particles.
         * @details Outputs the time evolution of radius, relative energy, and angular momentum
         *          for `nlowest` particles selected based on their initial angular momentum
         *          (either lowest absolute L or closest to a reference L, per `use_closest_to_Lcompare`).
         *          The specific particles are stored in the `chosen` array (by their original IDs).
         *          Written to `lowest_l_trajectories.dat`.
         */
        // Write trajectories for selected lowest-L particles if simulation was run
        if (!skip_file_writes)
        {
            char suffixed_filename[256];
            get_suffixed_filename("data/lowest_l_trajectories.dat", 1, suffixed_filename, sizeof(suffixed_filename));
            FILE *fp_lowest = fopen(suffixed_filename, "wb"); // Binary mode for fprintf_bin
            if (!fp_lowest)
            {
                fprintf(stderr, "Error: cannot open data/lowest_l_trajectories.dat\n");
                CLEAN_EXIT(1);
            }

            // Write data for Ntimes steps:
            for (int step = 0; step < Ntimes; step++)
            {
                double tval = step * dt;
                fprintf_bin(fp_lowest, "%f", tval);
                for (int p = 0; p < nlowest; p++)
                {
                    double rr = lowestL_r[p][step];
                    double Ecur = lowestL_E[p][step];
                    double lcur = lowestL_L[p][step];
                    fprintf_bin(fp_lowest, " %f %f %f", rr, Ecur, lcur);
                }
                fprintf_bin(fp_lowest, "\n");
            }
            fclose(fp_lowest);
        } // Close if (.skip_file_writes)
    }

    // Free lowest-L tracking arrays
    for (int p = 0; p < nlowest; p++)
    {
        free(lowestL_r[p]);
        free(lowestL_E[p]);
        free(lowestL_L[p]);
    }
    free(lowestL_r);
    free(lowestL_E);
    free(lowestL_L);



    int snapshot_steps[noutsnaps];
    // Determine the timesteps corresponding to the desired snapshot outputs
    for (int s = 0; s < noutsnaps; s++)
    {
        snapshot_steps[s] = (int)floor(
            s * (total_writes - 1) / (double)(noutsnaps - 1));
        log_message("DEBUG", "snapshot_steps[%d] = %d", s, snapshot_steps[s]);
    }

    if (g_doAllParticleData)
    {
        size_t total_outsnaps_size = (size_t)npts * noutsnaps; // Use size_t for large allocations
        log_message("INFO", "Allocating arrays for %zu particles across %d snapshots (total: %zu elements)",
                    (size_t)npts, noutsnaps, total_outsnaps_size);

        int *Rank_partdata_outsnaps = (int *)malloc(total_outsnaps_size * sizeof(int));
        float *R_partdata_outsnaps = (float *)malloc(total_outsnaps_size * sizeof(float));
        float *Vrad_partdata_outsnaps = (float *)malloc(total_outsnaps_size * sizeof(float));

        if (!Rank_partdata_outsnaps || !R_partdata_outsnaps || !Vrad_partdata_outsnaps)
        {
            printf("[ERROR] Failed to allocate memory for output arrays. Aborting.\n");
            fflush(stdout);
            CLEAN_EXIT(1);
        }

        char apd_filename_for_read[256];
        get_suffixed_filename("data/all_particle_data.dat", 1, apd_filename_for_read, sizeof(apd_filename_for_read));

        double fixed_bin_width = 0.0; // Initialized, set during first snapshot processing

        int start_index = 0; // Default start index for snapshot processing loop
        if (g_doRestart)
        {
            start_index = find_last_processed_snapshot(snapshot_steps, noutsnaps);
            if (start_index == -2) // Indicates all snapshots already processed
            {
                log_message("INFO", "All Rank files already exist for suffix '%s'. Skipping snapshot regeneration.", g_file_suffix);
                goto cleanup_partdata_outsnaps; // Jump past the processing loop
            }
            else if (start_index < 0) // Indicates error or no files found
            {
                start_index = 0; // Default to starting from the beginning
            }
            // Otherwise, start_index is the index of the first snapshot *to be* processed
            log_message("INFO", "Restart: Starting snapshot processing from index %d.", start_index);
        }

        printf("Starting parallel processing of %d snapshots (from %d to %d) using %d threads\n\n",
               noutsnaps - start_index, start_index, noutsnaps - 1, omp_get_max_threads());
        log_message("INFO", "Starting parallel processing with %d threads for %d snapshots (from index %d to %d)",
                    omp_get_max_threads(), noutsnaps - start_index, start_index, noutsnaps - 1);

#pragma omp parallel for schedule(static, 1) ordered
        for (int s = start_index; s < noutsnaps; s++)
        {
            int snap = snapshot_steps[s];

            float *tmpL_partdata_snap = (float *)malloc(npts * sizeof(float));
            int *tmpRank_partdata_snap = (int *)malloc(npts * sizeof(int));
            float *tmpR_partdata_snap = (float *)malloc(npts * sizeof(float));
            float *tmpV_partdata_snap = (float *)malloc(npts * sizeof(float));

            if (!tmpL_partdata_snap || !tmpRank_partdata_snap ||
                !tmpR_partdata_snap || !tmpV_partdata_snap)
            {
                fprintf(stderr, "Error: out of memory in parallel loop.\n");
                CLEAN_EXIT(1);
            }

            retrieve_all_particle_snapshot(
                apd_filename_for_read,
                snap,       // Which snap to read.
                npts,       // #particles
                block_size, // Block size.
                tmpL_partdata_snap,
                tmpRank_partdata_snap,
                tmpR_partdata_snap,
                tmpV_partdata_snap);

/**
 * Process snapshot data in order.
 * The #pragma omp ordered block ensures the copying into output arrays
 * (Rank_partdata_outsnaps, etc.) and status messages happen sequentially.
 */
#pragma omp ordered
            {
                // Display status message sequentially.
                printf("Processing snapshot %d...\n", snap);

                // Check array validity before copying
                if (!Rank_partdata_outsnaps || !R_partdata_outsnaps || !Vrad_partdata_outsnaps)
                {
                    log_message("ERROR", "Thread %d: Output arrays not properly allocated for snapshot %d",
                                omp_get_thread_num(), snap);
                }
                else if (!tmpRank_partdata_snap || !tmpR_partdata_snap || !tmpV_partdata_snap)
                {
                    log_message("ERROR", "Thread %d: Input arrays not properly allocated for snapshot %d",
                                omp_get_thread_num(), snap);
                }
                else
                {
                    // Copy data from thread-local buffers to the large shared output arrays.
                    for (int ii = 0; ii < npts; ii++)
                    {
                        if ((size_t)(s * npts + ii) < total_outsnaps_size)
                        {
                            Rank_partdata_outsnaps[s * npts + ii] = tmpRank_partdata_snap[ii];
                            R_partdata_outsnaps[s * npts + ii] = tmpR_partdata_snap[ii];
                            Vrad_partdata_outsnaps[s * npts + ii] = tmpV_partdata_snap[ii];
                        }
                        else
                        {
                            log_message("ERROR", "Thread %d: Index %zu out of bounds (%zu) in ordered copy for snapshot %d",
                                        omp_get_thread_num(), (size_t)(s * npts + ii), total_outsnaps_size, snap);
                            break;
                        }
                    }
                }
            }

            /** @brief Build "unsorted" arrays (Rank, Mass, R, Vrad, L) from local snapshot data. */
            log_message("INFO", "Thread %d: Allocating arrays for snapshot %d",
                        omp_get_thread_num(), snap);

            int *Rank_unsorted = (int *)malloc(npts * sizeof(int));
            double *Mass_unsorted = (double *)malloc(npts * sizeof(double));
            double *R_unsorted = (double *)malloc(npts * sizeof(double));
            double *Vrad_unsorted = (double *)malloc(npts * sizeof(double));
            double *L_unsorted = (double *)malloc(npts * sizeof(double));

            if (!Rank_unsorted || !Mass_unsorted || !R_unsorted || !Vrad_unsorted || !L_unsorted)
            {
                printf("[ERROR] Thread %d: Memory allocation failed for snapshot %d\n",
                       omp_get_thread_num(), snap);
                // Free any memory that was allocated.
                if (Rank_unsorted)
                    free(Rank_unsorted);
                if (Mass_unsorted)
                    free(Mass_unsorted);
                if (R_unsorted)
                    free(R_unsorted);
                if (Vrad_unsorted)
                    free(Vrad_unsorted);
                if (L_unsorted)
                    free(L_unsorted);
                continue; // Skip to next snapshot.
            }

            for (int ii = 0; ii < npts; ii++)
            {
                int rankval = tmpRank_partdata_snap[ii];
                double massv = (rankval + 1) * deltaM;

                Rank_unsorted[ii] = rankval;
                Mass_unsorted[ii] = massv;
                R_unsorted[ii] = (double)tmpR_partdata_snap[ii];
                Vrad_unsorted[ii] = (double)tmpV_partdata_snap[ii];
                L_unsorted[ii] = (double)tmpL_partdata_snap[ii];
            }

            // Sort by Rank
            // --- Prepare data structure for sorting by radius ---
            /** @brief Allocate `partarr` (array of struct PartData) for sorting snapshot data. */
            struct PartData *partarr = malloc(npts * sizeof(struct PartData));
            if (!partarr)
            {
                log_message("ERROR", "Thread %d: Failed to allocate partarr for snapshot %d",
                            omp_get_thread_num(), snap);
                // Free allocated memory.
                free(Rank_unsorted);
                free(Mass_unsorted);
                free(R_unsorted);
                free(Vrad_unsorted);
                free(L_unsorted);
                continue;
            }

            struct PartData *local_partarr = (struct PartData *)malloc(npts * sizeof(struct PartData));
            if (!local_partarr)
            {
                log_message("ERROR", "Thread %d: Failed to allocate local_partarr for snapshot %d",
                            omp_get_thread_num(), snap);
                free(partarr);
                continue;
            }

            for (int ii = 0; ii < npts; ii++)
            {
                // Check for invalid values to prevent sort issues.
                float rad_val = tmpR_partdata_snap[ii];
                float vrad_val = tmpV_partdata_snap[ii];
                float angmom_val = tmpL_partdata_snap[ii];

                // Handle potential NaN values before assigning to struct
                if (rad_val != rad_val) { // Check for NaN (safe with fast-math)
                    log_message("WARNING", "Thread %d: NaN radius at index %d for snapshot %d, replaced with 0",
                                omp_get_thread_num(), ii, snap);
                    rad_val = 0.0f;
                }
                if (vrad_val != vrad_val) vrad_val = 0.0f; // Check for NaN (safe with fast-math)
                if (angmom_val != angmom_val) angmom_val = 0.0f; // Check for NaN (safe with fast-math)

                local_partarr[ii].rank = tmpRank_partdata_snap[ii];
                local_partarr[ii].rad = rad_val;
                local_partarr[ii].vrad = vrad_val;
                local_partarr[ii].angmom = angmom_val;
                local_partarr[ii].original_index = ii; // Store original index before sort
            }

            memcpy(partarr, local_partarr, npts * sizeof(struct PartData));
            free(local_partarr);

            log_message("INFO", "Thread %d: Sorting partarr by radius for snapshot %d",
                        omp_get_thread_num(), snap);

            // Sort particle data by radius for snapshot processing
            sort_by_rad(partarr, npts);

            // If processing the first snapshot (index start_index), write lowest radius IDs
            if (s == start_index)
            {
                // Skip writing if file writes are disabled (e.g., restart post-processing only)
                if (!skip_file_writes)
                {
                    char fname_lowest_ids[256];
                    get_suffixed_filename("data/lowest_radius_ids.dat", 1, fname_lowest_ids, sizeof(fname_lowest_ids));
                    log_message("INFO", "Thread %d: Writing lowest radius IDs to %s",
                                omp_get_thread_num(), fname_lowest_ids);

                    FILE *id_file = fopen(fname_lowest_ids, "wb"); // Binary mode for fprintf_bin
                    if (id_file)
                    {
                        fprintf_bin(id_file, "# ID  Initial_Radius\n"); // Header (binary float for radius)

                        int num_tracked = (npts < 1000) ? npts : 1000; // Track up to 1000 particles

                        for (int i = 0; i < num_tracked; i++)
                        {
                            int original_id = partarr[i].original_index;
                            float initial_radius = partarr[i].rad;
                            fprintf_bin(id_file, "%d %f\n", original_id, initial_radius);
                        }
                        fclose(id_file);
                        log_message("INFO", "Thread %d: Successfully wrote %d lowest-radius particle IDs to %s",
                                    omp_get_thread_num(), num_tracked, fname_lowest_ids);
                    }
                    else
                    {
                        log_message("ERROR", "Thread %d: Failed to open %s for writing",
                                    omp_get_thread_num(), fname_lowest_ids);
                    }
                }
                else
                {
                    log_message("INFO", "Thread %d: Skipping lowest radius ID file creation in restart mode",
                                omp_get_thread_num());
                }
            }

            // Allocate arrays to store the data sorted by radius
            int *Rank_sorted = (int *)malloc(npts * sizeof(int));
            double *Mass_sorted = (double *)malloc(npts * sizeof(double));
            double *R_sorted = (double *)malloc(npts * sizeof(double));
            double *Vrad_sorted = (double *)malloc(npts * sizeof(double));
            double *L_sorted = (double *)malloc(npts * sizeof(double));

            if (!Rank_sorted || !Mass_sorted || !R_sorted || !Vrad_sorted || !L_sorted)
            {
                log_message("ERROR", "Thread %d: Failed to allocate sorted arrays for snapshot %d",
                            omp_get_thread_num(), snap);
                // Free allocated memory.
                if (Rank_sorted)
                    free(Rank_sorted);
                if (Mass_sorted)
                    free(Mass_sorted);
                if (R_sorted)
                    free(R_sorted);
                if (Vrad_sorted)
                    free(Vrad_sorted);
                if (L_sorted)
                    free(L_sorted);
                free(partarr);
                free(Rank_unsorted);
                free(Mass_unsorted);
                free(R_unsorted);
                free(Vrad_unsorted);
                free(L_unsorted);
                continue;
            }

            // Populate the sorted arrays using the sorted partarr
            for (int ii = 0; ii < npts; ii++)
            {
                partarr[ii].rank = ii; // Update rank based on sorted position
                int r_val = partarr[ii].rank; // Rank is now simply the index 'ii'
                Rank_sorted[ii] = r_val;
                Mass_sorted[ii] = (r_val + 1) * deltaM; // Mass enclosed up to this rank
                R_sorted[ii] = (double)partarr[ii].rad;
                Vrad_sorted[ii] = (double)partarr[ii].vrad;
                // Retrieve L from the unsorted array using the original index stored in partarr
                L_sorted[ii] = L_unsorted[partarr[ii].original_index];
            }

            /** @brief Set fixed bin width for density calculation using the first processed snapshot. */
            log_message("INFO", "Thread %d: Checking fixed bin width for snapshot %d",
                        omp_get_thread_num(), snap);
            if (s == start_index) // Only calculate on the first snapshot processed in this run
            {
                // Determine number of bins based on particle count (e.g., proportional to npts^(1/3))
                int num_bins_fixed = ceil(100 * pow(npts / 10000.0, 0.3333333333));
                if (num_bins_fixed < 1) num_bins_fixed = 1;
                // Use the maximum radius found in this first snapshot
                double r_max_first = R_sorted[npts - 1];
                if (r_max_first <= 0 || num_bins_fixed <= 0) {
                    log_message("ERROR", "Thread %d: Invalid r_max (%f) or bins (%d) for fixed_bin_width calc",
                                omp_get_thread_num(), r_max_first, num_bins_fixed);
                    // Set a default or handle error
                    fixed_bin_width = 1.0;
                } else {
                    fixed_bin_width = r_max_first / num_bins_fixed;
                }
                log_message("INFO", "Thread %d: Calculated fixed_bin_width=%f for snapshot %d (r_max=%f, bins=%d)",
                            omp_get_thread_num(), fixed_bin_width, snap, r_max_first, num_bins_fixed);
            }

            double *density_sorted = NULL;

            log_message("INFO", "Thread %d: Starting density calculation for snapshot %d",
                        omp_get_thread_num(), snap);
            {
                // Free previous allocation if it exists (safety measure)
                if (density_sorted != NULL)
                {
                    free(density_sorted);
                    density_sorted = NULL; // Avoid dangling pointer
                    log_message("INFO", "Thread %d: Freed previous density_sorted for snapshot %d",
                                omp_get_thread_num(), snap);
                }

                log_message("INFO", "Thread %d: Allocating density_sorted for snapshot %d",
                            omp_get_thread_num(), snap);
                density_sorted = malloc(npts * sizeof(double));
                if (!density_sorted)
                {
                    fprintf(stderr, "Error: Failed to allocate memory for density_sorted\n");
                    // Cleanup memory allocated within this snapshot's loop iteration
                    free(tmpL_partdata_snap); free(tmpRank_partdata_snap); free(tmpR_partdata_snap); free(tmpV_partdata_snap);
                    free(Rank_unsorted); free(Mass_unsorted); free(R_unsorted); free(Vrad_unsorted); free(L_unsorted);
                    free(partarr);
                    free(Rank_sorted); free(Mass_sorted); free(R_sorted); free(Vrad_sorted); free(L_sorted);
                    continue; // Proceed to the next snapshot
                }

                double bandwidth_factor = 0.08;
                log_message("INFO", "Thread %d: Checking R_sorted monotonicity for snapshot %d",
                            omp_get_thread_num(), snap);
                int r_violations = 0;
                for (int i = 1; i < npts; i++)
                {
                    if (R_sorted[i] <= R_sorted[i - 1])
                        r_violations++;
                }
                log_message("INFO", "Thread %d: Found %d radius violations for snapshot %d",
                            omp_get_thread_num(), r_violations, snap);

                double *R_filtered = NULL;
                double *Mass_filtered = NULL;
                int filtered_count = 0;

                if (r_violations > 0)
                {
                    // Allocate arrays to hold strictly monotonic radius/mass data
                    R_filtered = malloc(npts * sizeof(double));
                    Mass_filtered = malloc(npts * sizeof(double));
                    if (!R_filtered || !Mass_filtered) {
                        log_message("ERROR", "Thread %d: Failed to allocate filtered arrays for snapshot %d",
                            omp_get_thread_num(), snap);
                        continue;
                    }

                    // Copy only the points that maintain strict monotonicity
                    R_filtered[0] = R_sorted[0];
                    Mass_filtered[0] = Mass_sorted[0];
                    filtered_count = 1;

                    for (int i = 1; i < npts; i++)
                    {
                        if (R_sorted[i] > R_filtered[filtered_count - 1])
                        {
                            R_filtered[filtered_count] = R_sorted[i];
                            Mass_filtered[filtered_count] = Mass_sorted[i];
                            filtered_count++;
                        }
                    }
                }
                else
                {
                    // If no violations, point directly to the sorted data (no copy needed)
                    R_filtered = R_sorted;
                    Mass_filtered = Mass_sorted;
                    filtered_count = npts;
                }

                // Decimate the filtered data for spline interpolation efficiency
                int decimated_size = imin((int)ceil(npts / 100.0), 30000); // Limit max size
                if (decimated_size <= 0) decimated_size = 1; // Ensure at least one point
                double *R_decimated = malloc(decimated_size * sizeof(double));
                double *Mass_decimated = malloc(decimated_size * sizeof(double));
                if (!R_decimated || !Mass_decimated) {
                    log_message("ERROR", "Thread %d: Failed to allocate decimated arrays for snapshot %d",
                        omp_get_thread_num(), snap);
                    if (r_violations > 0) { free(R_filtered); free(Mass_filtered); }
                    continue;
                }

                if (filtered_count <= decimated_size)
                {
                    // Use all filtered points if fewer than desired decimated size
                    decimated_size = filtered_count;
                    memcpy(R_decimated, R_filtered, filtered_count * sizeof(double));
                    memcpy(Mass_decimated, Mass_filtered, filtered_count * sizeof(double));
                }
                else
                {
                    // Sample evenly from the filtered data
                    double step = (double)(filtered_count - 1) / (decimated_size - 1);
                    for (int i = 0; i < decimated_size; i++)
                    {
                        int idx = (int)round(i * step); // Use round for potentially better sampling
                        if (idx >= filtered_count)
                            idx = filtered_count - 1;
                        R_decimated[i] = R_filtered[idx];
                        Mass_decimated[i] = Mass_filtered[idx];
                    }
                }

                // Validate the range and finiteness of the decimated radius data for spline
                double min_r_check = (decimated_size > 0) ? R_decimated[0] : -1.0;
                double max_r_check = (decimated_size > 0) ? R_decimated[decimated_size - 1] : -1.0;
                // Check for invalid size, non-positive min, inverted range, or non-finite values (NaN or Inf)
                if (decimated_size <= 0 || min_r_check <= 0.0 || max_r_check <= min_r_check ||
                    ((min_r_check != min_r_check) || !(fabs(min_r_check) <= DBL_MAX)) || // Check !isfinite(min_r_check)
                    ((max_r_check != max_r_check) || !(fabs(max_r_check) <= DBL_MAX))    // Check !isfinite(max_r_check)
                   ) {
                    log_message("ERROR", "Thread %d: Invalid radius range [%f, %f] with %d points after decimation for snapshot %d",
                                omp_get_thread_num(), min_r_check, max_r_check, decimated_size, snap);
                    free(R_decimated); free(Mass_decimated);
                    if (r_violations > 0) { free(R_filtered); free(Mass_filtered); } // Free if allocated
                    free(density_sorted); density_sorted = NULL;
                    continue;
                }

                // Define the range and parameters for the uniform log-spaced grid
                double min_r = R_decimated[0];
                double max_r = R_decimated[decimated_size - 1] * 1.04; // Extend range slightly
                double log_min_r = log10(min_r);
                double log_max_r = log10(max_r);

                // Define size for the uniform grid (used for convolution)
                int grid_size = 131072; // Power of 2 often good for FFT, but direct used here
                if (grid_size <= 1) grid_size = 2; // Ensure at least 2 points
                double dlog = (log_max_r - log_min_r) / (grid_size - 1);

                // Allocate arrays for the uniform log-spaced grid
                double *r_grid = malloc(grid_size * sizeof(double));
                double *log_r_grid = malloc(grid_size * sizeof(double));
                double *mass_grid = malloc(grid_size * sizeof(double));
                double *density_grid = malloc(grid_size * sizeof(double));
                if (!r_grid || !log_r_grid || !mass_grid || !density_grid) {
                    log_message("ERROR", "Thread %d: Failed to allocate grid arrays for snapshot %d",
                        omp_get_thread_num(), snap);
                    // Free previously allocated resources
                    free(R_decimated); free(Mass_decimated);
                    if (r_violations > 0) { free(R_filtered); free(Mass_filtered); }
                    if (r_grid) free(r_grid);
                    if (log_r_grid) free(log_r_grid);
                    if (mass_grid) free(mass_grid);
                    if (density_grid) free(density_grid);
                    free(density_sorted); density_sorted = NULL;
                    continue;
                }

                // Populate the log-spaced grid coordinates
                for (int i = 0; i < grid_size; i++)
                {
                    log_r_grid[i] = log_min_r + i * dlog;
                    r_grid[i] = pow(10.0, log_r_grid[i]);
                }

                /** @brief Ensure r_grid is strictly monotonic for GSL spline init. */
                int rgrid_corrections_made = 0;
                for (int i = 1; i < grid_size; i++) {
                    if (r_grid[i] <= r_grid[i - 1]) {
                        // Ensure strict monotonicity with small absolute increment
                        r_grid[i] = r_grid[i - 1] + 1e-12;
                        rgrid_corrections_made++;
                    }
                }
                if (rgrid_corrections_made > 0) {
                     log_message("WARNING", "Thread %d made %d corrections to r_grid for monotonicity in snapshot %d", omp_get_thread_num(), rgrid_corrections_made, snap);
                }

                // Interpolate mass from decimated data onto the uniform log-spaced grid
                gsl_interp_accel *acc = NULL;
                gsl_spline *mass_spline = NULL;

#pragma omp critical(gsl_accel)
                {
                    acc = gsl_interp_accel_alloc();
                }

                if (!acc)
                {
                    log_message("ERROR", "Thread %d: Failed to allocate GSL interp accel for snapshot %d",
                                omp_get_thread_num(), snap);
                    continue;
                }

                /** @brief Ensure R_decimated is strictly monotonic for GSL spline init. */
                int corrections_made = 0;
                for (int i = 1; i < decimated_size; i++) {
                    if (R_decimated[i] <= R_decimated[i - 1]) {
                        // Ensure a strictly larger value using a small absolute increment
                        R_decimated[i] = R_decimated[i - 1] + 1e-12; // Add a small absolute value
                        corrections_made++;
                    }
                }
                if (corrections_made > 0) {
                     log_message("WARNING", "Thread %d made %d corrections to R_decimated for monotonicity in snapshot %d", omp_get_thread_num(), corrections_made, snap);
                }

#pragma omp critical(gsl_mass_spline)
                {
                    mass_spline = gsl_spline_alloc(gsl_interp_cspline, decimated_size);
                    if (mass_spline)
                    {
                        gsl_spline_init(mass_spline, R_decimated, Mass_decimated, decimated_size);
                    }
                }

                if (!mass_spline)
                {
                    log_message("ERROR", "Thread %d: Failed to allocate or initialize mass spline for snapshot %d",
                                omp_get_thread_num(), snap);
                    gsl_interp_accel_free(acc);
                    continue;
                }

                for (int i = 0; i < grid_size; i++)
                {
                    if (r_grid[i] < R_decimated[0])
                    {
                        mass_grid[i] = Mass_decimated[0]; // Extrapolate using first point
                    }
                    else if (r_grid[i] > R_decimated[decimated_size - 1])
                    {
                        mass_grid[i] = Mass_decimated[decimated_size - 1]; // Extrapolate using last point
                    }
                    else
                    {
#pragma omp critical(gsl_mass_eval)
                        {
                            // Use _e version for error checking if needed
                            mass_grid[i] = gsl_spline_eval(mass_spline, r_grid[i], acc);
                        }
                    }
                }

                // Calculate density = dM/dr on the uniform grid using central differences
                for (int i = 1; i < grid_size - 1; i++)
                {
                    // Avoid division by zero if grid points coincide
                    double dr = r_grid[i + 1] - r_grid[i - 1];
                    if (dr > 1e-15) {
                        double dM = mass_grid[i + 1] - mass_grid[i - 1];
                        density_grid[i] = dM / dr;
                    } else {
                        // Fallback for coincident points (e.g., use forward/backward difference or average)
                        density_grid[i] = (i > 1) ? density_grid[i-1] : 0.0; // Simple fallback
                    }

                    // Ensure non-negative density, apply floor
                    if (density_grid[i] < 1e-10)
                    {
                        density_grid[i] = 1e-10;
                    }
                }

                // Handle endpoints using one-sided difference or extrapolation
                if (grid_size >= 2) {
                   density_grid[0] = density_grid[1]; // Simple extrapolation
                   density_grid[grid_size - 1] = density_grid[grid_size - 2]; // Simple extrapolation
                } else if (grid_size == 1) {
                   density_grid[0] = 0.0; // Or some default
                }

                log_message("INFO", "Thread %d: Allocating memory for density smoothing arrays for snapshot %d",
                            omp_get_thread_num(), snap);
                double *density_smoothed_direct = malloc(npts * sizeof(double));
                double *direct_grid_result = malloc(grid_size * sizeof(double));

                if (!density_smoothed_direct || !direct_grid_result)
                {
                    log_message("ERROR", "Thread %d: Failed to allocate memory for smoothing arrays for snapshot %d",
                                omp_get_thread_num(), snap);
                    if (density_smoothed_direct)
                        free(density_smoothed_direct);
                    if (direct_grid_result)
                        free(direct_grid_result);
                    continue;
                }

                // Calculate smoothing kernel width (sigma) based on particle count
                double sigma = bandwidth_factor * pow(npts / 10000.0, -0.40);
                double sigma_log = sigma; // Apply smoothing in log-space

                // Perform Gaussian convolution on the uniform density grid
                if (!density_grid || !log_r_grid || !direct_grid_result)
                {
                    log_message("ERROR", "Thread %d: NULL arrays detected before gaussian_convolution for snapshot %d",
                                omp_get_thread_num(), snap);
                    // Clean up resources
                    if (density_grid) free(density_grid);
                    if (log_r_grid) free(log_r_grid);
                    if (direct_grid_result) free(direct_grid_result);
                    continue;
                }

                // Use critical section for FFTW thread safety if FFT is used internally
#pragma omp critical(gaussian_convolution)
                {
                    gaussian_convolution(density_grid, grid_size, log_r_grid, sigma_log, direct_grid_result);
                }

                // Create a GSL spline from the smoothed density on the uniform grid
                gsl_spline *density_spline = NULL;

#pragma omp critical(gsl_alloc)
                {
                    density_spline = gsl_spline_alloc(gsl_interp_cspline, grid_size);
                }

                if (!density_spline)
                {
                    log_message("ERROR", "Thread %d: Failed to allocate GSL spline for snapshot %d",
                                omp_get_thread_num(), snap);
                    // Clean up and skip.
                    free(density_smoothed_direct);
                    free(direct_grid_result);
                    continue;
                }

#pragma omp critical(gsl_init)
                {
                    gsl_spline_init(density_spline, r_grid, direct_grid_result, grid_size);
                }

                // Temporarily disable GSL default error handler to manage errors locally
                gsl_error_handler_t *old_handler = gsl_set_error_handler_off();

                // Interpolate smoothed density back onto the original sorted particle radii (R_sorted)
                for (int i = 0; i < npts; i++)
                {
                    double r_val = R_sorted[i]; // Target radius

                    if (r_val <= r_grid[0])
                    {
                        density_smoothed_direct[i] = direct_grid_result[0]; // Extrapolate below range
                    }
                    else if (r_val >= r_grid[grid_size - 1])
                    {
                        density_smoothed_direct[i] = direct_grid_result[grid_size - 1]; // Extrapolate above range
                    }
                    else
                    {
                        // For values within range, use spline interpolation with error checking.
                        int status = 0;
                        double result = 0.0;

// Try to evaluate with error checking.
#pragma omp critical(gsl_eval)
                        {
                            status = gsl_spline_eval_e(density_spline, r_val, acc, &result);
                        }

                        if (status != GSL_SUCCESS)
                        {
                            // Interp_errors++; // Commented out since variable is not used.

                            // Fall back to linear interpolation.
                            int idx_low = 0;

                            // Binary search to find the lower index.
                            int left = 0;
                            int right = grid_size - 1;

                            while (left <= right)
                            {
                                int mid = left + (right - left) / 2;

                                if (r_grid[mid] <= r_val && (mid == grid_size - 1 || r_grid[mid + 1] > r_val))
                                {
                                    idx_low = mid;
                                    break;
                                }
                                else if (r_grid[mid] > r_val)
                                {
                                    right = mid - 1;
                                }
                                else
                                {
                                    left = mid + 1;
                                }
                            }

                            int idx_high = idx_low + 1;

                            // Safety check.
                            if (idx_high >= grid_size)
                            {
                                idx_high = grid_size - 1;
                                idx_low = idx_high - 1;
                            }

                            // Linear interpolation.
                            double t = (r_val - r_grid[idx_low]) / (r_grid[idx_high] - r_grid[idx_low]);
                            result = direct_grid_result[idx_low] * (1.0 - t) + direct_grid_result[idx_high] * t;
                        }

                        density_smoothed_direct[i] = result;
                    }
                }

                // Restore the original error handler.
                gsl_set_error_handler(old_handler);

                // Copy the smoothed density values to the output array.
                memcpy(density_sorted, density_smoothed_direct, npts * sizeof(double));

                // Free resources.
                gsl_spline_free(mass_spline);
                gsl_spline_free(density_spline);
                gsl_interp_accel_free(acc);
                free(r_grid);
                free(log_r_grid);
                free(mass_grid);
                free(density_grid);
                free(direct_grid_result);
                free(density_smoothed_direct);
                free(R_decimated);
                free(Mass_decimated);

                // Free filtered arrays if they were allocated.
                if (r_violations > 0)
                {
                    free(R_filtered);
                    free(Mass_filtered);
                }
            }
            // =========================================================================
            // END DENSITY CALCULATION
            // =========================================================================

            // Density calculation complete using fixed binning from R_sorted min to max with (0,0) anchor point

            if (g_doDynPsi)
            {

                double *rrA = (double *)malloc(npts * sizeof(double));
                double *M_rA = (double *)malloc(npts * sizeof(double));

                for (int ii = 0; ii < npts; ii++)
                {
                    M_rA[ii] = Mass_sorted[ii];
                    rrA[ii] = R_sorted[ii];
                }

                double *inv_rA = (double *)malloc(npts * sizeof(double));
                for (int ii = 0; ii < npts; ii++)
                {
                    inv_rA[ii] = 1.0 / rrA[ii];
                }
                double *suffix_sum_invrA = (double *)malloc(npts * sizeof(double));
                suffix_sum_invrA[npts - 1] = inv_rA[npts - 1];
                for (int jj = npts - 2; jj >= 0; jj--)
                {
                    suffix_sum_invrA[jj] = suffix_sum_invrA[jj + 1] + inv_rA[jj];
                }

                double *psiAarr = (double *)malloc(npts * sizeof(double));
                for (int ii = 0; ii < npts; ii++)
                {
                    double rr = rrA[ii];
                    double M_r = M_rA[ii];
                    double outer_sum = 0.0;
                    if (ii < npts - 1)
                        outer_sum = (deltaM)*suffix_sum_invrA[ii + 1];
                    double psiA = G_CONST * ((M_r / rr) + outer_sum);
                    psiA *= VEL_CONV_SQ;
                    psiAarr[ii] = psiA;
                }

                // Compute PsiA(0).
                double sum_invr = 0.0;
                for (int ii = 0; ii < npts; ii++)
                {
                    sum_invr += inv_rA[ii];
                }
                double PsiA0 = G_CONST * deltaM * sum_invr * VEL_CONV_SQ; // PsiA at r=0.

                // No offset needed for PsiA, we keep psiAarr as is.
                // If Psi_theory_0 or offset is mentioned, we do nothing for PsiA.
                // For PsiB we do what original code says (it mentioned offset?).

                // Write PsiA and PsiB to files (as originally)
                char fname_PsiA[256];
                char base_filename[256];
                snprintf(base_filename, sizeof(base_filename), "data/Psi_methodA_t%05d.dat", snap);
                get_suffixed_filename(base_filename, 1, fname_PsiA, sizeof(fname_PsiA));

                FILE *fA = fopen(fname_PsiA, "wb"); // Binary mode for fprintf_bin
                // Write the zero point first.
                double rr_zero = 0.0;
                fprintf_bin(fA, "%f %f\n", rr_zero, PsiA0);
                for (int ii = 0; ii < npts; ii++)
                {
                    double rr = rrA[ii];
                    fprintf_bin(fA, "%f %f\n", rr, psiAarr[ii]);
                }
                fclose(fA);

                // Now create a new array with npts+1 points for PsiA to include (0,PsiA0) for spline
                double *rrA_spline = (double *)malloc((npts + 1) * sizeof(double));
                double *psiAarr_spline = (double *)malloc((npts + 1) * sizeof(double));
                rrA_spline[0] = 0.0;
                psiAarr_spline[0] = PsiA0;
                for (int ii = 0; ii < npts; ii++)
                {
                    rrA_spline[ii + 1] = rrA[ii];
                    psiAarr_spline[ii + 1] = psiAarr[ii];
                }

                // Start sketchy sorting hack.

                // Ensure strictly increasing values for spline initialization:
                for (int i = 1; i <= npts; i++)
                {
                    if (rrA_spline[i] <= rrA_spline[i - 1])
                    {
                        rrA_spline[i] = rrA_spline[i - 1] + 1e-5;
                    }
                }

                // End sketchy sorting hack.

                // Create a spline for PsiA(r) including the (0,PsiA0) point
                gsl_interp_accel *PsiAinterp = gsl_interp_accel_alloc();
                gsl_spline *splinePsiA = gsl_spline_alloc(gsl_interp_linear, npts + 1);

                gsl_spline_init(splinePsiA, rrA_spline, psiAarr_spline, npts + 1);

                // Compute PsiA for all unsorted and sorted radii using this spline
                double *PsiA_unsorted = (double *)malloc(npts * sizeof(double));
                for (int ii = 0; ii < npts; ii++)
                {
                    double rr_uf = R_unsorted[ii];
                    PsiA_unsorted[ii] = gsl_spline_eval(splinePsiA, rr_uf, PsiAinterp);
                }

                double *PsiA_sorted = (double *)malloc(npts * sizeof(double));
                for (int ii = 0; ii < npts; ii++)
                {
                    double rr_sf = R_sorted[ii];
                    PsiA_sorted[ii] = gsl_spline_eval(splinePsiA, rr_sf, PsiAinterp);
                }

                // =========================================================================
                // ENERGY CALCULATION FOR PARTICLE DISTRIBUTIONS
                // =========================================================================
                //
                // Calculate total energy for each particle using the formula:
                // E = Ψ - L²/(2r²) - v²/2
                //
                // Components:
                // - Ψ: Gravitational potential (from spline interpolation)
                // - L: Angular momentum
                // - r: Radial position
                // - v: Radial velocity
                //
                // Computed separately for both unsorted and sorted distributions
                double *E_unsorted = (double *)malloc(npts * sizeof(double));
                for (int ii = 0; ii < npts; ii++)
                {
                    double rr = R_unsorted[ii];      // Radius
                    double vrad = Vrad_unsorted[ii]; // Radial velocity
                    double l = L_unsorted[ii];       // Angular momentum

                    // Calculate total energy using the orbital energy equation
                    E_unsorted[ii] = PsiA_unsorted[ii]           // Potential energy (Ψ)
                                     - (l * l / (2.0 * rr * rr)) // Rotational energy (L²/2r²)
                                     - 0.5 * vrad * vrad;        // Kinetic energy (v²/2)
                }

                double *E_sorted = (double *)malloc(npts * sizeof(double));
                for (int ii = 0; ii < npts; ii++)
                {
                    double rr = R_sorted[ii];      // Radius
                    double vrad = Vrad_sorted[ii]; // Radial velocity
                    double l = L_sorted[ii];       // Angular momentum

                    // Same energy formula applied to sorted distribution
                    E_sorted[ii] = PsiA_sorted[ii]             // Potential energy (Ψ)
                                   - (l * l / (2.0 * rr * rr)) // Rotational energy (L²/2r²)
                                   - 0.5 * vrad * vrad;        // Kinetic energy (v²/2)
                }

                if (g_doDebug)
                {
                    // =========================================================================
                    // DEBUG ENERGY COMPUTATION - DYNAMIC ANALYSIS
                    // =========================================================================
                    //
                    // Computes dynamic energy for tracked particle before rewriting files.
                    // Uses PsiA_unsorted to calculate energy for the specific debug ID.
                    // This information is used to validate energy conservation during simulation.

                    {
                        /**
                         * TRACKED PARTICLE ENERGY ANALYSIS
                         *
                         * Extract and compute dynamic energy components for a single tracked particle
                         * identified by DEBUG_PARTICLE_ID. Used to validate energy conservation
                         * and compare with theoretical predictions. The components are stored for
                         * post-simulation analysis and visualization.
                         */

                        // Retrieve the particle ID to track from the global constant
                        int debug_id = DEBUG_PARTICLE_ID;

                        // Extract physical properties from the unsorted arrays (original ordering)
                        double r_val = R_unsorted[debug_id];       // Radius
                        double v_val = Vrad_unsorted[debug_id];    // Radial velocity
                        double l_val = L_unsorted[debug_id];       // Angular momentum
                        double psiA_val = PsiA_unsorted[debug_id]; // Potential (newly computed)

                        // Calculate energy components:
                        // Total energy = potential - kinetic
                        double E_dyn = psiA_val - 0.5 * (v_val * v_val + (l_val * l_val) / (r_val * r_val));

                        // Kinetic energy = radial + rotational components
                        double K_dyn = 0.5 * (v_val * v_val + (l_val * l_val) / (r_val * r_val));

                        // Calculate simulation time corresponding to this snapshot
                        // Ensure consistent time basis with approximate energy calculation
                        double sim_time = (double)(snapshot_steps[s]) * dtwrite * dt;

                        // Store energy components in global arrays for later analysis
                        store_debug_dynE_components(
                            s,        // Snapshot index
                            E_dyn,    // Total energy
                            K_dyn,    // Kinetic energy
                            psiA_val, // Potential energy
                            sim_time, // Simulation time
                            r_val     // Radius
                        );

                        // =========================================================================
                        // END DEBUG ENERGY COMPUTATION
                        // =========================================================================
                    }
                }

                log_message("INFO", "Thread %d: Starting to write rank files for snapshot %d",
                            omp_get_thread_num(), snap);

                if (g_doDynRank)
                {
                    log_message("INFO", "Thread %d: About to write unsorted Rank file for snapshot %d",
                                omp_get_thread_num(), snap);

                    // Unsorted file write
                    char fname_unsorted[256];
                    char base_filename[256];
                    snprintf(base_filename, sizeof(base_filename), "data/Rank_Mass_Rad_VRad_unsorted_t%05d.dat", snap);
                    get_suffixed_filename(base_filename, 1, fname_unsorted, sizeof(fname_unsorted));

                    log_message("INFO", "Thread %d: Opening %s for writing",
                                omp_get_thread_num(), fname_unsorted);

                    FILE *fun_final = fopen(fname_unsorted, "wb"); // Binary mode for fprintf_bin

                    if (!fun_final)
                    {
                        log_message("ERROR", "Thread %d: Failed to open %s for writing",
                                    omp_get_thread_num(), fname_unsorted);
                        // Skip this file but continue with other operations.
                        goto skip_unsorted_write;
                    }
                    for (int ii = 0; ii < npts; ii++)
                    {
                        fprintf_bin(fun_final, "%d %f %f %f %f %f %f\n",
                                    Rank_unsorted[ii],
                                    Mass_unsorted[ii],
                                    R_unsorted[ii],
                                    Vrad_unsorted[ii],
                                    PsiA_unsorted[ii],
                                    E_unsorted[ii],
                                    L_unsorted[ii]);
                    }
                    fclose(fun_final);
                    log_message("INFO", "Thread %d: Successfully wrote unsorted Rank file for snapshot %d",
                                omp_get_thread_num(), snap);

                skip_unsorted_write:

                    // Sorted file write
                    log_message("INFO", "Thread %d: About to write sorted Rank file for snapshot %d",
                                omp_get_thread_num(), snap);

                    char fname_sorted[256];
                    snprintf(base_filename, sizeof(base_filename), "data/Rank_Mass_Rad_VRad_sorted_t%05d.dat", snap);
                    get_suffixed_filename(base_filename, 1, fname_sorted, sizeof(fname_sorted));

                    log_message("INFO", "Thread %d: Opening %s for writing",
                                omp_get_thread_num(), fname_sorted);

                    FILE *fsort_final = fopen(fname_sorted, "wb"); // Binary mode for fprintf_bin

                    if (!fsort_final)
                    {
                        log_message("ERROR", "Thread %d: Failed to open %s for writing",
                                    omp_get_thread_num(), fname_sorted);
                        // Skip this file but continue.
                        goto skip_sorted_write;
                    }

                    for (int ii = 0; ii < npts; ii++)
                    {
                        fprintf_bin(fsort_final, "%d %f %f %f %f %f %f %f\n",
                                    Rank_sorted[ii],
                                    Mass_sorted[ii],
                                    R_sorted[ii],
                                    Vrad_sorted[ii],
                                    PsiA_sorted[ii],
                                    E_sorted[ii],
                                    L_sorted[ii],
                                    density_sorted[ii]);
                    }
                    fclose(fsort_final);
                    log_message("INFO", "Thread %d: Successfully wrote sorted Rank file for snapshot %d",
                                omp_get_thread_num(), snap);

                skip_sorted_write:; // Empty statement needed after label.
                } // END if (g_doDynRank).

                log_message("INFO", "Thread %d: Completed writing rank files for snapshot %d, freeing memory",
                            omp_get_thread_num(), snap);

                free(inv_rA);
                free(suffix_sum_invrA);
                free(rrA);
                free(M_rA);
                free(psiAarr);
                free(rrA_spline);
                free(psiAarr_spline);
                free(PsiA_unsorted);
                free(PsiA_sorted);
                free(E_unsorted);
                free(E_sorted);
                gsl_spline_free(splinePsiA);
                gsl_interp_accel_free(PsiAinterp);
            }
            // Free temporary arrays

            // =========================================================================
            // CLEANUP SECTION - PARTICLE DATA PROCESSING COMPLETE
            // =========================================================================
            //
            // Release all memory allocated during particle data processing phase.
            // Memory is freed in the reverse order of allocation to prevent memory leaks.
            // This includes thread-local arrays and shared memory structures.
            // Free thread-local arrays
            // Free all local arrays
            free(Rank_unsorted);
            free(Mass_unsorted);
            free(R_unsorted);
            free(Vrad_unsorted);
            free(L_unsorted);

            free(Rank_sorted);
            free(Mass_sorted);
            free(R_sorted);
            free(Vrad_sorted);
            free(L_sorted);

            free(partarr);

            // Free density_sorted array which was allocated during density calculation
            if (density_sorted != NULL) {
                free(density_sorted);
                log_message("INFO", "Thread %d: Freed density_sorted for snapshot %d",
                           omp_get_thread_num(), snap);
            }

            // Free the initial local buffers
            log_message("INFO", "Thread %d: FINAL cleanup for snapshot %d",
                        omp_get_thread_num(), snap);

            free(tmpL_partdata_snap);
            free(tmpRank_partdata_snap);
            free(tmpR_partdata_snap);
            free(tmpV_partdata_snap);

            log_message("INFO", "Thread %d COMPLETED processing of snapshot index %d (snapshot number %d)",
                        omp_get_thread_num(), s, snap);
        }

        log_message("INFO", "All snapshot processing completed. Cleaning up.");

    cleanup_partdata_outsnaps:
        // Free the memory for particle data arrays
        free(Rank_partdata_outsnaps);
        free(R_partdata_outsnaps);
        free(Vrad_partdata_outsnaps);
    }
    {
        log_message("INFO", "Final operations using file suffix: %s", g_file_suffix);

        for (int s = 0; s < noutsnaps; s++)
        {
            snapshot_steps[s] = (int)floor(
                s * (total_writes - 1) / (double)(noutsnaps - 1));
        }
    }

    free(inverse_map);

    for (i = 0; i < num_traj_particles; i++)
    {
        free(trajectories[i]);
        free(energies[i]);
        free(velocities_arr[i]);
        free(mu_arr[i]);
        free(E_arr[i]);
        free(L_arr[i]);
    }
    free(trajectories);
    free(energies);
    free(velocities_arr);
    free(mu_arr);
    free(E_arr);
    free(L_arr);

    free(E_i_arr);
    free(L_i_arr);

    free(r_initial);
    free(rv_initial);
    free(v_initial);
    free(l_initial);
    free(r_final);
    free(v_final);


    free(mass);
    free(radius);
    if (radius_monotonic_grid_nfw != NULL) {
        free(radius_monotonic_grid_nfw);
        radius_monotonic_grid_nfw = NULL;
    }
    free(Psivalues);
    free(nPsivalues);
    free(innerintegrandvalues);
    free(Evalues);
    if (w != NULL) {
        gsl_integration_workspace_free(w);
    }
    // Free GSL RNG resources
    if (g_rng != NULL) {
        gsl_rng_free(g_rng);
        g_rng = NULL;
    }

    // End method_select == 3 block.

    // Move the freeing of particles *outside* the if-block.
    for (i = 0; i < 5; i++)
    {
        free(particles[i]);
    }
    free(particles);

    if (g_doDebug)
    {
        finalize_debug_energy_output(); // Ensures all data is collected first
    }

    gsl_spline_free(splinemass);
    gsl_spline_free(splinePsi);
    gsl_spline_free(splinerofPsi);
    gsl_interp_accel_free(enclosedmass);
    gsl_interp_accel_free(Psiinterp);
    gsl_interp_accel_free(rofPsiinterp);
    gsl_interp_free(g_main_fofEinterp);
    gsl_interp_accel_free(g_main_fofEacc);

    free_local_snap_arrays();
    cleanup_all_particle_data();

#ifdef _OPENMP
    // Clean up FFTW threads only if they were initialized.
    fftw_cleanup_threads();

    if (g_enable_sidm_scattering) {
        printf("Total SIDM scattering events during simulation: %lld\n", g_total_sidm_scatters);
        log_message("INFO", "Total SIDM scattering events: %lld", g_total_sidm_scatters);
    }
#endif

    // Free per-thread GSL RNG resources
    if (g_rng_per_thread != NULL) {
        for (int i_rng = 0; i_rng < g_max_omp_threads_for_rng; ++i_rng) {
            if (g_rng_per_thread[i_rng] != NULL) {
                gsl_rng_free(g_rng_per_thread[i_rng]);
            }
        }
        free(g_rng_per_thread);
        g_rng_per_thread = NULL;
        log_message("INFO", "Freed per-thread GSL RNGs.");
    }

    // Cleanup for the conditionally declared persistent sort buffer.
    if (g_sort_columns_buffer != NULL) {
        for (int i = 0; i < g_sort_columns_buffer_npts; i++) {
            if (g_sort_columns_buffer[i]) free(g_sort_columns_buffer[i]);
        }
        free(g_sort_columns_buffer);
        g_sort_columns_buffer = NULL; // Mark as freed.
        g_sort_columns_buffer_npts = 0; // Reset size.
    }

    // Free global NFW spline resources if they were allocated for NFW profile
    if (g_nfw_splinemass_for_force) {
        gsl_spline_free(g_nfw_splinemass_for_force);
        g_nfw_splinemass_for_force = NULL;
        log_message("INFO", "Freed global NFW mass spline for force calculation.");
    }
    if (g_nfw_enclosedmass_accel_for_force) {
        gsl_interp_accel_free(g_nfw_enclosedmass_accel_for_force);
        g_nfw_enclosedmass_accel_for_force = NULL;
        log_message("INFO", "Freed global NFW mass spline accelerator.");
    }

    // Free particle scatter state array
    free(g_particle_scatter_state);

    return 0;
} // End main function.
