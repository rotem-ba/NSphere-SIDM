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
 #include <stdarg.h>
 #include <time.h>
 #include <sys/stat.h>
 #include "globals.h"
 #include "particle_array_ops.h"

/**
 * @brief Writes a formatted message to the log file with timestamp and severity level.
 *
 * Parameters
 * ----------
 * level : const char*
 *     Severity level (e.g., "INFO", "WARNING", "ERROR").
 * format : const char*
 *     Printf-style format string.
 * ... :
 *     Variable arguments for the format string.
 *
 * Returns
 * -------
 * None
 *
 * @note Creates the "log" directory if it doesn't exist.
 * @warning Prints an error to stderr if the log file cannot be opened.
 *          Logging only occurs if the global `g_enable_logging` flag is set.
 * @see g_enable_logging
 */
void log_message(const char *level, const char *format, ...)
{
    // Only write to log file if logging is enabled
    if (g_enable_logging)
    {
        // Create log directory if it doesn't exist
        struct stat st = {0};
        if (stat("log", &st) == -1)
        {
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
            mkdir("log"); // Windows
#else
            mkdir("log", 0755); // Unix-like systems
#endif
        }

        // Always use a single log file regardless of file suffix
        const char *log_filename = "log/nsphere.log";

        FILE *logfile = fopen(log_filename, "a");
        if (logfile)
        {
            time_t now;
            time(&now);
            char timestamp[64];
            strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", localtime(&now));

            // Include file suffix in log entries if available
            if (g_file_suffix[0] != '\0')
            {
                fprintf(logfile, "[%s] [%s] [%s] ", timestamp, level, g_file_suffix);
            }
            else
            {
                fprintf(logfile, "[%s] [%s] ", timestamp, level);
            }

            va_list args;
            va_start(args, format);
            vfprintf(logfile, format, args);
            va_end(args);

            fprintf(logfile, "\n");
            fclose(logfile);
        }
        else
        {
            // Print an error message to stderr if the log file cannot be opened
            fprintf(stderr, "Warning: Failed to open log file '%s'\n", log_filename);
        }
    }
}

/**
 * @brief Display final parameter values used for the simulation run
 */
void print_input_parameters() {
    printf("Parameter values requested:\n\n");
    printf("  Number of Particles:          %d\n", npts);
    printf("  Number of Time Steps:         %d\n", Ntimes);
    printf("  Number of Dynamical Times:    %d\n", tfinal_factor);
    printf("  Number of Output Snapshots:   %d\n", nout);
    printf("  Steps Between Writes:         %d\n", dtwrite);
    printf("  Tidal Stripping Fraction:     %.5f\n", tidal_fraction);
    printf("  Integration Method:           %d (%s)\n", method_select, method_name);
    printf("  Sorting Algorithm:            %d (%s)\n", display_sort, get_sort_description(g_defaultSortAlg));
    // printf("  Filename Tag:                 %s\n", filename_tag[0] ? filename_tag : "[none]");
    printf("  SIDM Scattering:              %s\n", g_enable_sidm_scattering ? "Enabled via --sidm" : "Disabled (Default)");
    printf("  SIDM Execution Mode:          %s\n", g_sidm_execution_mode == 1 ? "Parallel (Default)" : "Serial");
    printf("  SIDM Opacity Kappa:           %.1f cm^2/g (Default: 50.0, User set: %s)\n", g_sidm_kappa, g_sidm_kappa_provided ? "Yes" : "No");

}

/**
 * @brief Display density parameters
 */
void print_density_params(){
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
}

/**
 * @brief Display logging status based on g_enable_logging flag.
 */
void print_logging_status(){
    if (g_enable_logging) {
        printf("  Logging:                      Enabled (log/nsphere.log)\n\n");
        log_message("INFO", "Simulation started with %d particles, %d timesteps, %d dynamical times", npts, Ntimes, tfinal_factor);
    } else {
        printf("  Logging:                      Disabled\n\n");
    }
}

/**
 * @brief Display check for SIDM + parallel mode without OpenMP.
 */
void print_SIDM_OpenMP_status(){
    if (g_enable_sidm_scattering && g_sidm_execution_mode == 1) {
        #ifndef _OPENMP
            printf("Warning: SIDM parallel mode is enabled by default, but OpenMP is not available in this build.\n");
            printf("         SIDM will run serially. Use '--sidm-mode serial' to suppress this warning.\n\n");
            log_message("WARNING", "SIDM parallel mode requested but OpenMP not available, will run serially.");
            g_sidm_execution_mode = 0; // Force serial if no OpenMP
        #endif
    }
}

/**
 * @brief Display warning that OpenMP is not available
 */
void warn_no_OpenMP(){
    /** @warning OpenMP section: Warns user when compiled without OpenMP support. */
    printf("WARNING: OpenMP is NOT ENABLED in this build!\n");
    printf("This will result in significantly reduced performance.\n");
    printf("For better performance, please install OpenMP and recompile with -fopenmp flag.\n\n\n");

    log_message("WARNING", "OpenMP not available - running in single-threaded mode");
}
