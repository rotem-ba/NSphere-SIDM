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
 #include <stdio.h>
 #include <stdarg.h>
 #include <time.h>
 #include <sys/stat.h>
 #include "globals.h"

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
