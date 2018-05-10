//
// Created by tommaso on 05/05/18.
//

#include "logger.h"
#include "timing.h"
#include <stdio.h>
#include <stdarg.h>
#include <memory.h>

/*
 * Logging machinery - how to:
 * 1) Make sure to call openLogFile() once at the beginning of main(), before calling any logmsg().
 * 2) To log a message both to log file and to console, call logmsg() as follows:
 *      logEvent(time, message, [arg1, arg2, ...]);
 *   time must be the current simulation time, message a printf-friendly string,
 *   then optional args in prtf style can be passed.
 * 3) Make sure to call closeLogfile() before exiting the main (this closes the file at OS level).
 */

// TODO: logfile location should be configurable to allow for easy running of tests in batches
static long LOGGER_START_TIME = 0;
static char LOG_FILE_FOLDER[512] = "./";
static char* LOG_FILE_NAME = "sim.log";
static char LOG_FILE_FULL_PATH[512] = "";
static FILE* LOG_FILE;

void setLoggerStartTime()
{
    LOGGER_START_TIME = getCurrentTimeMillis();
}

void setLoggerOutputFolder(const char *outputFolder)
{
    strcpy(LOG_FILE_FOLDER, outputFolder);
    sprintf(LOG_FILE_FULL_PATH, "%s/%s", LOG_FILE_FOLDER, LOG_FILE_NAME);
}

void openLogFile()
{
    if (strlen(LOG_FILE_FULL_PATH) == 0)
    {
        sprintf(LOG_FILE_FULL_PATH, "%s/%s", LOG_FILE_FOLDER, LOG_FILE_NAME);
    }
    LOG_FILE = fopen(LOG_FILE_FULL_PATH, "w");
}

void logRawString(char *fmt, ...)
{
    // Newline at the end of the message is included.
    double timestamp = getTimeSpentSeconds(LOGGER_START_TIME, getCurrentTimeMillis());
    va_list args;
    va_start(args,fmt);
    printf("[%06.3f] ", timestamp);
    vprintf(fmt, args);
    va_end(args);
    va_start(args,fmt);
    fprintf(LOG_FILE, "[%06.3f] ", timestamp);
    vfprintf(LOG_FILE, fmt, args);
    va_end(args);
}

void logEvent(double t, char *fmt, ...)
{
    // Newline at the end of the message is included.
    double timestamp = getTimeSpentSeconds(LOGGER_START_TIME, getCurrentTimeMillis());
    va_list args;
    va_start(args,fmt);
    printf("[%06.3f] ", timestamp);
    printf("[%012.9f] ", t);
    vprintf(fmt, args);
    printf("\n");
    va_end(args);
    va_start(args,fmt);
    fprintf(LOG_FILE, "[%06.3f] ", timestamp);
    fprintf(LOG_FILE, "[%012.9f] ", t);
    vfprintf(LOG_FILE, fmt, args);
    fprintf(LOG_FILE, "\n");
    va_end(args);
}

void logMsg(char *fmt, ...)
{
    // Newline at the end of the message is included.
    double timestamp = getTimeSpentSeconds(LOGGER_START_TIME, getCurrentTimeMillis());
    va_list args;
    va_start(args,fmt);
    printf("[%06.3f] ", timestamp);
    printf("---> ");
    vprintf(fmt, args);
    printf("\n");
    va_end(args);
    va_start(args,fmt);
    fprintf(LOG_FILE, "[%06.3f] ", timestamp);
    fprintf(LOG_FILE, "---> ");
    vfprintf(LOG_FILE, fmt, args);
    fprintf(LOG_FILE, "\n");
    va_end(args);
}

void closeLogFile()
{
    fclose(LOG_FILE);
}

//eof
