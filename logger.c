//
// Created by tommaso on 05/05/18.
//

#include "logger.h"
#include "timing.h"
#include "helper.h"
#include <stdio.h>
#include <stdarg.h>
#include <memory.h>
#include <sys/stat.h>

/*
 * Logging machinery - how to:
 * 1) Make sure to call openLogFile() once at the beginning of main(), before calling any logmsg().
 * 2) To log a message both to log file and to console, call logmsg() as follows:
 *      logEvent(time, message, [arg1, arg2, ...]);
 *   time must be the current simulation time, message a printf-friendly string,
 *   then optional args in prtf style can be passed.
 * 3) Make sure to call closeLogfile() before exiting the main (this closes the file at OS level).
 */

static const char *DEBUG_STR[] = {
        FOREACH_DEBUG(GENERATE_STRING)
};

static long LOGGER_START_TIME = 0;
static DebugLevel DEBUG_LEVEL = INFO; // Default level
static char LOG_FILE_FOLDER[512] = "./";
static char LOG_FILE_NAME[128] = "sim.log";
static char LOG_FILE_FULL_PATH[512] = "";
static FILE* LOG_FILE;

void setLoggerStartTime()
{
    LOGGER_START_TIME = getCurrentTimeMillis();
}

void setLoggerOutputFolder(const char *outputFolder)
{
    strcpy(LOG_FILE_FOLDER, outputFolder);
    // Create folder if not on filesystem
    struct stat st = {0};
    if (stat(outputFolder, &st) == -1) {
        mkdir(outputFolder, 0700);
    }
    sprintf(LOG_FILE_FULL_PATH, "%s/%s", LOG_FILE_FOLDER, LOG_FILE_NAME);
}

void setLoggerFileName(const char *outputFile)
{
    strcpy(LOG_FILE_NAME, outputFile);
    sprintf(LOG_FILE_FULL_PATH, "%s/%s", LOG_FILE_FOLDER, LOG_FILE_NAME);
}

void setLoggerDebugLevel(DebugLevel debugLevel)
{
    DEBUG_LEVEL = debugLevel;
}

DebugLevel getLoggerDebugLevel()
{
    return DEBUG_LEVEL;
}

void openLogFile()
{
    if (strlen(LOG_FILE_FULL_PATH) == 0)
    {
        sprintf(LOG_FILE_FULL_PATH, "%s/%s", LOG_FILE_FOLDER, LOG_FILE_NAME);
    }
    LOG_FILE = fopen(LOG_FILE_FULL_PATH, "w");
    if (LOG_FILE == NULL)
    {
        char buf[512];
        sprintf(buf, "Cannot open file %s", LOG_FILE_FULL_PATH);
        ERROR(buf);
    }
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

void logEvent(DebugLevel eventDebugLevel, double t, char *fmt, ...)
{
    // Newline at the end of the message is included.
    // If this trace is too low level for current debug level, skip it.
    if (eventDebugLevel < DEBUG_LEVEL)
        return;
    //
    double timestamp = getTimeSpentSeconds(LOGGER_START_TIME, getCurrentTimeMillis());
    va_list args;
    va_start(args,fmt);
    printf("RUN_t: [%06.3f] SIM_t: [%012.9f] %s: ", timestamp, t, DEBUG_STR[eventDebugLevel]);
    //printf("RUNTIME: [%06.3f] %s: ", timestamp, DEBUG_STR[eventDebugLevel]);
    vprintf(fmt, args);
    //printf(" t = %012.9f", t);
    printf("\n");
    va_end(args);
    va_start(args,fmt);
    fprintf(LOG_FILE, "RUN_t: [%06.3f] SIM_t: [%012.9f] %s: ", timestamp, t, DEBUG_STR[eventDebugLevel]);
    //fprintf(LOG_FILE, "RUNTIME: [%06.3f] %s: ", timestamp, DEBUG_STR[eventDebugLevel]);
    vfprintf(LOG_FILE, fmt, args);
    //fprintf(LOG_FILE, " t = %012.9f", t);
    fprintf(LOG_FILE, "\n");
    va_end(args);
}

void logTestEvent(DebugLevel eventDebugLevel, const char *testName, char *fmt, ...)
{
    // Newline at the end of the message is included.
    // If this trace is too low level for current debug level, skip it.
    if (eventDebugLevel < DEBUG_LEVEL)
        return;
    //
    double timestamp = getTimeSpentSeconds(LOGGER_START_TIME, getCurrentTimeMillis());
    va_list args;
    va_start(args,fmt);
    printf("[%06.3f] %s: [%s] ", timestamp, DEBUG_STR[eventDebugLevel], testName);
    vprintf(fmt, args);
    printf("\n");
    va_end(args);
    va_start(args,fmt);
    fprintf(LOG_FILE, "[%06.3f] %s: [%s] ", timestamp, DEBUG_STR[eventDebugLevel], testName);
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

void flushLogFile()
{
    fflush(LOG_FILE);
}


void closeLogFile()
{
    fclose(LOG_FILE);
}

//eof
