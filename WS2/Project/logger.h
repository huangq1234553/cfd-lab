//
// Created by tommaso on 05/05/18.
//

#ifndef SIM_LOGGER_H
#define SIM_LOGGER_H

// This allows for automatically getting strings of debug levels (see https://stackoverflow.com/a/10966395 )
// NOTE: order is important for correctly managing incremental levels of debug.
#define FOREACH_DEBUG(DLEVEL) \
    DLEVEL(DEBUG) \
    DLEVEL(INFO) \
    DLEVEL(PRODUCTION) \
    DLEVEL(WARNING) \
    DLEVEL(ERROR) \

#define GENERATE_ENUM(ENUM) ENUM,
#define GENERATE_STRING(STRING) #STRING,

typedef enum DebugLevel
{
    FOREACH_DEBUG(GENERATE_ENUM)
} DebugLevel;

void setLoggerStartTime();

void setLoggerOutputFolder(const char *outputFolder);

void setLoggerFileName(const char *outputFile);

void setLoggerDebugLevel(DebugLevel debugLevel);

DebugLevel getLoggerDebugLevel();

void openLogFile();

void logRawString(char *fmt, ...);

void logEvent(DebugLevel eventDebugLevel, double t, char *fmt, ...);

void logTestEvent(DebugLevel eventDebugLevel, char *testName, char *fmt, ...);

void logMsg(char *fmt, ...);

void closeLogFile();

#endif //SIM_LOGGER_H
