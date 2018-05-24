#include <stdbool.h>
#include <zconf.h>
#include "testing.h"
#include "../helper.h"
#include "../logger.h"

void assertEqual(double value, double expectation)
{
    if (value != expectation)
    {
        logTestEvent(ERROR,"assertEqual: %f != %f", value, expectation);
        ERROR("AssertionError");
    }
}

int expectEqual(double value, double expectation, char *errorMsgFmt, ...)
{
    if (value != expectation)
    {
        char buf[1024], buf2[512];
        va_list localArgs;
        va_start(localArgs, errorMsgFmt);
        vsprintf(buf,errorMsgFmt,localArgs);
        va_end(localArgs);
        sprintf(buf2, "%.1f != %.1f", value, expectation);
        strcat(buf,buf2);
//        logTestEvent(WARNING,"expectEqual: %f != %f", value, expectation);
        logTestEvent(WARNING, buf);
        return 1;
    }
    else
        return 0; // Return if fails
}

//eof


