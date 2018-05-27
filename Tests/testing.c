#include <stdbool.h>
#include <zconf.h>
#include "testing.h"
#include "../helper.h"
#include "../logger.h"

bool acceptablyEqual(double a, double b)
{
    double eps = 1e-8;
    return (fabs(a-b) < eps);
}

void assertEqual(double value, double expectation, const char *testName)
{
    if (!acceptablyEqual(value, expectation))
    {
        logTestEvent(ERROR, testName, "assertEqual: %f != %f", value, expectation);
        THROW_ERROR("AssertionError");
    }
}

int expectEqual(double value, double expectation, const char *testName, char *errorMsgFmt, ...)
{
    if (!acceptablyEqual(value, expectation))
    {
        char buf[1024], buf2[512];
        va_list localArgs;
        va_start(localArgs, errorMsgFmt);
        vsprintf(buf, errorMsgFmt, localArgs);
        va_end(localArgs);
        sprintf(buf2, "%.1f != %.1f", value, expectation);
        strcat(buf, buf2);
//        logTestEvent(WARNING,"expectEqual: %f != %f", value, expectation);
        logTestEvent(WARNING, testName, buf);
        return 1;
    }
    else
    {
        return 0;
    } // Return if fails
}

//eof


