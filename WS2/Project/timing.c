//
// Created by tommaso on 10/05/18.
//

#include <sys/time.h>
#include <stdlib.h>
#include "timing.h"

double getTimeSpentSeconds(long startTimeMillis, long endTimeMillis)
{ return (endTimeMillis - startTimeMillis) / 1000.0; }

long getCurrentTimeMillis()
{
    struct timeval timecheck;
    gettimeofday(&timecheck, NULL);
    long currentTimeMillis = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;
    return currentTimeMillis;
}
