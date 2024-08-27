#include <sys/time.h>
#include <ctime>
#include "Bench.h"
#include <iostream>

/* Taken from SO
http://stackoverflow.com/questions/1861294/how-to-calculate-execution-time-of-a-code-snippet-in-c. Answer
from Andreas Bonini */

/* Remove if already defined */
typedef long long int64;
typedef unsigned long long uint64;

/* Returns the amount of milliseconds elapsed since the UNIX epoch.*/

uint64
Bench::GetTimeMs64() const
{
    /* Linux */
    struct timeval tv;

    gettimeofday(&tv, NULL);

    uint64 ret = tv.tv_usec;
    /* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
    ret /= 1000;

    /* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
    ret += (tv.tv_sec * 1000);

    return ret;
}

uint64
Bench::TimeElapsed() const
{
    uint64 start = GetTimeMs64();
    fFunction();
    return GetTimeMs64() - start;
}

long double
Bench::TimeRepetitions(unsigned iterations_)
{
    fIterations = iterations_;
    fMeanTime = 0;
    for (unsigned i = 0; i < iterations_; i++)
    {
        double timeElapsed = TimeElapsed();
        fMeanTime += timeElapsed;

        if (!i || timeElapsed < fMinTime)
            fMinTime = timeElapsed;
    }

    fMeanTime = (double)(fMeanTime) / iterations_;
    return fMeanTime;
}

void Bench::PrintResult() const
{
    std::cout << "Test " << fName << " Results: " << std::endl;
    std::cout << "Average of " << fMeanTime << "ms over " << fIterations << " iterations" << std::endl;
    std::cout << "Min execution time " << fMinTime << "ms" << std::endl;
}
