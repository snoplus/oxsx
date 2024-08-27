/*****************************************************************************/
/* A class that times the exectution of a function a certain number of times */
/*****************************************************************************/
#ifndef __BENCH__
#define __BENCH__
#include <string>

class Bench
{
public:
    Bench(void (*func_)(), const std::string &name_ = " ") : fFunction(func_),
                                                             fMeanTime(-1), fName(name_), fMinTime(-1) {}

    typedef long long int64;
    typedef unsigned long long uint64;

    long double TimeRepetitions(unsigned iterations_);
    void PrintResult() const;

private:
    void (*fFunction)();
    uint64 GetTimeMs64() const;
    uint64 TimeElapsed() const;
    double fMeanTime;
    double fMinTime;
    unsigned fIterations;
    std::string fName;
};
#endif
