/****************************************/
/* Abstract base class for an optimiser */
/****************************************/
#ifndef __OXSX_OPTIMISER__
#define __OXSX_OPTIMISER__
#include <vector>

class TestStatistic;
class FitResult;

class Optimiser
{
public:
    virtual ~Optimiser() {}
    virtual const FitResult &Optimise(TestStatistic *) = 0;
};
#endif
