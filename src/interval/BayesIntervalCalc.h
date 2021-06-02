#ifndef __OXSX_BAYESINTERVALCALC__
#define __OXSX_BAYESINTERVALCALC__
#include <Histogram.h>

class BayesIntervalCalc{
 public:
    static double UpperBound(Histogram lh_, double cl_); // needs to be a 1D histogram
    static double UpperBound(double expectedCounts_, int observedCounts_, double cl_, double tolerance_);
    static double OneSidedUpperInterval(double expectedCounts_, int observedCounts_, double upperEdge_);
};
#endif
