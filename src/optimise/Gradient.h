#ifndef __OXSX_GRADIENT__
#define __OXSX_GRADIENT__
#include <ParameterDict.h>

template <typename Stat>
class Gradient
{
public:
   Gradient(Stat &st_, double epsilon_) : fStat(st_), fEpsilon(epsilon_) {}
   ParameterDict operator()(const ParameterDict &x0_);

private:
   Stat &fStat;
   double fEpsilon;
};

#include <Gradient.hpp>
#endif
