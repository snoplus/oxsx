#ifndef __OXSX_MC_SAMPLER__
#define __OXSX_MC_SAMPLER__
#include <ParameterDict.h>
class MCSampler{
 public:
    virtual ~MCSampler() {}
    virtual ParameterDict Draw(const ParameterDict& currentStep_) = 0;
    virtual double CorrectAccParam(double in_) = 0;
};

#endif
