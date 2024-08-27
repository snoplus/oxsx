#ifndef __OXSX_FIT_PARAMETER__
#define __OXSX_FIT_PARAMETER__

class FitParameter
{
public:
    virtual ~FitParameter() {}
    virtual double Get() const = 0;
    virtual void Set(double) = 0;
};

#endif
