#ifndef __OXSX_DOUBLE_POINTER__
#define __OXSX_DOUBLE_POINTER__
#include <FitParameter.h>

class DoubleParameter : public FitParameter
{
public:
   DoubleParameter(double &ptr_) : fPtr(ptr_) {}
   double Get() const { return fPtr; }
   void Set(double d_) { fPtr = d_; }

private:
   double &fPtr;
};
#endif
