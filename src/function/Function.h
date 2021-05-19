#ifndef __OXSX_FUNCTION__
#define __OXSX_FUNCTION__
#include <FitComponent.h>
#include <vector>

class Function : public FitComponent{
 public:
    Function() {}
    virtual ~Function() {}
    virtual Function* Clone() const = 0;

    virtual double operator()(const std::vector<double>& vals_) const = 0;
    virtual size_t GetNDims() const = 0;
};
#endif
