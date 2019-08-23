#ifndef __OXSX_PDF__
#define __OXSX_PDF__
#include <Function.h>
// Remember to implement MakeFittable!!
class PDF : public Function{
 public:
    virtual ~PDF() {}

    virtual double Integral(const std::vector<double>& mins_,
                            const std::vector<double>& maxs_) const = 0;
    virtual double Integral(const std::vector<double>& mins_,
			    const std::vector<double>& maxs_,
    			    const double& bincentre) const = 0;
    virtual std::vector<double> Sample() const = 0;
};
#endif
