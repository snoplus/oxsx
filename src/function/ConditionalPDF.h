#ifndef __OXSX_CONDITIONAL_PDF__
#define __OXSX_CONDITIONAL_PDF__
#include <FitComponent.h>
#include <vector>

class ConditionalPDF : public FitComponent{
 public:
  virtual ~ConditionalPDF() {};
  virtual ConditionalPDF* Clone() const = 0;
  virtual double ConditionalProbability(const std::vector<double>& x_, 
    const std::vector<double>& x2_) = 0;
  
  virtual double Integral(const std::vector<double>& mins_,
    const std::vector<double>& maxs_,
    const std::vector<double>& x2_) const = 0;

  virtual std::vector<double> Sample(const std::vector<double>& x2_) const = 0;
};
#endif
