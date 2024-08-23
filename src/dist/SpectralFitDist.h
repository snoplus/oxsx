#ifndef __OXSX_SPECTRAL_FIT_PDF__
#define __OXSX_SPECTRAL_FIT_PDF__
#include <BinnedED.h>
#include <FitComponent.h>

class SpectralFitDist : public BinnedED, public FitComponent
{
public:
   SpectralFitDist(const std::string &name_,
                   const AxisCollection &axes_);

   // Make this fittable with each bin content adjustable
   void SetParameter(const std::string &name_, double value);
   double GetParameter(const std::string &name_) const;

   void SetParameters(const ParameterDict &);
   ParameterDict GetParameters() const;
   size_t GetParameterCount() const;

   std::set<std::string> GetParameterNames() const;
   void RenameParameter(const std::string &old_, const std::string &new_);

   std::string GetName() const;
   void SetName(const std::string &name_);
   // it already has get/set name from BinnedED
private:
   std::vector<std::string> fBinNames;
   std::string fName;
};
#endif
