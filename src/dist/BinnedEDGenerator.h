#ifndef __OXSX_BINNED_ED_GENERATOR__
#define __OXSX_BINNED_ED_GENERATOR__
#include <BinnedED.h>

class BinnedEDGenerator
{
public:
   void SetPdfs(const std::vector<BinnedED> &pdfs_) { fPdfs = pdfs_; }
   std::vector<BinnedED> GetPdfs() const { return fPdfs; }

   void SetRates(const std::vector<double> rates_) { fRates = rates_; }
   std::vector<double> GetRates() const { return fRates; }

   BinnedED ExpectedRatesED(const std::string &name = "") const;
   BinnedED PoissonFluctuatedED(const std::string &name = "") const;

private:
   std::vector<BinnedED> fPdfs;
   std::vector<double> fRates;
   size_t RandomBin(size_t pdfIndex_) const;
};
#endif
