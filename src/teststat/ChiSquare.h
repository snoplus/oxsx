/********************************************************************/
/* A routine to calculate a binned chi squared test.                */
/* NOTE: this can only make sense if all pdfs have the same binning */
/********************************************************************/
#ifndef __OXSX_CHI_SQUARED__
#define __OXSX_CHI_SQUARED__
#include <SystematicManager.h>
#include <BinnedEDManager.h>
#include <TestStatistic.h>
#include <ComponentManager.h>

// FIXME::Enforce equal binning on the pdfs
class DataSet;
class ChiSquare : public TestStatistic
{
public:
   ChiSquare() : fCalculatedDataDist(false) {}

   void SetDataSet(DataSet *d);
   DataSet *GetDataSet();

   // Fit Component interface
   double Evaluate();
   void RegisterFitComponents();
   size_t GetParameterCount() const;
   void SetParameters(const ParameterDict &params_);
   ParameterDict GetParameters() const;
   std::set<std::string> GetParameterNames() const;

private:
   bool fCalculatedDataDist;
   BinnedED fDataDist;
   void BinData();

   DataSet *fDataSet;
   BinnedEDManager fPdfManager;
   SystematicManager fSystematicManager;
   ComponentManager fComponentManager;
};
#endif
