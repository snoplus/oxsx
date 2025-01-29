/***************************************************************************************************/
/* A collection of BinnedPDFs with normalisations for calculating the probability of a given obs   */
/* given some rates and some pdfs                                                                  */
/***************************************************************************************************/
#ifndef __OXSX_BINNED_ED_MANAGER__
#define __OXSX_BINNED_ED_MANAGER__
#include <BinnedED.h>
#include <FitComponent.h>
#include <ParameterManager.h>
#include <vector>

class Event;
class SystematicManager;
class BinnedEDShrinker;

enum NormFittingStatus
{
   FALSE = 0,
   DIRECT = 1,
   INDIRECT = 2
};
class BinnedEDManager : public FitComponent
{
public:
   BinnedEDManager() : fAllNormsDirFittable(true), fNPdfs(0), fName("norms") {}

   void AddPdf(const BinnedED &, const NormFittingStatus norm_fitting_status = INDIRECT);
   void AddPdfs(const std::vector<BinnedED> &,
                const std::vector<NormFittingStatus> *norm_fitting_statuses = nullptr);

   double Probability(const Event &);
   double BinProbability(size_t);

   const std::vector<double> &GetNormalisations() const;
   void SetNormalisations(const std::vector<double> &normalisations_);

   void ApplySystematics(const SystematicManager &sysMan_);
   void ApplyShrink(const BinnedEDShrinker &);

   const BinnedED &GetOriginalPdf(size_t index_) const;
   unsigned GetNPdfs() const;
   size_t GetNDims() const;

   void AssertDimensions(const std::vector<std::string> &observables);
   void ReassertNorms(bool calcing_binprob = false);

   // Make a fittable component - i.e. rescale the binned pdfs inside to fit
   void SetParameter(const std::string &name_, double value);
   double GetParameter(const std::string &name_) const;

   void SetParameters(const ParameterDict &);
   ParameterDict GetParameters() const;
   size_t GetParameterCount() const;

   std::set<std::string> GetParameterNames() const;
   void RenameParameter(const std::string &old_, const std::string &new_);

   std::string GetName() const;
   void SetName(const std::string &);

private:
   ParameterManager fParameterManager;
   std::vector<BinnedED> fOriginalPdfs;
   std::vector<BinnedED> fWorkingPdfs;
   std::vector<double> fNormalisations;
   std::vector<NormFittingStatus> fAllowNormsFittable;
   std::vector<double> fFittableNorms;
   bool fAllNormsDirFittable;
   int fNPdfs;
   size_t fNDims;

   std::string fName; // component name

   void RegisterParameters();
};
#endif
