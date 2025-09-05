#ifndef __BinnedNLLH__
#define __BinnedNLLH__
#include <TestStatistic.h>
#include <BinnedEDManager.h>
#include <BinnedEDShrinker.h>
#include <map>
#include <vector>

class DataSet;
class BinnedNLLH : public TestStatistic
{
public:
   BinnedNLLH() : fAlreadyShrunk(false) {}

   void SetPdfManager(const BinnedEDManager &);

   void SetNormalisations(const std::vector<double> &norms_);
   std::vector<double> GetNormalisations() const;

   void BinData();

   void SetDataDist(const BinnedED &);
   BinnedED GetDataDist() const;

   void SetBarlowBeeston(const bool);

   void SetBuffer(const std::string &dim_, unsigned lower_, unsigned upper_);
   std::pair<unsigned, unsigned> GetBuffer(const std::string &dim_) const;

   void AddPdf(const BinnedED &pdf, const NormFittingStatus norm_fitting_status = INDIRECT);
   void AddPdf(const BinnedED &pdf, const std::vector<std::string> &syss_, const NormFittingStatus norm_fitting_status = INDIRECT);
   void AddPdf(const BinnedED &pdf, const int &rate_, const NormFittingStatus norm_fitting_status = INDIRECT);
   void AddPdf(const BinnedED &pdf, const std::vector<std::string> &syss_, const int &rate_, const NormFittingStatus norm_fitting_status = INDIRECT);
   void AddPdfs(const std::vector<BinnedED> &pdfs, const std::vector<NormFittingStatus> *norm_fitting_statuses = nullptr);
   void AddPdfs(const std::vector<BinnedED> &pdfs, const std::vector<std::vector<std::string>> &syss_, const std::vector<NormFittingStatus> *norm_fitting_statuses = nullptr);
   void AddPdfs(const std::vector<BinnedED> &pdfs, const std::vector<int> &rates_, const std::vector<NormFittingStatus> *norm_fitting_statuses = nullptr);
   void AddPdfs(const std::vector<BinnedED> &pdfs, const std::vector<std::vector<std::string>> &syss_, const std::vector<int> &rates_, const std::vector<NormFittingStatus> *norm_fitting_statuses = nullptr);

   void SetBufferAsOverflow(bool b_); // true by default
   bool GetBufferAsOverflow() const;

   virtual void RegisterFitComponents();
   virtual void SetParameters(const ParameterDict &params_);
   virtual ParameterDict GetParameters() const;
   virtual std::set<std::string> GetParameterNames() const;
   virtual size_t GetParameterCount() const;
   virtual double Evaluate();

private:
   BinnedEDManager fPdfManager;
   BinnedEDShrinker fPdfShrinker;

   BinnedED fDataDist;
   bool fCalculatedDataDist;
   bool fAlreadyShrunk;

   std::vector<unsigned int> fGenRates;
   bool fUseBarlowBeeston = false;
};
#endif
