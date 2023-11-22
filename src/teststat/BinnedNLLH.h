#ifndef __BinnedNLLH__
#define __BinnedNLLH__
#include <TestStatistic.h>
#include <BinnedEDManager.h>
#include <SystematicManager.h>
#include <BinnedEDShrinker.h>
#include <ComponentManager.h>
#include <DataSet.h>
#include <CutCollection.h>
#include <CutLog.h>
#include <QuadraticConstraint.h>
#include <map>
#include <vector>

class DataSet;
class BinnedNLLH : public TestStatistic{
 public:
    BinnedNLLH() : fDataSet(NULL), fSignalCutEfficiency(1), fCalculatedDataDist(false), fAlreadyShrunk(false), fDebugMode(false) {}

    void   SetPdfManager(const BinnedEDManager&);
    void   SetSystematicManager(const SystematicManager&);

    void   AddPdf(const BinnedED&, const NormFittingStatus norm_fitting_status=DIRECT);
    void   AddSystematic(Systematic* sys_);
    void   AddSystematic(Systematic* sys_, const std::string& group_ );

    void   AddSystematics(const std::vector<Systematic*>);
    void   AddSystematics(const std::vector<Systematic*>, const std::vector<std::string>&);

    void   SetConstraint(const std::string& paramName_, double mean_, double sigma_);
    void   SetConstraint(const std::string& paramName_, double mean_, double sigma_lo_, double sigma_hi_);
    
    void SetNormalisations(const std::vector<double>& norms_);
    std::vector<double> GetNormalisations() const;

    void  BinData();

    void SetDataDist(const BinnedED&);
    BinnedED GetDataDist() const;

    void SetBarlowBeeston(const bool);

    void SetDataSet(DataSet*);
    DataSet* GetDataSet();

    void SetBuffer(const std::string& dim_, unsigned lower_, unsigned upper_);
    std::pair<unsigned, unsigned> GetBuffer(const std::string& dim_) const;

    void AddPdf(const BinnedED& pdf, const std::vector<std::string>& syss_, const NormFittingStatus norm_fitting_status=DIRECT);
    void AddPdf(const BinnedED& pdf, const int& rate_, const NormFittingStatus norm_fitting_status=DIRECT);
    void AddPdf(const BinnedED& pdf, const std::vector<std::string>& syss_, const int& rate_, const NormFittingStatus norm_fitting_status=DIRECT);
    void AddPdfs(const std::vector<BinnedED>& pdfs, const std::vector<NormFittingStatus>* norm_fitting_statuses=nullptr);
    void AddPdfs(const std::vector<BinnedED>& pdfs, const std::vector<std::vector<std::string> >& syss_, const std::vector<NormFittingStatus>* norm_fitting_statuses=nullptr);
    void AddPdfs(const std::vector<BinnedED>& pdfs, const std::vector<int>& rates_, const std::vector<NormFittingStatus>* norm_fitting_statuses=nullptr);
    void AddPdfs(const std::vector<BinnedED>& pdfs, const std::vector<std::vector<std::string> >& syss_, const std::vector<int>& rates_, const std::vector<NormFittingStatus>* norm_fitting_statuses=nullptr);

    void SetBufferAsOverflow(bool b_); // true by default
    bool GetBufferAsOverflow() const;

    void AddCut(const Cut&);
    void SetCuts(const CutCollection&);

    double GetSignalCutEfficiency() const;
    void   SetSignalCutEfficiency(double);

    CutLog GetSignalCutLog() const;
    void   SetSignalCutLog(const CutLog&);

    bool GetDebugMode() const { return fDebugMode; }
    void SetDebugMode(bool mode) { fDebugMode = mode; }

    // Test statistic interface
    void RegisterFitComponents(); 
    void SetParameters(const ParameterDict&);
    ParameterDict GetParameters() const;
    size_t GetParameterCount() const;
    std::set<std::string> GetParameterNames() const;
    double Evaluate();

 private:
    BinnedEDManager      fPdfManager;
    SystematicManager    fSystematicManager;
    BinnedEDShrinker     fPdfShrinker;
    DataSet*             fDataSet;
    CutCollection        fCuts;
    std::map<std::string, QuadraticConstraint> fConstraints;

    double  fSignalCutEfficiency;
    CutLog  fSignalCutLog;

    BinnedED         fDataDist;
    bool             fCalculatedDataDist;
    bool             fAlreadyShrunk;
    ComponentManager fComponentManager;

    std::vector<int> fGenRates;
    bool fUseBarlowBeeston;    
    
    bool fDebugMode;
};
#endif
