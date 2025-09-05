/****************************************************************************/
/* Test Statistics must implement the Evaluate method for optimisation,     */
/* as well as registering any components that should be adjusted in the fit */
/*     using AddFitComponent() in the method RegisterFitComponents()        */
/****************************************************************************/
#ifndef __OXSX_TEST_STATISTIC__
#define __OXSX_TEST_STATISTIC__
#include <set>
#include <ParameterDict.h>
#include <ConstraintManager.h>
#include <SystematicManager.h>
#include <ComponentManager.h>
#include <DataSet.h>
#include <CutCollection.h>
#include <CutLog.h>
#include <ConstraintManager.h>
#include <string>
#include <stddef.h>

class TestStatistic
{
public:
    TestStatistic() : fCalculatedDataDist(false), fDataSet(NULL), fSignalCutEfficiency(1), fDebugMode(false) {}

    virtual ~TestStatistic();

    void SetSystematicManager(const SystematicManager &);
    SystematicManager &GetSystematicManager();
    const SystematicManager &GetSystematicManager() const;

    void SetConstraintManager(const ConstraintManager &);
    ConstraintManager &GetConstraintManager();
    const ConstraintManager &GetConstraintManager() const;

    void SetComponentManager(const ComponentManager &);
    ComponentManager &GetComponentManager();
    const ComponentManager &GetComponentManager() const;

    void AddSystematic(Systematic *sys_);
    void AddSystematic(Systematic *sys_, const std::string &group_);
    void AddSystematics(const std::vector<Systematic *>);
    void AddSystematics(const std::vector<Systematic *>, const std::vector<std::string> &);

    void SetConstraint(const std::string &paramName_, double mean_, double sigma_);
    void SetConstraint(const std::string &paramName_, double mean_, double sigma_lo_, double sigma_hi_);
    void SetConstraint(const std::string &paramName_1, double mean_1, double sigma_1, const std::string &paramName_2, double mean_2, double sigma_2, double correlation);
    void SetConstraint(const std::string &paramName_1, const std::string &paramName_2, double ratiomean_, double ratiosigma_);

    void SetDataSet(DataSet *);
    DataSet *GetDataSet();

    void SetCalculatedDataDist(bool);
    bool GetCalculatedDataDist();

    void AddCut(const Cut &);
    void SetCuts(const CutCollection &);
    CutCollection GetCuts();

    double GetSignalCutEfficiency() const;
    void SetSignalCutEfficiency(double);

    CutLog GetSignalCutLog() const;
    void SetSignalCutLog(const CutLog &);

    bool GetDebugMode() const { return fDebugMode; }
    void SetDebugMode(bool mode);

    virtual double Evaluate() = 0;
    virtual void SetParameters(const ParameterDict &params_);
    virtual ParameterDict GetParameters() const;
    virtual std::set<std::string> GetParameterNames() const;
    virtual size_t GetParameterCount() const;

    // Set up all the components for a fit
    virtual void RegisterFitComponents() = 0;

private:
    SystematicManager fSystematicManager;
    ConstraintManager fConstraintManager;
    ComponentManager fComponentManager;

    bool fCalculatedDataDist;
    DataSet *fDataSet;

    CutCollection fCuts;
    double fSignalCutEfficiency;
    CutLog fSignalCutLog;

    bool fDebugMode;
};
#endif
