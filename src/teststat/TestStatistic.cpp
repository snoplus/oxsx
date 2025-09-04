#include <TestStatistic.h>

TestStatistic::~TestStatistic() = default;

void TestStatistic::SetSystematicManager(const SystematicManager &man_)
{
    fSystematicManager = man_;
}

SystematicManager &TestStatistic::GetSystematicManager()
{
    return fSystematicManager;
}

const SystematicManager &TestStatistic::GetSystematicManager() const
{
    return fSystematicManager;
}

void TestStatistic::SetConstraintManager(const ConstraintManager &man_)
{
    fConstraintManager = man_;
}

ConstraintManager &TestStatistic::GetConstraintManager()
{
    return fConstraintManager;
}

const ConstraintManager &TestStatistic::GetConstraintManager() const
{
    return fConstraintManager;
}

void TestStatistic::SetComponentManager(const ComponentManager &man_)
{
    fComponentManager = man_;
}

ComponentManager &TestStatistic::GetComponentManager()
{
    return fComponentManager;
}

const ComponentManager &TestStatistic::GetComponentManager() const
{
    return fComponentManager;
}

void TestStatistic::AddSystematic(Systematic *sys_)
{
    fSystematicManager.Add(sys_);
}

void TestStatistic::AddSystematic(Systematic *sys_, const std::string &group_)
{
    fSystematicManager.Add(sys_, group_);
}

void TestStatistic::AddSystematics(const std::vector<Systematic *> systematics_)
{
    for (const auto &systematic : systematics_)
    {
        AddSystematic(systematic);
    }
}

void TestStatistic::AddSystematics(const std::vector<Systematic *> sys_, const std::vector<std::string> &groups_)
{
    if (groups_.size() != sys_.size())
        throw DimensionError(Formatter() << "TestStatistic:: #sys_ != #group_");
    for (size_t i = 0; i < sys_.size(); i++)
        AddSystematic(sys_.at(i), groups_.at(i));
}

void TestStatistic::SetConstraint(const std::string &paramName_, double mean_, double sigma_)
{
    fConstraintManager.SetConstraint(paramName_, mean_, sigma_);
}

void TestStatistic::SetConstraint(const std::string &paramName_, double mean_, double sigma_lo_, double sigma_hi_)
{
    fConstraintManager.SetConstraint(paramName_, mean_, sigma_hi_, sigma_lo_);
}

void TestStatistic::SetConstraint(const std::string &paramName_1, double mean_1, double sigma_1,
                                  const std::string &paramName_2, double mean_2, double sigma_2, double correlation)
{
    fConstraintManager.SetConstraint(paramName_1, mean_1, sigma_1, paramName_2, mean_2, sigma_2, correlation);
}

void TestStatistic::SetConstraint(const std::string &paramName_1, const std::string &paramName_2, double ratiomean_, double ratiosigma_)
{
    fConstraintManager.SetConstraint(paramName_1, paramName_2, ratiomean_, ratiosigma_);
}

void TestStatistic::SetDataSet(DataSet *dataSet_)
{
    fDataSet = dataSet_;
    fCalculatedDataDist = false;
}

DataSet *
TestStatistic::GetDataSet()
{
    return fDataSet;
}

void TestStatistic::SetCalculatedDataDist(bool calculatedDataDist_)
{
    fCalculatedDataDist = calculatedDataDist_;
}

bool TestStatistic::GetCalculatedDataDist()
{
    return fCalculatedDataDist;
}

void TestStatistic::AddCut(const Cut &cut_)
{
    fCuts.AddCut(cut_);
}

void TestStatistic::SetCuts(const CutCollection &cuts_)
{
    fCuts = cuts_;
}

CutCollection TestStatistic::GetCuts()
{
    return fCuts;
}

double
TestStatistic::GetSignalCutEfficiency() const
{
    return fSignalCutEfficiency;
}

void TestStatistic::SetSignalCutEfficiency(double eff_)
{
    fSignalCutEfficiency = eff_;
}

CutLog
TestStatistic::GetSignalCutLog() const
{
    return fSignalCutLog;
}

void TestStatistic::SetSignalCutLog(const CutLog &lg_)
{
    fSignalCutLog = lg_;
}

void TestStatistic::SetDebugMode(bool mode_)
{
    fDebugMode = mode_;
    GetConstraintManager().SetDebugMode(mode_);
}

void TestStatistic::SetParameters(const ParameterDict &params_)
{
    try
    {
        fComponentManager.SetParameters(params_);
    }
    catch (const ParameterError &e_)
    {
        throw ParameterError(std::string("TestStatistic::") + e_.what());
    }
}

ParameterDict
TestStatistic::GetParameters() const
{
    return fComponentManager.GetParameters();
}

size_t
TestStatistic::GetParameterCount() const
{
    return fComponentManager.GetTotalParameterCount();
}

std::set<std::string>
TestStatistic::GetParameterNames() const
{
    return fComponentManager.GetParameterNames();
}