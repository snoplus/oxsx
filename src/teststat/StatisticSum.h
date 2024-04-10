/* 
Note that this object is formed without copying the statistics to add
This means is should exist in the same scope as those summed.
*/

#ifndef __OXSX_SUMMED_STATISTIC__
#define __OXSX_SUMMED_STATISTIC__
#include <TestStatistic.h>
#include <vector>

class StatisticSum : public TestStatistic{
 public:
    StatisticSum() {}
    StatisticSum(TestStatistic&, TestStatistic&);
    void AddStat(TestStatistic&);

    // Test Statistic Interface
    virtual double Evaluate();
    virtual void   SetParameters(const ParameterDict& params_);
    virtual ParameterDict GetParameters() const;
    virtual size_t GetParameterCount() const;

    virtual std::set<std::string> GetParameterNames() const;
    
    // Set up all the components for a fit
    virtual void RegisterFitComponents();
    
 private:
    std::vector<TestStatistic*> fStats;
};

StatisticSum operator + (TestStatistic&, TestStatistic&);
StatisticSum operator + (StatisticSum&, TestStatistic&);
StatisticSum Sum(std::vector<TestStatistic*>&);
#endif
