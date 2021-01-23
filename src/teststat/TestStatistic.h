/****************************************************************************/
/* Test Statistics must implement the Evaluate method for optimisation,     */
/* as well as registering any components that should be adjusted in the fit */
/*     using AddFitComponent() in the method RegisterFitComponents()        */
/****************************************************************************/
#ifndef __OXSX_TEST_STATISTIC__
#define __OXSX_TEST_STATISTIC__
#include <set>
#include <ParameterDict.h>
#include <string>
#include <stddef.h>

class TestStatistic{
 public:
    virtual ~TestStatistic() {}

    virtual double Evaluate() = 0;    
    virtual void   SetParameters(const ParameterDict& params_) = 0;
    virtual ParameterDict GetParameters() const = 0;
    virtual size_t GetParameterCount() const = 0;

    virtual std::set<std::string> GetParameterNames() const = 0;
    
    // Set up all the components for a fit
    virtual void RegisterFitComponents() = 0;
};
#endif
