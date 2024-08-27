#ifndef __OXSX_FIT_COMPONENT__
#define __OXSX_FIT_COMPONENT__
#include <vector>
#include <set>
#include <string>
#include <ParameterDict.h>

class FitComponent
{
public:
    virtual ~FitComponent() {}

    virtual void SetParameter(const std::string &name_, double value) = 0;
    virtual double GetParameter(const std::string &name_) const = 0;

    virtual void SetParameters(const ParameterDict &) = 0;
    virtual ParameterDict GetParameters() const = 0;
    virtual size_t GetParameterCount() const = 0;

    virtual std::set<std::string> GetParameterNames() const = 0;
    virtual void RenameParameter(const std::string &old_, const std::string &new_) = 0;

    virtual std::string GetName() const = 0;
    virtual void SetName(const std::string &) = 0;
};
#endif
