#include <ParameterManager.h>
#include <DoubleParameter.h>
#include <iostream>
#include <ContainerTools.hpp>
#include <Exceptions.h>

using ContainerTools::CompareKeys;
using ContainerTools::GetKeys;
using ContainerTools::HasKey;
using ContainerTools::ToString;

ParameterManager::~ParameterManager()
{
    for (std::map<std::string, FitParameter *>::iterator it = fParamPtrs.begin(); it != fParamPtrs.end(); ++it)
        delete it->second;
}

void ParameterManager::SetParameters(const ParameterDict &params_)
{
    std::set<std::string> names = GetParameterNames();
    try
    {
        for (std::set<std::string>::const_iterator it = names.begin(); it != names.end(); ++it)
            fParamPtrs.at(*it)->Set(params_.at(*it));
    }
    catch (const std::out_of_range &e_)
    {
        throw ParameterError("Parameters missing from set dict! \n " + CompareKeys(fParamPtrs, params_, "Parameter Names", "SetDict"));
    }
}

ParameterDict
ParameterManager::GetParameters() const
{
    ParameterDict returnD;
    for (std::map<std::string, FitParameter *>::const_iterator it = fParamPtrs.begin(); it != fParamPtrs.end(); ++it)
    {
        returnD[it->first] = it->second->Get();
    }
    return returnD;
}

std::set<std::string>
ParameterManager::GetParameterNames() const
{
    return GetKeys(fParamPtrs);
}

size_t
ParameterManager::GetParameterCount() const
{
    return fParamPtrs.size();
}

void ParameterManager::Add(FitParameter *p_, const std::string &n_)
{
    fParamPtrs[n_] = p_;
}

void ParameterManager::AddDouble(double &d_, const std::string &n_)
{
    fParamPtrs[n_] = new DoubleParameter(d_);
}

void ParameterManager::Clear()
{
    for (std::map<std::string, FitParameter *>::iterator it = fParamPtrs.begin(); it != fParamPtrs.end(); ++it)
        delete it->second;
    fParamPtrs.clear();
}

void ParameterManager::RenameParameter(const std::string &old_, const std::string &new_)
{
    if (!HasKey(fParamPtrs, old_))
        throw ParameterError("Can't rename " + old_ + " : parameter doesn't exist! Available: \n" + ToString(GetKeys(fParamPtrs)));

    if (HasKey(fParamPtrs, new_))
        throw ParameterError("Can't rename " + old_ + " to " + new_ + " : parameter already exists! Available: \n" + ToString(GetKeys(fParamPtrs)));
    // erase the old key and create a new one
    fParamPtrs[new_] = fParamPtrs[old_];
    fParamPtrs.erase(old_);
}

double
ParameterManager::GetParameter(const std::string &name_) const
{
    try
    {
        return fParamPtrs.at(name_)->Get();
    }
    catch (const std::out_of_range &)
    {
        throw ParameterError("GetParameter :  " + name_ + "does not exist");
    }
}

void ParameterManager::SetParameter(const std::string &name_, double value_)
{
    try
    {
        return fParamPtrs.at(name_)->Set(value_);
    }
    catch (const std::out_of_range &)
    {
        throw ParameterError("SetParameter : " + name_ + "does not exist");
    }
}
