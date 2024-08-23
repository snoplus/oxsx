#include <GaussianFitter.h>
#include <Gaussian.h>
#include <iostream>
#include <sstream>
#include <Formatter.hpp>
#include <ContainerTools.hpp>
#include <Exceptions.h>
#include <Formatter.hpp>
#include <algorithm>

using ContainerTools::ToString;

GaussianFitter::GaussianFitter(Gaussian *gaus, const size_t &nDims_)
{
    fOrignalFunc = gaus;
    std::stringstream ss;
    for (size_t i = 0; i < nDims_; ++i)
    {
        ss << "means" << "_" << i;
        fMeansNames.push_back(ss.str());
        ss.str("");
        ss << "stddevs" << "_" << i;
        fStdDevsNames.push_back(ss.str());
        ss.str("");
    }
}

GaussianFitter::GaussianFitter(Gaussian *gaus, const std::vector<std::string> &meanNames_, const std::vector<std::string> &stdDevNames_)
{
    if (meanNames_.size() != stdDevNames_.size())
        throw OXSXException(Formatter() << "GaussianFitter:: #meanName != #stdDevNames");
    fOrignalFunc = gaus;
    std::stringstream ss;
    for (const auto &mean_name : meanNames_)
        fMeansNames.push_back(mean_name);
    for (const auto &stdDev_name : stdDevNames_)
        fStdDevsNames.push_back(stdDev_name);
}

void GaussianFitter::RenameParameter(const std::string &old_, const std::string &new_)
{
    std::vector<std::string>::iterator it;
    if (find(fMeansNames.begin(), fMeansNames.end(), old_) == fMeansNames.end() && find(fStdDevsNames.begin(), fStdDevsNames.end(), old_) == fStdDevsNames.end())
        throw NotFoundError(Formatter() << "GaussianFitter:: When attempting to renaming the parameter " << old_ << ", it wasn't found. Available names: " << ToString(GetParameterNames()));

    it = find(fMeansNames.begin(), fMeansNames.end(), old_);
    while (it != fMeansNames.end())
    {
        *it = new_;
        it = find(it++, fMeansNames.end(), old_);
    }

    it = find(fStdDevsNames.begin(), fStdDevsNames.end(), old_);
    while (it != fStdDevsNames.end())
    {
        *it = new_;
        it = find(it++, fStdDevsNames.end(), old_);
    }
}

std::vector<std::string>
GaussianFitter::GetMeanNames() const
{
    return fMeansNames;
}

std::vector<std::string>
GaussianFitter::GetStdDevNames() const
{
    return fStdDevsNames;
}

void GaussianFitter::SetParameter(const std::string &name_, double value_)
{
    std::vector<std::string>::iterator it;
    it = find(fMeansNames.begin(), fMeansNames.end(), name_);
    while (it != fMeansNames.end())
    {
        fOrignalFunc->SetMean(it - fMeansNames.begin(), value_);
        it = find(++it, fMeansNames.end(), name_);
    }
    it = find(fStdDevsNames.begin(), fStdDevsNames.end(), name_);
    while (it != fStdDevsNames.end())
    {
        fOrignalFunc->SetStDev(it - fStdDevsNames.begin(), value_);
        it = find(++it, fStdDevsNames.end(), name_);
    }
}

double
GaussianFitter::GetParameter(const std::string &name_) const
{
    // BL: If n parameters have the same name (either across both means
    // and stddevs or not) the value will be the same. So the value of the
    // first instance is sufficient.

    std::vector<std::string>::const_iterator it;
    it = find(fMeansNames.begin(), fMeansNames.end(), name_);
    if (it == fMeansNames.end())
    {
        it = find(fStdDevsNames.begin(), fStdDevsNames.end(), name_);
        if (it == fStdDevsNames.end())
            throw NotFoundError(Formatter() << "GaussianFitter:: Parameter : " << name_ << " was not known to the GaussianFitter. Available names: " << ToString(GetParameterNames()));
        return fOrignalFunc->GetStDev(it - fStdDevsNames.end());
    }
    return fOrignalFunc->GetMean(it - fStdDevsNames.end());
}

void GaussianFitter::SetParameters(const ParameterDict &ps_)
{
    for (ParameterDict::const_iterator i = ps_.begin(); i != ps_.end(); ++i)
    {
        SetParameter(i->first, i->second);
    }
}

ParameterDict
GaussianFitter::GetParameters() const
{
    std::vector<double> means = fOrignalFunc->GetMeans();
    std::vector<double> stddevs = fOrignalFunc->GetStdDevs();
    std::vector<double> values;

    values.reserve(means.size() + stddevs.size()); // preallocate memory
    values.insert(values.end(), means.begin(), means.end());
    values.insert(values.end(), stddevs.begin(), stddevs.end());

    std::vector<std::string> names;
    names.reserve(fMeansNames.size() + fStdDevsNames.size()); // preallocate memory
    names.insert(names.end(), fMeansNames.begin(), fMeansNames.end());
    names.insert(names.end(), fStdDevsNames.begin(), fStdDevsNames.end());

    return ContainerTools::CreateMap(names, values);
}

size_t
GaussianFitter::GetParameterCount() const
{
    return fMeansNames.size() + fStdDevsNames.size();
}

std::set<std::string>
GaussianFitter::GetParameterNames() const
{
    std::set<std::string> names;
    for (size_t i = 0; i < fMeansNames.size(); ++i)
    {
        names.insert(fMeansNames.at(i));
        names.insert(fStdDevsNames.at(i));
    }
    return names;
}
