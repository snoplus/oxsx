#include <SpectralFitDist.h>
#include <Exceptions.h>
#include <Formatter.hpp>
#include <ContainerTools.hpp>
#include <sstream>
#include <algorithm>
#include <iostream>

using ContainerTools::CreateMap;
using ContainerTools::ToString;

// JD: this will probably be too slow one day - there is overkill checking here.
// you can just set all of the bin contents in one go, but noone needs this code yet

SpectralFitDist::SpectralFitDist(const std::string &name_,
                                 const AxisCollection &axes_) : BinnedED(name_, axes_)
{
    fBinNames.reserve(axes_.GetNBins());
    for (size_t i = 0; i < axes_.GetNBins(); i++)
        fBinNames.push_back(std::string(Formatter() << name_ << "_bin_" << i));

    fName = name_;
}

void SpectralFitDist::SetParameter(const std::string &name_, double value)
{
    std::vector<std::string>::iterator it = std::find(fBinNames.begin(), fBinNames.end(), name_);
    if (it == fBinNames.end())
        throw ParameterError("Can't set " + name_ + ", parameters are called " + ToString(fBinNames));

    *it = value;
}

double
SpectralFitDist::GetParameter(const std::string &name_) const
{
    std::vector<std::string>::const_iterator it = std::find(fBinNames.begin(), fBinNames.end(), name_);
    if (it == fBinNames.end())
        throw ParameterError("Can't get " + name_ + ", parameters are called: \n " + ToString(fBinNames));
    return fHistogram.GetBinContent(it - fBinNames.begin());
}

void SpectralFitDist::SetParameters(const ParameterDict &pd_)
{
    for (size_t i = 0; i < fBinNames.size(); i++)
        try
        {
            fHistogram.SetBinContent(i, pd_.at(fBinNames.at(i)));
        }
        catch (const std::out_of_range &e_)
        {
            throw ParameterError("Set Parameters: Couldn't set bin content " + fBinNames.at(i));
        }
}

ParameterDict
SpectralFitDist::GetParameters() const
{
    return CreateMap(fBinNames, fHistogram.GetBinContents());
}

size_t
SpectralFitDist::GetParameterCount() const
{
    return fHistogram.GetNBins();
}

std::set<std::string>
SpectralFitDist::GetParameterNames() const
{
    return std::set<std::string>(fBinNames.begin(), fBinNames.end());
}

void SpectralFitDist::RenameParameter(const std::string &old_, const std::string &new_)
{
    std::vector<std::string>::iterator it = std::find(fBinNames.begin(), fBinNames.end(), old_);
    if (it == fBinNames.end())
        throw ParameterError("Can't rename " + old_ + " it doesn't exist, parameters are called: \n " + ToString(fBinNames));
    *it = new_;
}

void SpectralFitDist::SetName(const std::string &name_)
{
    fName = name_;
}

std::string
SpectralFitDist::GetName() const
{
    return fName;
}
