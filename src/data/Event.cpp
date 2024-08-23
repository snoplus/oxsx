#include <Event.h>
#include <ObsSet.h>
#include <Exceptions.h>
#include <ContainerTools.hpp>
#include <algorithm>

using ContainerTools::CreateMap;

Event::Event(const std::vector<double> &obs_)
{
    fObservations = obs_;
    fNObservables = obs_.size();
    fObsNames = NULL;
}

std::vector<double>
Event::GetData() const
{
    return fObservations;
}

std::map<std::string, double>
Event::GetDataMap() const
{
    return CreateMap(GetObservableNames(), fObservations);
}

double
Event::GetDatum(size_t index_) const
{
    try
    {
        return fObservations.at(index_);
    }

    catch (const std::out_of_range &e_)
    {
        throw NotFoundError(Formatter() << "Event::Attempted access on non-existent observable " << index_);
    }
}

void Event::SetDatum(size_t index_, double val_)
{
    try
    {
        fObservations[index_] = val_;
    }

    catch (const std::out_of_range &e_)
    {
        throw NotFoundError(Formatter() << "Event::Attempted access on non-existent observable " << index_);
    }
}

void Event::SetDatum(const std::string &name_, double val_)
{

    // work out which index that is
    std::vector<std::string>::const_iterator it = std::find(fObsNames->begin(), fObsNames->end(), name_);
    if (it == fObsNames->end())
        throw NotFoundError("Event::GetDatum  event doesn't contain data for observable " + name_ + "!");

    SetDatum(it - fObsNames->begin(), val_);
}

double
Event::GetDatum(const std::string &name_) const
{

    // work out which index that is
    std::vector<std::string>::const_iterator it = std::find(fObsNames->begin(), fObsNames->end(), name_);
    if (it == fObsNames->end())
        throw NotFoundError("Event::GetDatum  event doesn't contain data for observable " + name_ + "!");

    return GetDatum(it - fObsNames->begin());
}

std::vector<double>
Event::ToObsSet(const ObsSet &rep_) const
{
    std::vector<double> vals;
    vals.reserve(rep_.GetNObservables());

    const std::vector<std::string> &obsNames = rep_.GetNames();
    for (size_t i = 0; i < obsNames.size(); i++)
        vals.push_back(GetDatum(obsNames.at(i)));

    return vals;
}

void Event::SetObservableNames(std::vector<std::string> const *p_)
{
    fObsNames = p_;
}

const std::vector<std::string> &
Event::GetObservableNames() const
{
    if (!fObsNames)
        throw NULLPointerAccessError("Event::GetObservableNames -> names have not been set!");
    return *fObsNames;
}
