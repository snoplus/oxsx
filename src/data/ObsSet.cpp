#include <ObsSet.h>
#include <Exceptions.h>
#include <iostream>

ObsSet::ObsSet(const std::vector<std::string> &names_)
{
    SetNames(names_);
}

ObsSet::ObsSet(const std::string &name_)
{
    SetNames(std::vector<std::string>(1, name_));
}

const std::vector<std::string> &
ObsSet::GetNames() const
{
    return fNames;
}

void ObsSet::SetNames(const std::vector<std::string> &names_)
{
    fNames = names_;
    fInverseMap.clear();
    for (size_t i = 0; i < fNames.size(); i++)
        fInverseMap[fNames.at(i)] = i;
}

size_t
ObsSet::GetIndex(const std::string &name_) const
{
    return fInverseMap.at(name_);
}

size_t
ObsSet::GetNObservables() const
{
    return fNames.size();
}

std::vector<size_t>
ObsSet::GetRelativeIndices(const ObsSet &otherRep_) const
{
    // given a longer represetation containing the obs of this rep, where are the shared indices
    // in the big rep

    if (otherRep_.GetNObservables() < GetNObservables())
        throw RepresentationError("Tried to take relatice indices in compatible representations ");

    std::vector<size_t> relativeIndices;
    relativeIndices.reserve(GetNObservables());
    for (size_t i = 0; i < fNames.size(); i++)
        relativeIndices.push_back(otherRep_.GetIndex(fNames.at(i)));

    return relativeIndices;
}

bool ObsSet::operator==(const ObsSet &other_) const
{
    return fNames == other_.fNames;
}

bool ObsSet::operator!=(const ObsSet &other_) const
{
    return !(*this == other_);
}
