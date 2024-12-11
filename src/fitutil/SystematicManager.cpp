#include <SystematicManager.h>
#include <Exceptions.h>
#include <Formatter.hpp>
#include <ContainerTools.hpp>

using ContainerTools::GetKeys;
using ContainerTools::GetValues;

void SystematicManager::Construct()
{
    // Don't do anything if there are no systematics
    if (!CountNSystematics())
        return;

    // Construct the response matrices.
    for (std::map<std::string, std::vector<Systematic *>>::const_iterator group = fGroups.begin(); group != fGroups.end(); ++group)
    {
        // If "" is empty get on with it.
        if (group->first == "" && group->second.size() == 0)
            continue;
        for (size_t i = 0; i < group->second.size(); i++)
            fGroups[group->first].at(i)->Construct();
    }

    // This loop should construct all fGroups.
    for (std::map<std::string, std::vector<Systematic *>>::const_iterator group = fGroups.begin(); group != fGroups.end(); ++group)
    {
        // If "" is empty get on with it.
        if (group->first == "" && group->second.size() == 0)
            continue;

        SparseMatrix resp = fGroups[group->first].at(0)->GetResponse();
        for (size_t i = 1; i < group->second.size(); i++)
        {
            SparseMatrix holder = fGroups[group->first].at(i)->GetResponse();
            resp = holder * resp;
        }
        fTotalReponses[group->first] = resp;
    }
}

const SparseMatrix &
SystematicManager::GetTotalResponse(const std::string &groupName_) const
{
    // should you construct?
    try
    {
        return fTotalReponses.at(groupName_);
    }
    catch (const std::out_of_range &e_)
    {
        throw NotFoundError(Formatter() << "SystematicManager :: name : " << groupName_ << " not found in systematic map");
    }
}

void SystematicManager::Add(Systematic *sys_, const std::string &groupName_)
{
    if (groupName_ == "")
    {
        std::cout << "Warning::Adding systematic with empty group name. Will apply to all events." << std::endl;
        Add(sys_);
        return;
    }
    fGroups[groupName_].push_back(sys_);
    fNGroups = fGroups.size();
}

void SystematicManager::Add(Systematic *sys_)
{
    checkAllOtherGroups(sys_);
    fGroups[""].push_back(sys_);
    fNGroups = fGroups.size();
}

void SystematicManager::checkAllOtherGroups(const Systematic *syss_)
{
    std::set<std::string> allname;
    allname.insert(syss_->GetName());
    for (size_t j = 0; j < fGroups.at("").size(); ++j)
    {
        if (allname.find(fGroups.at("").at(j)->GetName()) == allname.end())
            allname.insert(fGroups.at("").at(j)->GetName());
        else
            throw NotFoundError(Formatter() << "SystematicManager :: Systematic: " << fGroups.at("").at(j)->GetName() << " is in the global systematic group more than once.");
    }
    for (std::map<std::string, std::vector<Systematic *>>::const_iterator i = fGroups.begin(); i != fGroups.end(); ++i)
    {
        if (i->first == "")
            continue;
        for (std::vector<Systematic *>::const_iterator sys = i->second.begin(); sys != i->second.end(); ++sys)
        {
            if (allname.find((*sys)->GetName()) == allname.end())
            {
                allname.insert((*sys)->GetName());
            }
            else
                throw NotFoundError(Formatter() << "SystematicManager :: Systematic: " << (*sys)->GetName() << " is in the another group and therefore can't be added to the global systematics.");
        }
    }
}

void SystematicManager::UniqueSystematics(const std::vector<std::string> &syss_)
{
    std::set<std::string> allname;
    for (size_t i = 0; i < syss_.size(); ++i)
    {
        if (fGroups.find(syss_.at(i)) == fGroups.end())
            throw NotFoundError(
                Formatter() << "SystematicManager:: Systematic group " << syss_.at(i) << "  is not known to the SystematicManager.");
        std::vector<Systematic *> group = fGroups[syss_.at(i)];
        for (size_t j = 0; j < group.size(); ++j)
        {
            if (allname.find(group.at(j)->GetName()) == allname.end())
                allname.insert(group.at(j)->GetName());
            else
                throw NotFoundError(Formatter() << "SystematicManager :: Systematic: " << group.at(i)->GetName() << " is in more than one group.");
        }
    }
    // Check the "" group only after everything else has been checked and
    // allnames has been populated.
    if (syss_ != std::vector<std::string>(1, ""))
    {
        std::vector<Systematic *> group = fGroups[""];
        for (size_t j = 0; j < group.size(); ++j)
        {
            if (allname.find(group.at(j)->GetName()) == allname.end())
                allname.insert(group.at(j)->GetName());
            else
                throw NotFoundError(Formatter() << "SystematicManager :: Systematic: " << group.at(j)->GetName() << " is in the global (\"\") group, and therefore would be applied twice.");
        }
    }
}

void SystematicManager::AddDist(const BinnedED &pdf, const std::vector<std::string> &syss_)
{
    // Check whether all systematics in syss_ are unique.
    UniqueSystematics(syss_);
    fEDGroups[pdf.GetName()] = syss_;
}
void SystematicManager::AddDist(const BinnedED &pdf, const std::string &syss_)
{
    AddDist(pdf, std::vector<std::string>(1, syss_));
}

void SystematicManager::DistortEDs(const std::vector<BinnedED> &OrignalEDs_,
                                   std::vector<BinnedED> &WorkingEDs_,
                                   std::vector<double> *norms) const
{
    for (size_t j = 0; j < WorkingEDs_.size(); j++)
    {
        const std::string name = OrignalEDs_.at(j).GetName();

        if (fEDGroups.find(name) == fEDGroups.end())
            continue;

        // Copy EDs because WorkingED_s have to have the same binning as OrignalEDs_s.
        WorkingEDs_[j] = OrignalEDs_.at(j);

        // If we haven't explicitly added the "" group to this dist's fEDGroups, apply that group first
        if (std::find(fEDGroups.at(name).begin(), fEDGroups.at(name).end(), "") == fEDGroups.at(name).end() && fGroups.at("").size())
        {
            WorkingEDs_[j].SetBinContents(GetTotalResponse("").operator()(WorkingEDs_.at(j).GetBinContents()));
        }

        // Apply everything else.
        for (size_t i = 0; i < fEDGroups.at(name).size(); ++i)
        {
            // Check that the group has systematics in it.
            std::string groupName = fEDGroups.at(name).at(i);
            // If "" is empty the do nothing.
            if (groupName == "" && fGroups.at("").size() == 0)
                continue;
            if (!fGroups.at(groupName).size())
                throw LogicError(Formatter() << "SystematicManager:: ED " << name
                                             << " has a systematic group of zero size acting on it");

            // Now apply this group to the WorkingEDs
            WorkingEDs_[j].SetBinContents(GetTotalResponse(groupName).operator()(WorkingEDs_.at(j).GetBinContents()));
        }
        if (norms != nullptr)
        {
            norms->operator[](j) = WorkingEDs_.at(j).Integral();
        }
    }
}

// Getters and Setters

size_t
SystematicManager::GetNSystematicsInGroup(const std::string &name_) const
{
    try
    {
        return fGroups.find(name_)->second.size();
    }
    catch (const std::out_of_range &e_)
    {
        throw NotFoundError(Formatter() << "SystematicManager :: name : " << name_ << " not known to the SystematicManager.");
    }
}

const std::set<std::string>
SystematicManager::GetGroupNames() const
{
    return GetKeys(fGroups);
}

const std::vector<std::string>
SystematicManager::GetSystematicsNamesInGroup(const std::string &name) const
{
    try
    {
        std::vector<Systematic *> sys = fGroups.at(name);
        std::vector<std::string> names;
        for (size_t i = 0; i < sys.size(); ++i)
        {
            names.push_back(sys.at(i)->GetName());
        }
        return names;
    }
    catch (const std::out_of_range &e_)
    {
        throw NotFoundError(Formatter() << "SystematicManager :: name : " << name << " not known to the SystematicManager.");
    }
}

const std::vector<Systematic *> &
SystematicManager::GetSystematicsInGroup(const std::string &name) const
{
    try
    {
        return fGroups.at(name);
    }
    catch (const std::out_of_range &e_)
    {
        throw NotFoundError(Formatter() << "SystematicManager :: name : " << name << " not known to the SystematicManager.");
    }
}

const std::vector<Systematic *> &
SystematicManager::GetGlobalSystematics() const
{
    return fGroups.at("");
}

const std::vector<std::string>
SystematicManager::GetGroup(const std::string &name) const
{
    // Check group exist.
    if (fGroups.find(name) == fGroups.end())
        throw NotFoundError(Formatter() << "SystematicManager:: Group name " << name << " not known to the SystematicManager.");

    std::vector<std::string> names;
    for (std::vector<Systematic *>::const_iterator sys = fGroups.at(name).begin(); sys != fGroups.at(name).end(); ++sys)
        names.push_back((*sys)->GetName());
    return names;
}

const std::map<std::string, std::vector<Systematic *>> &
SystematicManager::GetSystematicsGroup() const
{
    return fGroups;
}

size_t
SystematicManager::CountNSystematics() const
{
    size_t NSystematics = 0;
    for (std::map<std::string, std::vector<Systematic *>>::const_iterator group = fGroups.begin(); group != fGroups.end(); ++group)
        NSystematics += group->second.size();
    return NSystematics;
}

size_t
SystematicManager::GetNSystematics() const
{
    return CountNSystematics();
}

size_t
SystematicManager::GetNGroups() const
{
    return fNGroups;
}
