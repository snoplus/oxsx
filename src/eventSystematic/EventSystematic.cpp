#include <EventSystematic.h>

void EventSystematic::SetInObservables(const std::vector<std::string> &names_)
{
    fInObservables = ObsSet(names_);
}

const std::vector<std::string> &
EventSystematic::GetInObservables() const
{
    return fInObservables.GetNames();
}

void EventSystematic::SetOutObservables(const std::vector<std::string> &names_)
{
    fOutObservables = ObsSet(names_);
}

const std::vector<std::string> &
EventSystematic::GetOutObservables() const
{
    return fOutObservables.GetNames();
}
