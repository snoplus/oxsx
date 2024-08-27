#include <EventSystematicManager.h>
#include <EventSystematic.h>
#include <Event.h>

void EventSystematicManager::Clear()
{
    fSystematics.clear();
}

void EventSystematicManager::Add(EventSystematic *sys_)
{
    fSystematics.push_back(sys_);
}

const std::vector<EventSystematic *> &
EventSystematicManager::GetSystematics() const
{
    return fSystematics;
}

size_t
EventSystematicManager::GetNSystematics() const
{
    return fSystematics.size();
}

Event EventSystematicManager::ApplySystematics(const Event &event_) const
{
    if (!fSystematics.size())
        return event_;

    Event modified = event_;
    for (size_t i = 0; i < fSystematics.size(); i++)
        modified = fSystematics.at(i)->operator()(modified);
    return modified;
}
