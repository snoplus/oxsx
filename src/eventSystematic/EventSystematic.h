#ifndef __OXSX_EVENT_SYSTEMATIC__
#define __OXSX_EVENT_SYSTEMATIC__
#include <ObsSet.h>
#include <FitComponent.h>
#include <vector>

class Event;
class EventSystematic : public FitComponent{
 public:
    virtual Event operator()(const Event&) = 0;

    void   SetInObservables(const std::vector<std::string>&);
    const std::vector<std::string>& GetInObservables() const;

    void   SetOutObservables(const std::vector<std::string>&);
    const std::vector<std::string>& GetOutObservables() const;

    
 protected:   
    ObsSet  fInObservables;
    ObsSet  fOutObservables;
};
#endif
