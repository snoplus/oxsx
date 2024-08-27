#ifndef __EVENT_SYSTEMATIC_MANAGER__
#define __EVENT_SYSTEMATIC_MANAGER__
#include <vector>
#include <stddef.h>

class EventSystematic;
class Event;
class EventSystematicManager
{
public:
   Event ApplySystematics(const Event &) const;

   void Clear();
   void Add(EventSystematic *);
   const std::vector<EventSystematic *> &GetSystematics() const;
   size_t GetNSystematics() const;

private:
   std::vector<EventSystematic *> fSystematics;
};

#endif
