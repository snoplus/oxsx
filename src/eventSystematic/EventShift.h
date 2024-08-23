#ifndef __OXSX_EVENT_SHIFT__
#define __OXSX_EVENT_SHIFT__
#include <EventSystematic.h>
#include <FitComponent.h>
#include <vector>

class EventShift : public EventSystematic
{
public:
   EventShift(const std::string &name_) : fShift(0), fName(name_), fParamName("shift") {}
   Event operator()(const Event &);

   void SetShift(double);
   double GetShift() const;

   // FitComponent Interface - just the shift, set directly
   void SetParameter(const std::string &name_, double value);
   double GetParameter(const std::string &name_) const;

   void SetParameters(const ParameterDict &);
   ParameterDict GetParameters() const;
   size_t GetParameterCount() const;

   std::set<std::string> GetParameterNames() const;
   void RenameParameter(const std::string &old_, const std::string &new_);

   std::string GetName() const;
   void SetName(const std::string &);

private:
   double fShift;
   std::string fName;
   std::string fParamName;
};
#endif
