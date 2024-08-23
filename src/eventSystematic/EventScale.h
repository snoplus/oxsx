#ifndef __OXSX_EVENT_SCALE__
#define __OXSX_EVENT_SCALE__
#include <EventSystematic.h>
#include <FitComponent.h>
#include <vector>

class EventScale : public EventSystematic
{
public:
   EventScale(const std::string &name_) : fScale(1), fName(name_), fParamName("scale") {}
   virtual Event operator()(const Event &);

   void SetScale(double);
   double GetScale() const;

   // Fit Component Interface - just the scale parameter, set directly
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
   double fScale;
   std::string fName;
   std::string fParamName;
};
#endif
