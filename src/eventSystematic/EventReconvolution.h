#ifndef __OXSX_EVENT_RECONVOLUTION__
#define __OXSX_EVENT_RECONVOLUTION__
#include <EventSystematic.h>
#include <FitComponent.h>

class Event;
class EventReconvolution : public EventSystematic
{
public:
   EventReconvolution(const std::string &name_) : fCorrection(1), fName(name_), fParamName("correction") {}
   Event operator()(const Event &);

   double GetCorrection() const;
   void SetCorrection(double);

   // Fit Component Interface - just the correction, set directly
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
   double fCorrection;
   std::string fName;
   std::string fParamName;
};
#endif
