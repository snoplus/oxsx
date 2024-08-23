#ifndef __OXSX_BOX_CUT__
#define __OXSX_BOX_CUT__
#include <Event.h>
#include <Cut.h>
#include <string>

class BoxCut : public Cut
{
public:
   BoxCut(const std::string &name_, const std::string &obs_, double lower_, double upper_)
       : fName(name_), fObs(obs_), fUpperLim(upper_), fLowerLim(lower_) {}

   virtual Cut *Clone() const;

   virtual bool PassesCut(const Event &ev_) const;

   void SetLowerLimit(double lower_) { fLowerLim = lower_; }
   double GetLowerLimit() const { return fLowerLim; }

   void SetUpperLimit(double upper_) { fUpperLim = upper_; }
   double GetUpperLimit() const { return fUpperLim; }

   virtual std::string GetName() const;

private:
   std::string fName;
   std::string fObs;
   double fUpperLim;
   double fLowerLim;
};
#endif
