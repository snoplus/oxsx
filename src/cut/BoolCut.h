#ifndef __OXSX_BOOL_CUT__
#define __OXSX_BOOL_CUT__
#include <Cut.h>
#include <stddef.h>
#include <string>

class BoolCut : public Cut
{
public:
   BoolCut(const std::string &name_, const std::string &obs_, double value_) : fName(name_), fObs(obs_), fVal(value_) {}
   virtual bool PassesCut(const Event &ev_) const;
   virtual Cut *Clone() const;
   virtual std::string GetName() const;

private:
   std::string fName;
   std::string fObs;
   size_t fDim;
   double fVal;
};
#endif
