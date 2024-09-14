#ifndef __OXSX_METROPOLIS_SAMPLER__
#define __OXSX_METROPOLIS_SAMPLER__
#include <MCSampler.h>
#include <ParameterDict.h>
#include <set>

class MetropolisSampler : public MCSampler
{
public:
   const ParameterDict &GetSigmas() const;
   void SetSigmas(const ParameterDict &);

   void Fix(const std::string &param_);
   void Release(const std::string &param_);

   ParameterDict Draw(const ParameterDict &current_);
   inline double CorrectAccParam() { return 0; }

private:
   ParameterDict fSigmas;
   std::set<std::string> fFixedParameters;
};
#endif
