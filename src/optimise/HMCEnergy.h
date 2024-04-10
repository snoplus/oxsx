#ifndef __OXSX_HMC_ENERGY__
#define __OXSX_HMC_ENERGY__

template<typename LogLikelihood>
class HMCEnergy{
 public:
    HMCEnergy(LogLikelihood& lh_) : fLLikelihood(lh_) {}
    double Evaluate();
    void   SetParameters(const ParameterDict& p_);
    std::set<std::string> GetParameterNames() const;

 private:
    LogLikelihood& fLLikelihood;
};
#include <HMCEnergy.hpp>
#endif
