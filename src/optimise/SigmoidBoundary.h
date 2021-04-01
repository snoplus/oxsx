#ifndef __OXSX_SIGMOID_BOUNDARY__
#define __OXSX_SIGMOID_BOUNDARY__
#include <ParameterDict.h>

template<typename Energy>
class SigmoidBoundary{
 public:
 SigmoidBoundary(Energy& stat_, const ParameterDict& minima_,
                 const ParameterDict& maxima_, double mu_ = 1000)
                : fStat(stat_), fMu(mu_), fMinima(minima_),
                  fMaxima(maxima_) {}
    
    void SetParameters(const ParameterDict& d_);
    double Evaluate();
    std::set<std::string> GetParameterNames() const {return fStat.GetParameterNames();}
    
 private:
    Energy&  fStat;
    double   fMu;
    ParameterDict fParams;
    ParameterDict fMinima;
    ParameterDict fMaxima;
};
#include <SigmoidBoundary.hpp>
#endif
