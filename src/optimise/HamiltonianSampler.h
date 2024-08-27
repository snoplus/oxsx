#ifndef __OXSX__HAMILTONIAN_SAMPLER__
#define __OXSX__HAMILTONIAN_SAMPLER__
#include <ParameterDict.h>
#include <MCSampler.h>
#include <Gradient.h>
#include <vector>

template <typename StatType>
class HamiltonianSampler : public MCSampler
{
public:
    HamiltonianSampler(StatType &stat_, double epsilon_, int nSteps_);

    ParameterDict Draw(const ParameterDict &thisStep_);

    void SetMasses(const ParameterDict &m_) { fMasses = m_; }
    ParameterDict GetMasses() const { return fMasses; }

    void SetMinima(const ParameterDict &m_);
    ParameterDict GetMinima() const;

    void SetMaxima(const ParameterDict &m_);
    ParameterDict GetMaxima() const;

    double KineticEnergy(const ParameterDict &momenta_) const;
    double CorrectAccParam();

private:
    typedef Gradient<StatType> GradType;
    ParameterDict fMasses;
    ParameterDict fMinima;
    ParameterDict fMaxima;

    double fCorr;

    int fNSteps;
    GradType fDiff;
    double fEpsilon;
};
#include <HamiltonianSampler.hpp>
#endif
