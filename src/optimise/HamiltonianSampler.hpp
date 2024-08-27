#include <Rand.h>
#include <ContainerTools.hpp>
#include <LeapFrog.h>
#include <iostream>
#include <cmath>

template <typename StatType>
HamiltonianSampler<StatType>::HamiltonianSampler(StatType &stat_, double epsilon_, int nSteps_)
    : fEpsilon(epsilon_), fNSteps(nSteps_), fDiff(GradType(stat_, 0.05 * epsilon_)), fCorr(0)
{

    typedef std::set<std::string> StringSet;
    StringSet names = stat_.GetParameterNames();
    for (StringSet::iterator it = names.begin(); it != names.end(); ++it)
        fMasses[*it] = 1;
}

template <typename StatType>
ParameterDict
HamiltonianSampler<StatType>::Draw(const ParameterDict &thisStep_)
{
    // step 1: sample from masses to create momenta
    ParameterDict momenta;

    for (ParameterDict::const_iterator it = thisStep_.begin(); it != thisStep_.end(); ++it)
    {
        momenta[it->first] = Rand::Gaus(0, sqrt(fMasses[it->first]));
    }

    fCorr = KineticEnergy(momenta);

    // step 2: Hamiltonian dynamics, including reflections or not
    ParameterDict nextStep = thisStep_;
    if (fMinima.size() && fMaxima.size())
        LeapFrog::Hamiltonian(nextStep, momenta, fMasses, fEpsilon, fDiff, fNSteps, fMinima, fMaxima);
    else
        LeapFrog::Hamiltonian(nextStep, momenta, fMasses, fEpsilon, fDiff, fNSteps);

    fCorr -= KineticEnergy(momenta);

    return nextStep;
}

template <typename StatType>
void HamiltonianSampler<StatType>::SetMinima(const ParameterDict &m_)
{
    fMinima = m_;
}

template <typename StatType>
ParameterDict
HamiltonianSampler<StatType>::GetMinima() const
{
    return fMinima;
}

template <typename StatType>
void HamiltonianSampler<StatType>::SetMaxima(const ParameterDict &m_)
{
    fMaxima = m_;
}

template <typename StatType>
ParameterDict
HamiltonianSampler<StatType>::GetMaxima() const
{
    return fMaxima;
}

template <typename StatType>
double
HamiltonianSampler<StatType>::KineticEnergy(const ParameterDict &momenta_) const
{
    double kinEnergy = 0;
    for (ParameterDict::const_iterator it = momenta_.begin(); it != momenta_.end(); ++it)
        kinEnergy += it->second * it->second / 2 / fMasses.at(it->first);

    return kinEnergy;
}

template <typename StatType>
double
HamiltonianSampler<StatType>::CorrectAccParam()
{
    return fCorr;
}
