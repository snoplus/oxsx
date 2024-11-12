#include <MetropolisSampler.h>
#include <Rand.h>

const ParameterDict &
MetropolisSampler::GetSigmas() const
{
    return fSigmas;
}

void MetropolisSampler::SetSigmas(const ParameterDict &sig_)
{
    fSigmas = sig_;
}

void MetropolisSampler::Fix(const std::string &name_)
{
    fFixedParameters.insert(name_);
}

void MetropolisSampler::Release(const std::string &name_)
{
    fFixedParameters.erase(name_);
}

ParameterDict
MetropolisSampler::Draw(const ParameterDict &current_)
{
    ParameterDict newStep;
    for (ParameterDict::const_iterator it = current_.begin(); it != current_.end(); ++it)
    {
        // If the parameter is fixed, don't change it
        if (fFixedParameters.count(it->first))
            newStep[it->first] = it->second;
        else
            newStep[it->first] = it->second + Rand::Gaus(0, fSigmas.at(it->first));
    }
    return newStep;
}
