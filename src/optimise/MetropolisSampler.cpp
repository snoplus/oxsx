#include <MetropolisSampler.h>
#include <Rand.h>

const ParameterDict&
MetropolisSampler::GetSigmas() const{
    return fSigmas;
}

void
MetropolisSampler::SetSigmas(const ParameterDict& sig_){
    fSigmas = sig_;
}


ParameterDict
MetropolisSampler::Draw(const ParameterDict& current_){
    ParameterDict newStep;
    for(ParameterDict::const_iterator it = current_.begin(); it != current_.end(); ++it){
        newStep[it->first] = it->second + Rand::Gaus(0, fSigmas.at(it->first));
    }
    return newStep;    
}

