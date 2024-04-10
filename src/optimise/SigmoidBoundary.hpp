#include <cmath>

template<typename Energy>
void
SigmoidBoundary<Energy>::SetParameters(const ParameterDict& d_){
    fStat.SetParameters(d_);
    fParams = d_;    
}

template<typename Energy>
double
SigmoidBoundary<Energy>::Evaluate(){
    // the target
    double correctedEnergy = fStat.Evaluate();
    // the soft boundaries
    for(ParameterDict::iterator it = fMinima.begin(); it != fMinima.end();
        ++it){
        correctedEnergy += log(1 + exp(-fMu * (fParams[it->first] - it->second)));
    }

    
    for(ParameterDict::iterator it = fMaxima.begin(); it != fMaxima.end();
        ++it){
        correctedEnergy += log(1 + exp(+fMu * (fParams[it->first] - it->second)));
    }
    return correctedEnergy;
}
