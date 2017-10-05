
template<typename T>
void
LeapFrog::Hamiltonian(ParameterDict& qs_, ParameterDict& ps_,
                      const ParameterDict& masses_,
                      double epsilon_, T& differ_){
    // step 1
    ParameterDict grad = differ_(qs_);
    for(ParameterDict::iterator it = ps_.begin(); it != ps_.end(); ++it)
        it->second -= epsilon_/2 * grad[it->first];

    // step 2
    for(ParameterDict::iterator it = qs_.begin(); it != qs_.end(); ++it)
        it->second += epsilon_ * ps_[it->first]/masses_.at(it->first);

    // step 3
    grad = differ_(qs_);
    for(ParameterDict::iterator it = ps_.begin(); it != ps_.end(); ++it)
        it->second -= epsilon_/2 * grad[it->first];


}


template<typename T>
void
LeapFrog::Hamiltonian(ParameterDict& qs_, ParameterDict& ps_,
                      const ParameterDict& masses_,
                      double epsilon_, T& differ_, int nSteps_){
    // step 1 take a half step in momentum
    ParameterDict grad = differ_(qs_);
    for(ParameterDict::iterator it = ps_.begin(); it != ps_.end(); ++it)
        it->second -= epsilon_/2 * grad[it->first];

    // step 2 alternate full steps for position and momentum
    for(int i = 0; i < nSteps_; i++){
        // Full step for position
        for(ParameterDict::iterator it = qs_.begin(); it != qs_.end(); ++it)
            it->second += epsilon_ * ps_[it->first]/masses_.at(it->first);
        
        // Full step for momentum, unless this is the last step
        if(i != nSteps_ -1){
            grad = differ_(qs_);
            for(ParameterDict::iterator it = ps_.begin(); it != ps_.end(); ++it)
                it->second -= epsilon_/2 * grad[it->first];
        }
        
    }

    // step 3 make a half step for momentum
    grad = differ_(qs_);
    for(ParameterDict::iterator it = ps_.begin(); it != ps_.end(); ++it)
        it->second -= epsilon_/2 * grad[it->first];    
}

