#ifndef __OXSX__LEAP_FROG__
#define __OXSX__LEAP_FROG__
#include <vector>
#include <ParameterDict.h>

class LeapFrog{
public:
    // one step
    template<typename T>
    static void Hamiltonian(ParameterDict& qs_, ParameterDict& ps_,
                            const ParameterDict& masses_,
                            double epsilon_, T&);

    // several
    template<typename T>
    static void Hamiltonian(ParameterDict& qs_, ParameterDict& ps_,
                            const ParameterDict& masses_,
                            double epsilon_, T&, int nSteps_);
    

    
    //  reflections
    template<typename T>
    static void Hamiltonian(ParameterDict& qs_, ParameterDict& ps_,
                            const ParameterDict& masses_,
                            double epsilon_, T&, int nSteps_,
                            const ParameterDict& minima_, const ParameterDict& maxima_);
    
    // look for reflections
    static double IntersectionTime(double wall_, double intPos_, double vel_);

    static void FirstIntersectionTime(const ParameterDict& initPos_,
                                      const ParameterDict& momenta_, const ParameterDict& masses_,
                                      const ParameterDict& minima_, const ParameterDict& maxima_,
                                      double epsilon_,
                                      std::string& dim_, double& time_
                                      );
};

#include <LeapFrog.hpp>
#endif
