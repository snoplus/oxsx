set(optimise_src
    optimise/AutoCorrelationCalc.cpp
    optimise/FitResult.cpp
    optimise/GridSearch.cpp
    optimise/LeapFrog.cpp
    optimise/MCMC.cpp
    optimise/MCMCSamples.cpp
    optimise/MetropolisSampler.cpp
    optimise/Minuit.cpp
    optimise/MinuitFCN.cpp
    PARENT_SCOPE
)

set(optimise_headers
    AutoCorrelationCalc.h
    FitResult.h
    Gradient.h
    Gradient.hpp
    GridSearch.h
    HamiltonianSampler.h
    HamiltonianSampler.hpp
    HMCEnergy.h
    HMCEnergy.hpp
    LeapFrog.h
    LeapFrog.hpp
    MCMC.h
    MCMCSamples.h
    MCSampler.h
    MetropolisSampler.h
    Minuit.h
    MinuitFCN.h
    Optimiser.h
    SigmoidBoundary.h
    SigmoidBoundary.hpp
)

file(COPY ${optimise_headers} DESTINATION ${OXSX_INCLUDE_DIR})