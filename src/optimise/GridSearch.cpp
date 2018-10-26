#include <GridSearch.h>
#include <TestStatistic.h>
#include <iostream>
#include <sstream>
#include <Histogram.h>
#include <FitResult.h>
#include <Exceptions.h>
#include <Formatter.hpp>
#include <ContainerTools.hpp>

using ContainerTools::ToString;
using ContainerTools::HasSameKeys;
using ContainerTools::GetKeys;

void 
GridSearch::SetMinima(const ParameterDict& minima_){
    fMinima = minima_;
}

void 
GridSearch::SetMaxima(const ParameterDict& maxima_){
    fMaxima = maxima_;
}

void
GridSearch::SetStepSizes(const ParameterDict& steps_){
    fStepSizes = steps_;
}

ParameterDict 
GridSearch::GetMaxima() const{
    return fMaxima;
}

ParameterDict 
GridSearch::GetMinima() const{
    return fMinima;
}

ParameterDict
GridSearch::GetStepSizes() const{
    return fStepSizes;
}

const FitResult&
GridSearch::Optimise(TestStatistic* testStat_){
    // list of rates followed by list of systematics
    testStat_->RegisterFitComponents();

    // use this map to set the parameters, change the values in place
    ParameterDict setParams = testStat_ -> GetParameters();
    
    // check initialisation
    size_t nParams = testStat_ -> GetParameterCount();
    if( !HasSameKeys(fStepSizes, fMaxima)
        || !HasSameKeys(fStepSizes, fMinima)
        )
        throw LogicError(Formatter()
                         << "Grid Search initialisation error "
                         << " minima, maxima, stepsize parameters dont't match:\n"
                         << "Minima for :\n" << ToString(GetKeys(fMinima)) << "\n"
                         << "Maxima for :\n" << ToString(GetKeys(fMaxima)) << "\n"
                         << "StepSizes for :\n" << ToString(GetKeys(fStepSizes)) << "\n"
                         );

    // Prepare best fit
    ParameterDict bestFit;   
    fMinVal = 0;

    // Initialise LH space and count the number of grid steps
    unsigned maxSteps = 1;
    ParameterDict gridCounts;

    for(ParameterDict::iterator it = fMinima.begin(); it != fMinima.end(); ++it){
        size_t axisCounts = 1 + static_cast<size_t>((fMaxima[it->first] - it->second) / fStepSizes[it->first]);
        maxSteps *= axisCounts;
        gridCounts[it->first] = axisCounts;
    }    

    // start at min value
    fParamVals = fMinima;
    
    if(maxSteps == 1){
        std::cout << "Warning: Grid Search has only one grid point"
                  << std::endl;
    }
	
    // count interations
    unsigned stepCount = 0;
    int oneTenth = maxSteps/10;

    if(oneTenth < 10)
        oneTenth = 1;
    
    while(Increment(fStepSizes.begin(), fStepSizes.end())){
        // calculate the new value
        // if bigger, grab this as new best fit

        if(!(stepCount++ % oneTenth)){
            std::cout << stepCount << " / " << maxSteps  << "\t"
                      <<  100 * double(stepCount)/maxSteps 
                      << "%" << std::endl;
        }
        
        testStat_ -> SetParameters(fParamVals);
        double currentVal = testStat_ -> Evaluate();
        
        // if maximising, take the negative of the test stat
        if(fMaximising)
            currentVal *= -1;

        if (currentVal < fMinVal || !fMinVal){
            fMinVal = currentVal;
            bestFit = fParamVals;
        }
       
    }
    fFitResult.SetBestFit(bestFit);
    fFitResult.SetExtremeVal(fMinVal);
    return fFitResult;
}


bool 
GridSearch::Increment(ParameterDict::iterator it_, 
                      ParameterDict::iterator end_){
    const std::string paramToChange = it_->first;
    // nothing doing
    if(fMinima == fMaxima)
        return false;
    fParamVals[paramToChange] += fStepSizes.at(paramToChange);
    
    // wrap around past the maximum
    if (fParamVals[paramToChange] > fMaxima.at(paramToChange)){
        fParamVals[paramToChange] = fMinima.at(paramToChange);

        // if its the last index no rippling to do
        if (it_ == end_--)
            return false;

        // ripple up to next index
        if (Increment(++it_, end_)) // propagate the false down
            return true;
        else
            return false;

    }
    return true;
}

FitResult
GridSearch::GetFitResult() const{
    return fFitResult;
}

bool
GridSearch::GetMaximising() const{
    return fMaximising;
}

void
GridSearch::SetMaximising(bool b_){
    fMaximising = b_;
}
