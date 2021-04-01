#include <MCMC.h>
#include <TestStatistic.h>
#include <Exceptions.h>
#include <HistTools.h>
#include <Formatter.hpp>
#include <Combinations.hpp>
#include <ContainerTools.hpp>
#include <Rand.h>
#include <iostream>
#include <math.h> //isinf

void
MCMC::SetHistogramAxes(const AxisCollection& ax_){
    fSamples.SetHistogramAxes(ax_);
}
AxisCollection
MCMC::GetHistogramAxes() const{
    return fSamples.GetHistogramAxes();
}

unsigned
MCMC::GetBurnIn() const{
    return fSamples.GetBurnIn();
}


void
MCMC::SetBurnIn(unsigned burnIn_){
    fSamples.SetBurnIn(burnIn_);
}


unsigned
MCMC::GetThinFactor() const{
    return fSamples.GetThinFactor();
}


void
MCMC::SetThinFactor(unsigned thinFac_){
    fSamples.SetThinFactor(thinFac_);
}


unsigned
MCMC::GetMaxIter() const{
    return fMaxIter;
}


void
MCMC::SetMaxIter(unsigned maxIter_){
    fMaxIter = maxIter_;
}


const MCMCSamples&
MCMC::GetSamples() const{
    return fSamples;
}


ParameterDict
MCMC::GetMinima() const{
    return fMinima;
}



void
MCMC::SetMinima(const ParameterDict& minima_){
    fMinima = minima_;
}


ParameterDict
MCMC::GetMaxima() const{
    return fMaxima;
}


void
MCMC::SetMaxima(const ParameterDict& maxima_){
    fMaxima = maxima_;
}


void
MCMC::SetFlipSign(bool b_){
    fFlipSign = b_;
}


bool
MCMC::GetFlipSign() const{
    return fFlipSign;
}


bool
MCMC::GetTestStatLogged() const{
    return fTestStatLogged;
}


void
MCMC::SetTestStatLogged(bool b_){
    fTestStatLogged = b_;
}


bool 
MCMC::GetSaveFullHistogram() const{
    return fSamples.GetSaveFullHistogram();
}


void
MCMC::SetSaveFullHistogram(bool b_){
    fSamples.SetSaveFullHistogram(b_);
}


void
MCMC::SetInitialTrial(const ParameterDict& trial_){
    fInitialTrial = trial_;
}


ParameterDict
MCMC::GetInitialTrial() const{
    return fInitialTrial;
}


const FitResult&
MCMC::Optimise(TestStatistic* testStat_){
    testStat_->RegisterFitComponents();
    pTestStatistic = testStat_;
    fNDims = fMinima.size();

    // Check initialisation
    size_t nParams = testStat_ -> GetParameterCount();
    if(   (fInitialTrial.size() != nParams && fInitialTrial.size())
          || fMaxima.size() != nParams
          || fMinima.size() != nParams
          )
        throw LogicError(Formatter()
                         << "MCMC initialisation error - Got "
                         << fMinima.size() << " Minima, "
                         << fMaxima.size() << " Maxima and "
                         << " - Need one per fit parameter (" << nParams
                         << ")"
                         << "\n Fit Parameters: \n"
                         << ContainerTools::ToString(testStat_->GetParameterNames())
                         << "\n If initial trial is specified it should be one per fit parameter"
                         );

    // 1. Choose a random starting point or the user defined one
    std::set<std::string> parameterNames = testStat_->GetParameterNames();
    if(fInitialTrial.size())
        fCurrentStep = fInitialTrial;
    else{
        for(std::set<std::string>::iterator it = parameterNames.begin(); 
        it != parameterNames.end(); ++it){
        const std::string& name = *it;
        fCurrentStep[name] = fMinima[name] + Rand::Uniform() * (fMaxima[name] - fMinima[name]);
        }
    }

    std::cout << "MCMC::Initial Position @:" << std::endl;
    for(ParameterDict::iterator it = fCurrentStep.begin(); it != fCurrentStep.end(); ++it)
        std::cout << it->first << " : " << it->second << std::endl;
    std::cout << std::endl;

    pTestStatistic->SetParameters(fCurrentStep);
    fCurrentVal = pTestStatistic->Evaluate();
    if(fFlipSign)
        fCurrentVal *= -1;

    fBestFit = fCurrentStep;
    fMaxVal = fCurrentVal;
    
    // 2. Loop step through the space a fixed number of times
    for(unsigned i = 0; i < fMaxIter; i++){      
        if(!(i%100000) && i)
            std::cout << i << "  /  " << fMaxIter
                      << "\t" << fSamples.GetAcceptanceRate()
                      << "\t" << fSamples.GetAutoCorrelations()[1]
                      << "\t" << fSamples.GetAutoCorrelations()[10]
                      << "\t" << fSamples.GetAutoCorrelations()[100]
                      << std::endl;


        // b. Propose a new step according
        ParameterDict proposedStep = fSampler.Draw(fCurrentStep);

        // c. Decide whethere to step there or not
        bool accepted =  StepAccepted(proposedStep);

        // d. log
        fSamples.Fill(fCurrentStep, fCurrentVal, accepted);
    }
    std::cout << "MCMC:: acceptance rate = " << fSamples.GetAcceptanceRate()
              << std::endl;

    fFitResult.SetBestFit(fBestFit);
    fFitResult.SetValid(true);
    fFitResult.SetExtremeVal(fMaxVal);

    return fFitResult;
}

bool
MCMC::StepAccepted(const ParameterDict& proposedStep_){
    // dont step outside of the fit region 
    for (ParameterDict::const_iterator it  = fCurrentStep.begin();
         it != fCurrentStep.end(); ++it){
        if (proposedStep_.at(it->first) < fMinima.at(it->first) || proposedStep_.at(it->first) > fMaxima.at(it->first))
            return false;
    }

    pTestStatistic -> SetParameters(proposedStep_);
    double proposedVal = pTestStatistic -> Evaluate();
    
    if(fFlipSign){
        proposedVal = -proposedVal;
    }    
   
    if(fCurrentVal > fMaxVal){
        fMaxVal = fCurrentVal;
        fBestFit = fCurrentStep;
    }    

    long double acceptanceParam = 0;
    if(fTestStatLogged)
        acceptanceParam = exp(proposedVal - fCurrentVal + fSampler.CorrectAccParam(acceptanceParam));
    else
        acceptanceParam = proposedVal/fCurrentVal;


    bool accept = false;
    if (acceptanceParam > 1 || isinf(acceptanceParam))
        accept =  true;
    else if (acceptanceParam >= Rand::Uniform())
        accept = true;

    if(accept){
        fCurrentVal = proposedVal;
        fCurrentStep = proposedStep_;
    }

    return accept;
}
