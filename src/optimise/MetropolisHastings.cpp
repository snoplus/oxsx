#include <MetropolisHastings.h>
#include <TestStatistic.h>
#include <Exceptions.h>
#include <HistTools.h>
#include <Formatter.hpp>
#include <Combinations.hpp>
#include <ContainerTools.hpp>
#include <Rand.h>
#include <iostream>
#include <math.h> //sqrt

const ParameterDict&
MetropolisHastings::GetSigmas() const{
    return fSigmas;
}

void
MetropolisHastings::SetSigmas(const ParameterDict& sigmas_){
    fSigmas = sigmas_;
}

unsigned
MetropolisHastings::GetBurnIn() const{
    return fBurnIn;
}

void
MetropolisHastings::SetBurnIn(unsigned burnIn_){
    fBurnIn = burnIn_;
}

unsigned
MetropolisHastings::GetThinFactor() const{
    return fThinFactor;
}

void
MetropolisHastings::SetThinFactor(unsigned thinFac_){
    fThinFactor = thinFac_;
}

unsigned
MetropolisHastings::GetMaxIter() const{
    return fMaxIter;
}

void
MetropolisHastings::SetMaxIter(unsigned maxIter_){
    fMaxIter = maxIter_;
}

double
MetropolisHastings::GetRejectionRate() const{
    return fRejectionRate;
}

ParameterDict
MetropolisHastings::GetMinima() const{
    return fMinima;
}

void
MetropolisHastings::SetMinima(const ParameterDict& minima_){
    fMinima = minima_;
}


ParameterDict
MetropolisHastings::GetMaxima() const{
    return fMaxima;
}

void
MetropolisHastings::SetMaxima(const ParameterDict& maxima_){
    fMaxima = maxima_;
}


void
MetropolisHastings::SetFlipSign(bool b_){
    fFlipSign = b_;
}

bool
MetropolisHastings::GetFlipSign() const{
    return fFlipSign;
}

bool
MetropolisHastings::GetTestStatLogged() const{
    return fTestStatLogged;
}

void
MetropolisHastings::SetTestStatLogged(bool b_){
    fTestStatLogged = b_;
}


bool 
MetropolisHastings::GetSaveFullHistogram() const{
    return fSaveFullHistogram;
}

void
MetropolisHastings::SetSaveFullHistogram(bool b_){
    fSaveFullHistogram = b_;
}

void
MetropolisHastings::SetHistogramAxes(const AxisCollection& template_){
    fHistogramAxes = template_;
}
             
AxisCollection
MetropolisHastings::GetHistogramAxes() const{
    return fHistogramAxes;
}                                                        

void
MetropolisHastings::SetInitialTrial(const ParameterDict& trial_){
    fInitialTrial = trial_;
}

ParameterDict
MetropolisHastings::GetInitialTrial() const{
    return fInitialTrial;
}

const FitResult&
MetropolisHastings::Optimise(TestStatistic* testStat_){
    testStat_->RegisterFitComponents();
    pTestStatistic = testStat_;
    fNDims = fMinima.size();
    fRejectionRate = 0;
    fMaxVal = 0;

    // Check initialisation
    size_t nParams = testStat_ -> GetParameterCount();
    if(   (fInitialTrial.size() != nParams && fInitialTrial.size())
          || fMaxima.size() != nParams
          || fMinima.size() != nParams
          || fSigmas.size() != nParams
          )
        throw LogicError(Formatter()
                         << "Metropolis Hastings initialisation error - Got "
                         << fMinima.size() << " Minima, "
                         << fMaxima.size() << " Maxima and "
                         << fSigmas.size() << " Sigmas"
                         << " - Need one per fit parameter (" << nParams
                         << ")"
                         << "\n Fit Parameters: \n"
                         << ContainerTools::ToString(testStat_->GetParameterNames())
                         << "\n If initial trial is specified it should be one per fit parameter"
                         );

    InitialiseHistograms();

    // 1. Choose a random starting point or the user defined one
    ParameterDict currentStep;
    std::set<std::string> parameterNames = testStat_->GetParameterNames();
    if(fInitialTrial.size())
      currentStep = fInitialTrial;
    else{
        for(std::set<std::string>::iterator it = parameterNames.begin(); 
        it != parameterNames.end(); ++it){
        const std::string& name = *it;
            currentStep[name] = fMinima[name] + Rand::Uniform() * (fMaxima[name] - fMinima[name]);
        }
}

    std::cout << "Metropolis Hastings::Initial Position @:" << std::endl;
    for(ParameterDict::iterator it = currentStep.begin(); it != currentStep.end(); ++it)
        std::cout << it->first << " : " << it->second << std::endl;
    std::cout << std::endl;

    fSample.empty();
    //    fSample.reserve(fMaxIter - fBurnIn);
   
    // 2. Loop step through the space a fixed number of times and
    for(unsigned i = 0; i < fMaxIter; i++){
        // a. Save the point in question if you are past burn-in phase and according to thinning
        if (i > fBurnIn && !(i%fThinFactor)){
          //            fSample.push_back(currentStep);

            if(fSaveFullHistogram)
                fHist.Fill(currentStep);

            else{
                FillProjections(currentStep);
            }
        }
      
        if(!(i%10000))
            std::cout << i << "  /  " << fMaxIter 
              << "\t" << (i - fRejectionRate) / static_cast<double>(i) << std::endl;

        // b. Propose a new step according to a random jump distribution
        ParameterDict proposedStep = JumpDraw(currentStep);

        // c. Test for acceptance
        if (StepAccepted(currentStep, proposedStep))
            currentStep = proposedStep;
        else
            fRejectionRate++;
    }
    fRejectionRate/= static_cast<double>(fMaxIter);
    std::cout << "Metropolis Hastings:: acceptance rate = " << 1 - fRejectionRate << std::endl;

    fFitResult.SetBestFit(fBestFit);
    fFitResult.SetStatSample(fSample);
    fFitResult.SetStatSpace(fHist);
    fFitResult.SetValid(true);
    if(!fSaveFullHistogram)
      SaveProjections();
    fFitResult.SetExtremeVal(fMaxVal);

    return fFitResult;
}

bool
MetropolisHastings::StepAccepted(const ParameterDict& thisStep_, 
                                 const ParameterDict& proposedStep_){

    // dont step outside of the fit region 
    for (ParameterDict::const_iterator it  = thisStep_.begin();
         it != thisStep_.end(); ++it){
        if (proposedStep_.at(it->first) < fMinima.at(it->first) || proposedStep_.at(it->first) > fMaxima.at(it->first))
            return false;
    }
    
    
    pTestStatistic -> SetParameters(thisStep_);
    double thisVal = pTestStatistic -> Evaluate();
    

    pTestStatistic -> SetParameters(proposedStep_);
    double proposedVal = pTestStatistic -> Evaluate();

    if(fFlipSign){
        thisVal = -thisVal;
        proposedVal = -proposedVal;
    }

    if(thisVal > fMaxVal || fMaxVal == 0.0){
        fMaxVal = thisVal;
        fBestFit = thisStep_;
    }
    
    double acceptanceParam = 0;
    if(fTestStatLogged)
        acceptanceParam = exp(proposedVal - thisVal);

    else
        acceptanceParam = proposedVal/thisVal;
    
    if (acceptanceParam > 1)
        return true;

    else if (acceptanceParam >= Rand::Uniform())
        return true;

    else{
        return false;
    }
}

ParameterDict
MetropolisHastings::JumpDraw(const ParameterDict& thisStep_) const{
    ParameterDict newStep;
    for(ParameterDict::const_iterator it = thisStep_.begin(); it != thisStep_.end(); ++it){
        if(fMinima.at(it->first) == fMaxima.at(it->first))
            newStep[it->first] = fMinima.at(it->first);
        
        else
            newStep[it->first] = it->second + Rand::Gaus(0, fSigmas.at(it->first));
    }
    return newStep;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// std::vector<double>                                                                             //
// MetropolisHastings::JumpDraw(const std::vector<double>& thisStep_) const{                       //
//     double stepSize = Rand::Gaus(0, fSigma);                                                    //
//                                                                                                 //
//     // choose random direction on NParam -1 Sphere according to                                 //
//     // http://mathworld.wolfram.com/HyperspherePointPicking.html (Muller 1959, Marsaglia 1972)  //
//                                                                                                 //
//     std::vector<double> newPoint;                                                               //
//     double norm = 0;                                                                            //
//     for(size_t i = 0; i < fNDims; i++){                                                         //
//         double component = Rand::Gaus();                                                        //
//         newPoint.push_back(component);                                                          //
//         norm += component * component;                                                          //
//     }                                                                                           //
//                                                                                                 //
//     norm = sqrt(norm);                                                                          //
//     for(size_t i = 0; i < fNDims; i++){                                                         //
//         newPoint[i] *= stepSize/norm;                                                           //
//     }                                                                                           //
//                                                                                                 //
//     // Add to orginal vector and return                                                         //
//     for(size_t i = 0; i < fNDims; i++)                                                          //
//         newPoint[i] += thisStep_.at(i);                                                         //
//                                                                                                 //
//     return newPoint;                                                                            //
// }                                                                                               //
/////////////////////////////////////////////////////////////////////////////////////////////////////


void
MetropolisHastings::InitialiseHistograms(){
    // Create histogram(s)
    // If the binning is specified, use that
    AxisCollection histAxes;
    if(fHistogramAxes.GetNDimensions()){
        if(fHistogramAxes.GetNDimensions() != pTestStatistic->GetParameterCount())
            throw DimensionError("MetropolisHastings::InitialiseHistograms",
                                 pTestStatistic->GetParameterCount(),
                                 fHistogramAxes.GetNDimensions(),
                                 "Template axes - one axis for each parameter"
                                 );
        histAxes = fHistogramAxes;
    }
    
    // otherwise take a guess
    else{
        std::set<std::string> paramNames = pTestStatistic->GetParameterNames();
        for(std::set<std::string>::iterator it =  paramNames.begin();
            it != paramNames.end(); ++it){
            const std::string& name = *it;
            double max = fMaxima.at(name);
            if(max == fMinima.at(name))
                max += 0.01;
            histAxes.AddAxis(BinAxis(name, fMinima.at(name), max,
                                     int(pow(fMaxIter, 1./fMinima.size()))));
        }
    }
    
    // Set up the histogram, either with a big ND histogram
    if(fSaveFullHistogram){
        fHist = Histogram(histAxes);
    }
    
    
    // otherwise just save the visualisable projections    
    else{
        f2DProjNames = Combinations::AllCombsNoDiag(pTestStatistic->GetParameterNames());
        f2DProjections = HistTools::MakeAllHists(histAxes, f2DProjNames);
        f1DProjections = HistTools::MakeAllHists(histAxes, pTestStatistic->GetParameterNames());
    }
}


void
MetropolisHastings::FillProjections(const ParameterDict& params_){
    // 1D
    HistTools::FillAllHists(f1DProjections, params_);
 
    // 2D
    HistTools::FillAllHists(f2DProjections, params_);
}

void
MetropolisHastings::SaveProjections(){
  std::map<std::string, Histogram> map1d;
  for(size_t i = 0; i < f1DProjections.size(); i++){
    const Histogram& hist = f1DProjections.at(i);
    std::string name = hist.GetAxes().GetAxis(0).GetName();
    map1d[name + "_projection"] = hist;
  }

  std::map<std::string, Histogram> map2d;
  for(size_t i = 0; i < f2DProjections.size(); i++){
    const Histogram& hist = f2DProjections.at(i);
    std::string name0 = hist.GetAxes().GetAxis(0).GetName();
    std::string name1 = hist.GetAxes().GetAxis(1).GetName();

    map2d[name0 + "_" + name1 + "_projection"] = hist;
  }

  fFitResult.Set1DProjections(map1d);
  fFitResult.Set2DProjections(map2d);
  
}
