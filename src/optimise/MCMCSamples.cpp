#include <MCMCSamples.h>
#include <ParameterDict.h>
#include <HistTools.h>
#include <MCMC.h>
#include <ContainerTools.hpp>

void
MCMCSamples::SetSaveFullHistogram(bool b_){
    fSaveFullHistogram = b_;
}

bool
MCMCSamples::GetSaveFullHistogram() const{
    return fSaveFullHistogram;
}

void
MCMCSamples::FillProjections(const ParameterDict& params_){
    // 1D
    HistTools::FillAllHists(f1DProjections, params_);
 
    // 2D
    HistTools::FillAllHists(f2DProjections, params_);
}


double
MCMCSamples::GetRejectionRate() const{
    return 1 - GetAcceptanceRate();
}

double
MCMCSamples::GetAcceptanceRate() const{
    return (double)fAcceptedSteps/fTotalSteps;
}


void
MCMCSamples::SetHistogramAxes(const AxisCollection& ax_){
    fHistogramAxes = ax_;
}
AxisCollection
MCMCSamples::GetHistogramAxes() const{
    return fHistogramAxes;
}

void
MCMCSamples::InitialiseHistograms(){
    std::set<std::string> paramNames = ContainerTools::GetKeys(fMCMC->GetMinima());
    // Create histogram(s)
    // If the binning is specified, use that
    AxisCollection histAxes;
    if(fHistogramAxes.GetNDimensions()){
        if(fHistogramAxes.GetNDimensions() != fMCMC->GetMinima().size())
            throw DimensionError("MCMC::InitialiseHistograms",
                                 fMCMC->GetMinima().size(),
                                 fHistogramAxes.GetNDimensions(),
                                 "Template axes - one axis for each parameter"
                                 );
        histAxes = fHistogramAxes;
    }

    // otherwise take a guess
    else{
        for(std::set<std::string>::iterator it =  paramNames.begin();
            it != paramNames.end(); ++it){
            const std::string& name = *it;
            double max = fMCMC->GetMaxima().at(name);
            if(max == fMCMC->GetMinima().at(name))
                max += 0.01;
            histAxes.AddAxis(BinAxis(name, fMCMC->GetMinima().at(name), max,
                                     int(pow(fMCMC->GetMaxIter(), 1./fMCMC->GetMinima().size()))));
        }
    }
    // Set up the histogram, either with a big ND histogram
    if(fSaveFullHistogram){
        fHist = Histogram(histAxes);
    }
    
    
    // otherwise just save the visualisable projections    
    else{
        std::set<std::pair<std::string, std::string> > twoDProjNames;
        twoDProjNames = Combinations::AllCombsNoDiag(paramNames);
        f2DProjections = HistTools::MakeAllHists(histAxes, twoDProjNames);
        f1DProjections = HistTools::MakeAllHists(histAxes, paramNames);
    }

    fInitialised = true;
}

void
MCMCSamples::Fill(const ParameterDict& params_, double val_, bool accepted_){
    if(!fInitialised)
        InitialiseHistograms();

    // a. Save the point in question if you are past burn-in phase and according to thinning
    if (fTotalSteps > fBurnIn && !(fTotalSteps % fThinFactor)){
        //            fSample.push_back(currentStep);        
        if(fSaveFullHistogram)
            fHist.Fill(params_);
        else
            FillProjections(params_);

        fAutoCorrelator.Fill(val_);
    }
    if(accepted_)
        fAcceptedSteps++;
    
    fTotalSteps++;
}


const std::map<std::string, Histogram>&
MCMCSamples::Get1DProjections() const{
    return f1DProjections;
}
    
const std::map<std::string, Histogram>&
MCMCSamples::Get2DProjections() const{
    return f2DProjections;
}
                                                       

const Histogram&
MCMCSamples::GetHistogram() const{
    return fHist;   
}

const std::vector< std::vector<double> >&
MCMCSamples::GetRawSamples() const{
    return fSample;
}

int
MCMCSamples::GetBurnIn() const{
    return fBurnIn;
}

void
MCMCSamples::SetBurnIn(int b_){
    fBurnIn = b_;
}

int
MCMCSamples::GetThinFactor() const{
    return fThinFactor;
}

void
MCMCSamples::SetThinFactor(int tf_){
    fThinFactor = tf_;
}

void
MCMCSamples::Clear(){
    fAcceptedSteps = 0;
    fTotalSteps = 0;
    fSample.clear();
    f1DProjections.clear();
    f2DProjections.clear();
    fInitialised = false;
    fHist.Empty();
    fAutoCorrelator.Clear();
}

std::vector<double>
MCMCSamples::GetAutoCorrelations(){
    return fAutoCorrelator.Get();
}
