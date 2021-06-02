#ifndef __OXSX_MCMCSAMPLES__
#define __OXSX_MCMCSAMPLES__
#include <AxisCollection.h>
#include <Histogram.h>
#include <ParameterDict.h>
#include <AutoCorrelationCalc.h>
#include "TStopwatch.h"
#include <TTree.h>

class MCMC;
class MCMCSamples{
 public:
    MCMCSamples(MCMC* p_) : fSaveFullHistogram(false), fAcceptedSteps(0), fTotalSteps(0),
                            fBurnIn(1000), fThinFactor(1), fInitialised(false), fMCMC(p_),
                            fAutoCorrelator(1000),
                            fSaveChain(false)
                            {}


    int  GetBurnIn() const;
    void SetBurnIn(int);

    int  GetThinFactor() const;
    void SetThinFactor(int);
    
    void SetHistogramAxes(const AxisCollection&);
    AxisCollection GetHistogramAxes() const;

    bool GetSaveFullHistogram() const;
    void SetSaveFullHistogram(bool);

    bool GetSaveChain() const;
    void SetSaveChain(bool);

    double   GetRejectionRate() const;
    double   GetAcceptanceRate() const;

    TTree* GetChain() const;

    void Fill(const ParameterDict&, double val_, bool accepted_);

    const std::map<std::string, Histogram>& Get1DProjections() const;
    const std::map<std::string, Histogram>& Get2DProjections() const;

    const Histogram& GetHistogram() const;
    const std::vector< std::vector<double> >& GetRawSamples() const;

    void Clear();

    std::vector<double> GetAutoCorrelations();

 private:
    AutoCorrelationCalc fAutoCorrelator;
    MCMC*   fMCMC;
    bool    fInitialised;
    int     fBurnIn;
    int     fThinFactor;
    void    InitialiseHistograms();
    void    FillProjections(const ParameterDict&);
    
    // results
    bool      fSaveFullHistogram;
    AxisCollection fHistogramAxes;

    int fAcceptedSteps;
    int fTotalSteps;

    std::map<std::string, Histogram>           f1DProjections;
    std::map<std::string, Histogram>           f2DProjections;

    Histogram fHist;
    std::vector< std::vector<double> > fSample;

    TTree *fChain;
    bool fAccepted;
    double fCurrentVal;
    int fStepNumber;
    double fStepTime;
    TStopwatch stepClock;
    bool fSaveChain;
    std::vector<double> parvals;
};

#endif
