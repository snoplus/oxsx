#ifndef __OXSX_METROPOLIS_HASTINGS__
#define __OXSX_METROPOLIS_HASTINGS__
#include <Optimiser.h>
#include <FitResult.h>
#include <Histogram.h>
#include <set>

class TestStatistic;
class MetropolisHastings : public Optimiser{
 public:
     MetropolisHastings() : fBurnIn(3000), 
                            fThinFactor(1),
                            fMaxIter(100000), 
                            fTestStatLogged(false),
                            fFlipSign(false), 
                            fSaveFullHistogram(false),
                            pTestStatistic(NULL),
                            fMaxVal(0)
                            {}               
    
    const FitResult& Optimise(TestStatistic*); 

    unsigned GetBurnIn() const;
    void     SetBurnIn(unsigned);
    
    unsigned GetThinFactor() const;
    void     SetThinFactor(unsigned);
    
    unsigned GetMaxIter() const;
    void     SetMaxIter(unsigned);

    double   GetRejectionRate() const;

    ParameterDict GetMaxima() const;
    void SetMaxima(const ParameterDict&);

    ParameterDict GetMinima() const;
    void SetMinima(const ParameterDict&);

    const ParameterDict& GetSigmas() const;
    void   SetSigmas(const ParameterDict&);

    bool GetFlipSign() const;
    void SetFlipSign(bool);
    
    bool GetTestStatLogged() const;
    void SetTestStatLogged(bool b_);

    bool GetSaveFullHistogram() const;
    void SetSaveFullHistogram(bool);
    void SetHistogramAxes(const AxisCollection&);
    AxisCollection GetHistogramAxes() const;

    void SetInitialTrial(const ParameterDict&);
    ParameterDict GetInitialTrial() const;

 private:
    // configuration
    unsigned  fBurnIn;
    unsigned  fThinFactor;
    unsigned  fMaxIter;
    unsigned  fNDims;
    bool      fTestStatLogged;
    bool      fFlipSign;
    bool      fSaveFullHistogram;
    AxisCollection fHistogramAxes;

    ParameterDict fMaxima;
    ParameterDict fMinima;
    ParameterDict fSigmas;
    ParameterDict fInitialTrial;

    // internal copy
    TestStatistic* pTestStatistic;

    // results
    double fRejectionRate;

    std::vector<Histogram>           f1DProjections;
    std::vector<Histogram>           f2DProjections;
    std::set<std::pair<std::string, std::string> > f2DProjNames;

    Histogram fHist;
    std::vector< std::vector<double> > fSample;

    FitResult fFitResult;
    ParameterDict fBestFit;
    ParameterDict fSetParameters;
    double fMaxVal;
    
    // private functions
    ParameterDict JumpDraw(const ParameterDict& thisStep_) const; 
    inline double JumpProbRatio(const ParameterDict& thisStep_, 
                                const ParameterDict& proposedStep_) const {return 1;}

    bool   StepAccepted(const ParameterDict& thisStep_,
                        const ParameterDict& proposedStep_);
    
    void    InitialiseHistograms();
    void    FillProjections(const ParameterDict&);
    void    SaveProjections();
};
#endif
