#ifndef __OXSX_MCMC__
#define __OXSX_MCMC__
#include <Optimiser.h>
#include <FitResult.h>
#include <MCMCSamples.h>
#include <MCSampler.h>
#include <set>

class TestStatistic;
class MCMC : public Optimiser{
 public:
    MCMC(MCSampler& s_) : fMaxIter(100000), 
                          fTestStatLogged(false),
                          fFlipSign(false),
                          fMaxVal(0),
                          fCurrentVal(0.),
                          fSamples(this),
                          pTestStatistic(NULL),
                          fSampler(s_)
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

    const MCMCSamples& GetSamples() const;

 private:
    // configuration
    unsigned  fMaxIter;
    unsigned  fNDims;
    bool      fTestStatLogged;
    bool      fFlipSign;
    
    double fMaxVal;
    double fCurrentVal;
    FitResult fFitResult;
    ParameterDict fBestFit;
    ParameterDict fCurrentStep;
    
    MCMCSamples fSamples;
    
    ParameterDict fMaxima;
    ParameterDict fMinima;
    ParameterDict fInitialTrial;

    // internal copy
    TestStatistic* pTestStatistic;
    
    MCSampler& fSampler;
    
    bool   StepAccepted(const ParameterDict& proposedStep_);
};
#endif
