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
    MCMC(MCSampler& s_) : fSampler(s_),
                          fSamples(this),
                          fMaxIter(100000), 
                          fMaxVal(0), fFlipSign(false), 
                          fTestStatLogged(false), pTestStatistic(NULL),
                          fCurrentVal(0.),
                          fSaveChain(false)
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

    bool GetSaveChain() const;
    void SetSaveChain(bool);

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
    bool fSaveChain;
};
#endif
