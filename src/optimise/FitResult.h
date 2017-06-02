/****************************************************/
/* Output of a fit routine, used for limit setting. */
/****************************************************/
#ifndef __OXSX_FIT_RESULT__
#define __OXSX_FIT_RESULT__
#include <Histogram.h>
#include <cstdlib>
#include <vector>
#include <ParameterDict.h>
#include <string>
#include <map>

typedef std::map<std::string, Histogram> HistMap;

class FitResult{
 public:
    FitResult() : fStatSpace(NULL), fIsValid(true),  printPrecision(5) {}
    FitResult(const FitResult&); //deep copy
    FitResult operator=(const FitResult&); //deep copy
    ~FitResult(); // frees stat space

    const ParameterDict& GetBestFit() const;
    void SetBestFit(const ParameterDict&);
    
    void SetStatSpace(const Histogram&);
    const Histogram& GetStatSpace() const; 

    void SetStatSample(const std::vector<std::vector<double> >&);
    const std::vector<std::vector<double> >& GetStatSample() const;
    
    void SetValid(bool b_);
    bool GetValid() const;

    void SetPrintPrecision(const size_t&);
    size_t GetPrintPrecision() const;
    
    std::string AsString() const;
    void Print() const;
    void SaveAs(const std::string&) const;

    void Set1DProjections(const HistMap&);
    void Set2DProjections(const HistMap&);

    HistMap Get1DProjections() const;
    HistMap Get2DProjections() const;


 private:
    ParameterDict fBestFit;
    std::vector<std::vector<double> > fStatSample;
    Histogram*    fStatSpace;
    HistMap f1DProjections;
    HistMap f2DProjections;
   
    size_t printPrecision;
    bool fIsValid;
};
#endif

