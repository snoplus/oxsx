#ifndef __OXSX_AUTO_CORRELATION_CALC__
#define __OXSX_AUTO_CORRELATION_CALC__
#include <deque>
#include <vector>

class AutoCorrelationCalc{
 public:
 AutoCorrelationCalc(size_t nVals_) : fDiffs(std::vector<double>(nVals_, 0)),
                                      fNorms(std::vector<double>(nVals_, 0)),
                                      fQLen(nVals_)
                                      {}

    void Fill(double);
    double Mean() const;
    double Var() const;
    
    std::vector<double> Get() const;
    void Clear();
    
 private:
    std::deque<double>    fValues;
    std::vector<double>   fDiffs;
    std::vector<double>   fNorms;
    size_t fQLen;
};
#endif
