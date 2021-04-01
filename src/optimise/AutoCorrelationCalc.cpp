#include <numeric>
#include <AutoCorrelationCalc.h>
#include <iostream>
void
AutoCorrelationCalc::Fill(double v_){
    // shuffle everything along
    fValues.push_back(v_);
    if(fValues.size() > fQLen)
        fValues.pop_front();

    // update the estimates
    double mean = Mean();

    double newest = fValues.back();
    for(size_t i = 0; i < fValues.size(); i++){
        double older = fValues.rbegin()[i];
        fDiffs[i] += (newest - mean) * (older - mean);
        fNorms[i] += (newest - mean) * (newest - mean);
    }    
}

std::vector<double>
AutoCorrelationCalc::Get() const{
    std::vector<double> autocorrs;
    autocorrs.reserve(fValues.size());

    for(size_t i = 0; i < fDiffs.size(); i++){
        autocorrs.push_back(fDiffs[i]/fNorms[i]);
    }

    return autocorrs;
}

double
AutoCorrelationCalc::Mean() const{
    return std::accumulate(fValues.begin(), fValues.end(), 0.0)/fValues.size();
}

double
AutoCorrelationCalc::Var() const{
    double mean = Mean();
    double sumSq = std::inner_product(fValues.begin(), fValues.end(), fValues.begin(), 0);
    return sumSq - mean * mean;
}

void
AutoCorrelationCalc::Clear(){
    fDiffs = std::vector<double>(fQLen, 0);
    fNorms = std::vector<double>(fQLen, 0);
    fValues.clear();
}
