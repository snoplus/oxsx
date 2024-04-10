#include <CountingResult.h>
#include <Formatter.hpp>
#include <Comparison.hpp>
#include <utility>
#include <algorithm>
#include <map>
#include <numeric>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <cmath>

void
CountingResult::AddBackground(double expectedRate_, const std::string& name_,
                              double error_, const CutLog& log_){
    fExpectedRates.push_back(expectedRate_);
    fExpectedRateErrors.push_back(error_);
    fBackgroundNames.push_back(name_);
    fBackgroundLogs.push_back(log_);
}

void
CountingResult::SetDataLog(const CutLog& log_){
    fDataLog = log_;
}

void
CountingResult::SetObservedCounts(int counts_){
    fObservedCounts = counts_;
}

double
CountingResult::GetExpectedCounts() const{
    return std::accumulate(fExpectedRates.begin(), fExpectedRates.end(),
                           0.);
}

double
CountingResult::GetExpectedCountError() const{
  double sqSum = 0;
  for(size_t i = 0; i < fExpectedRateErrors.size(); i++)
    sqSum += fExpectedRateErrors.at(i) * fExpectedRateErrors.at(i);
  return sqrt(sqSum);
}

double
CountingResult::GetSignalEfficiency() const{
  return fSignalEfficiency;
}

int 
CountingResult::GetObservedCounts() const{
    return fObservedCounts;
}

void 
CountingResult::SetSignal(double sigEff_, const std::string& name_, const CutLog& log_){
    fSignalName = name_;
    fSignalLog  = log_;
    fSignalEfficiency = sigEff_;
    
}

void
CountingResult::Print(){
    std::cout << AsString() << std::endl;
}


const std::vector<double>& 
CountingResult::GetExpectedRates() const{
    return fExpectedRates;
}

const std::vector<double>&
CountingResult::GetExpectedRateErrors() const{
    return fExpectedRateErrors;
}

void 
CountingResult::SaveAs(const std::string& fileName_){
    time_t t = time(0);
    struct tm * now = localtime(&t);

    std::ofstream fs(fileName_.c_str());
    fs << "-----------------------------------"
       << "-----------------------------------" << std::endl
       << "OXSX Counting Result: "
       << now -> tm_year + 1900 << '-' 
       << now -> tm_mon + 1 << '-'
       << now -> tm_mday
       << "  "
       << now -> tm_hour
       << ":"
       << now -> tm_min
       << ":"
       << now -> tm_sec
       << std::endl
       << "-----------------------------------"
       << "-----------------------------------" << std::endl
       << std::endl
       << AsString();       
}
std::string
CountingResult::AsString(){
    Formatter resultString;
    resultString << "------Backgrounds:----------------------------------------------------\n";
    for(size_t i = 0; i < fBackgroundLogs.size(); i++){
        resultString << "------------------------\n"
                     << fBackgroundNames.at(i) << "\n"
                     << "------------------------\n"
                     << fBackgroundLogs.at(i).AsString()
                     << "\n\n";
        
    }

    resultString << "------Signal----------------------------------------------------------\n";
    resultString << "------------------------\n"
                 << fSignalName << "\n"
                 << "------------------------\n"
                 << fSignalLog.AsString()
                 << "\n\n";
    
    resultString << "------Test:-----------------------------------------------------------\n"
                 << fDataLog.AsString()
                 << "\n\n";
    
    resultString << "------Summary:--------------------------------------------------------\n";
    resultString << "------------------------\n" 
                 << std::left
                 << std::setw(25)
                 << "Background"
                 << std::left
                 << std::setw(10)
                 << "Counts";
    
    std::vector<std::pair<std::string, double> > sortedRatesByName;
    std::map<std::string, double> sortedErrorsByName;
    for(size_t i = 0; i < fExpectedRates.size(); i++){
      sortedRatesByName.push_back(std::pair<std::string, double>(fBackgroundNames.at(i), fExpectedRates.at(i)));
      sortedErrorsByName[fBackgroundNames.at(i)] = fExpectedRateErrors.at(i);
    }   
    std::sort(sortedRatesByName.begin(), sortedRatesByName.end(), PairSortBySecond<std::string, double, std::greater<double> >());

    for(size_t i = 0; i < sortedRatesByName.size(); i++){
        resultString << "\n"
                     << std::setw(25)
                     << std::left
                     << sortedRatesByName.at(i).first
                     << std::setw(10)
                     << std::left
                     << sortedRatesByName.at(i).second
                     << " +- " 
                     << sortedErrorsByName[sortedRatesByName.at(i).first];
    }
        
    resultString << "\n"
                 << std::setw(25)
                 << std::left
                 << "Total :" 
                 << std::setw(10)
                 << std::left
                 << GetExpectedCounts()
                 << " +- "
                 << GetExpectedCountError()
                 << "\n------------------------\n"
                 << std::setw(18)
                 << std::left
                 << "Signal Efficiency"
                 << std::setw(10)
                 << std::left
                 << fSignalEfficiency
                 << "\n------------------------\n"
                 << std::setw(18)
                 << std::left
                 << "Observed: "
                 << std::setw(10)
                 << std::left
                 << fObservedCounts
                 << "\n"
                 << "----------------------------------------------------------------------\n";
    return resultString;
}
