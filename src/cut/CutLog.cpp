#include <CutLog.h>
#include <Exceptions.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Exceptions.h>

void
CutLog::CalculateMeta(){
    fRemainderCounts.clear();
    fRemainderPercentages.clear();
    fCutPercentages.clear();

    fRemainderCounts.reserve(fCutCounts.size());
    fRemainderPercentages.reserve(fCutCounts.size());
    fCutPercentages.reserve(fCutCounts.size());

    int numRemaining = fNEvents;
    for(size_t i = 0; i < fCutCounts.size(); i++){
        numRemaining -= fCutCounts.at(i);
        fRemainderCounts.push_back(numRemaining);

        fCutPercentages.push_back(fCutCounts.at(i) * 100./fNEvents);
        fRemainderPercentages.push_back(numRemaining * 100./fNEvents);
    }
}

void
CutLog::LogPass(){
    fNEvents++;
}

void
CutLog::LogCut(size_t index_){
    if(index_ >= fCutCounts.size())
        throw NotFoundError(Formatter() << "CutLog::LogCut tried to log event as cut by non-existent cut #" << index_ );
    fCutCounts[index_]++;
    fNEvents++;
}

std::vector<int>
CutLog::GetCutCounts() const{
    return fCutCounts;
}

std::vector<double>
CutLog::GetCutPercentages() const{
    return fCutPercentages;
}

std::vector<std::string>
CutLog::GetCutNames() const{
    return fCutNames;
}

std::string
CutLog::AsString(){
    if (!fCutCounts.size())
        return "";

    if (fCutCounts.size() != fCutNames.size())
        throw LogicError(Formatter() << "CutLog::AsString() there are "
                         << fCutCounts.size()
                         << " cut counts and "
                         << fCutNames.size() << " names!");
    CalculateMeta();
    Formatter outString;
    outString << std::left
              << std::setw(15)
              << "Cut"
              << std::left
              << std::setw(15)
              << "# cut"
              << std::left
              << std::setw(9)
              << "(%)"
              << std::left
              << std::setw(15)
              << "# remaining"
              << std::left
              << std::setw(9)
              << "(%)";

    outString << "\n";
    for(size_t i = 0; i < fCutNames.size(); i++){
        outString << std::left
                  << std::setw(15)
                  << fCutNames.at(i)
                  << std::left
                  << std::setw(15)
                  << fCutCounts.at(i)
                  << std::left 
                  << std::setw(9)
                  << std::setprecision(3)
                  << fCutPercentages.at(i)                  
                  << std::left
                  << std::setw(15)
                  << fRemainderCounts.at(i)
                  << std::left
                  << std::setw(9)
                  << std::setprecision(3)
                  << fRemainderPercentages.at(i)
                  << "\n";


    }
    return outString;
}

void
CutLog::Print(){
    std::cout << AsString() << std::endl;
}

void
CutLog::SaveAs(const std::string& title_, 
	       const std::string& path_){
  std::ofstream fs;
  if (fs.is_open()){
    throw IOError("CutLog::SaveAs - couldn't open file" + path_);
  }
  fs.open(path_.c_str());
  fs << "-----------------------------------"
     << "-----------------------------------" 
     << std::endl
     << "Cut Log : " << title_  << "\n"
     << "-----------------------------------"
     << std::endl
     << AsString()
     << "-----------------------------------"
     << "-----------------------------------" 
     << std::endl;
  fs.close();  
}
