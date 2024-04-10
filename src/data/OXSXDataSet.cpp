#include <OXSXDataSet.h>
#include <Event.h>
#include <Exceptions.h>
#include <iostream>
void
OXSXDataSet::AddEntry(const Event& evData_){
    // fNObservables is zero at initialisation, after the first time it is set it shouldnt change
    if(!fNObservables)
        fNObservables = evData_.GetNObservables();

    else if (fNObservables != evData_.GetNObservables())
        throw DimensionError("OXSXDataSet::AddEntry", fNObservables, 
                             evData_.GetNObservables(), 
                             " observables in event");

    if(!evData_.GetNObservables())
        throw DimensionError("OXSXDataSet::Tried to add empty event to data set!");
    
    fData.push_back(evData_);
    fData.back().SetObservableNames(&fObservableNames);
}

Event
OXSXDataSet::GetEntry(size_t eventIndex_) const{
    try{
        return fData.at(eventIndex_);
    }
    catch(std::out_of_range&){
        throw NotFoundError("");
    }
}

unsigned
OXSXDataSet::GetNEntries() const{
    return fData.size();
}

void
OXSXDataSet::SetObservableNames(const std::vector<std::string>& names_){
    fObservableNames = names_;
}

std::vector<std::string>
OXSXDataSet::GetObservableNames() const{
    return fObservableNames;
}


OXSXDataSet
OXSXDataSet::operator+ (const OXSXDataSet& other_){
    for(size_t i = 0; i < other_.GetNEntries(); i++)
        AddEntry(other_.GetEntry(i));
    
    return *this;
}

void
OXSXDataSet::Reserve(int i){
  fData.reserve(i);
}

unsigned
OXSXDataSet::GetNObservables() const{
    return fNObservables;
}
