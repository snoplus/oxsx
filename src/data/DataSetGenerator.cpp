#include <DataSetGenerator.h>
#include <DataSet.h>
#include <OXSXDataSet.h>
#include <Exceptions.h>
#include <Event.h>
#include <Rand.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <iostream>

void
DataSetGenerator::SetCuts(const CutCollection& cuts_){
    fCuts=cuts_;
}

void
DataSetGenerator::AddCut(const Cut& cut_){
    fCuts.AddCut(cut_);
}

void
DataSetGenerator::SetDataSets(const std::vector<DataSet*> sets_){
    fDataSets = sets_;
    fEventIndicies.clear();
    fMaxs.clear();
    fEventIndicies.resize(sets_.size());
    fMaxs.resize(sets_.size(), -999);
}

void
DataSetGenerator::SetExpectedRates(const std::vector<double>& rates_){
    fExpectedRates = rates_;
}

OXSXDataSet
DataSetGenerator::ExpectedRatesDataSet(){
    if(fExpectedRates.size() != fDataSets.size())
        throw LogicError("Can't generate fake data: need one rate exactly for each data set");

    if(!fDataSets.size())
        throw LogicError("Can't generate fake data: no input data sets");

    OXSXDataSet dataSet;    
    dataSet.SetObservableNames(fDataSets.at(0)->GetObservableNames());
    for(size_t i = 0; i < fDataSets.size(); i++){
        unsigned expectedCounts = round(fExpectedRates.at(i));
	fBootstrap ? RandomDrawsWithReplacement(i, expectedCounts, dataSet) 
	  : RandomDrawsNoReplacement(i, expectedCounts, dataSet);
    }
            
    return dataSet;
}

OXSXDataSet
DataSetGenerator::PoissonFluctuatedDataSet(){
    if(fExpectedRates.size() != fDataSets.size())
        throw LogicError("Can't generate fake data: need one rate exactly for each data set");
    if(!fDataSets.size())
      throw LogicError("Can't generate fake data: no input data sets");

    OXSXDataSet dataSet;    
    dataSet.SetObservableNames(fDataSets.at(0)->GetObservableNames());
    for(size_t i = 0; i < fDataSets.size(); i++){
        int counts = Rand::Poisson(fExpectedRates.at(i));
        fBootstrap ? RandomDrawsWithReplacement(i, counts, dataSet) 
	  : RandomDrawsNoReplacement(i, counts, dataSet);
    }        
    return dataSet;
}


OXSXDataSet
DataSetGenerator::AllValidEvents(){
    OXSXDataSet dataSet;
    dataSet.SetObservableNames(fDataSets.at(0)->GetObservableNames());
    for(size_t i = 0; i < fDataSets.size(); i++){
      for(size_t j = 0; j < fDataSets.at(i)->GetNEntries(); j++){
        Event event_ = fDataSets.at(i)->GetEntry(j);
        dataSet.AddEntry(event_);
      } // events
    } // data sets
    return dataSet;
}


std::vector<OXSXDataSet*> 
DataSetGenerator::AllRemainingEvents(){
  std::vector<OXSXDataSet*> remainders(fDataSets.size(), NULL);
  for(size_t iDS = 0; iDS < fDataSets.size(); iDS++){   
    remainders[iDS] = new OXSXDataSet();
    remainders[iDS]->SetObservableNames(fDataSets.at(0)->GetObservableNames());
  }
  for(size_t iDS = 0; iDS < fDataSets.size(); iDS++){   
    size_t nEvents = fMaxs.at(iDS);
    if (nEvents==static_cast<size_t>(-999)) 
      nEvents = (fDataSets.at(iDS) -> GetNEntries()) - 1;

    if (nEvents==static_cast<size_t>(-1)) continue; //in case all events drawn and not replaced
    const std::vector<size_t>& eventIndices = fEventIndicies[iDS];

    for(size_t jEV = 0; jEV <= nEvents; jEV++){   
      remainders[iDS]->AddEntry(fDataSets.at(iDS)->GetEntry(eventIndices.at(jEV)));
    }
  }
  Reset();
  return remainders;
}

void
DataSetGenerator::Reset(){
  fEventIndicies.clear();
  fEventIndicies.resize(fDataSets.size());
  fMaxs = std::vector<size_t>(fDataSets.size(), -999);
}

void
DataSetGenerator::RandomDrawsNoReplacement(size_t handleIndex_, int nEvents_,
                                           OXSXDataSet& outData_
                                           ){

  // see http://stackoverflow.com/questions/196017/unique-non-repeating-random-numbers-in-o1#196065
  // credit to Robert Gamble http://stackoverflow.com/users/25222/robert-gamble
  std::vector<size_t>& eventIndices = fEventIndicies[handleIndex_];
  size_t& max = fMaxs[handleIndex_];
  DataSet* origData = fDataSets.at(handleIndex_);
  if(static_cast<int>(origData->GetNEntries()) < nEvents_)
    throw NotFoundError(Formatter() << "DataSetGenerator::RandomDrawsNoReplacement() asked for "
                        << nEvents_ << " but only have " << origData -> GetNEntries() 
                        << "events!");

  if(!eventIndices.size()){
    eventIndices.reserve(origData->GetNEntries());
    for(size_t i = 0; i < origData -> GetNEntries(); i++)
      eventIndices.push_back(i);
    max = eventIndices.size() - 1; // the effective end of the array
  }

  size_t cache = -1; // for swapping
  size_t draw  = -1; // the random draw 

  outData_.Reserve(outData_.GetNEntries() + nEvents_);
  for(int i = 0; i < nEvents_; i++){
    if(!(i%10000))
      std::cout << i << "/" << nEvents_ << std::endl;

    if (max==static_cast<size_t>(-999))
      max = eventIndices.size() -1;
    // draw
    draw = Rand::Shoot(max);

    // swap
    cache = eventIndices[draw];
    eventIndices[draw] = eventIndices[max];
    eventIndices[max]  = cache;

    // return 
    if(fCuts.PassesCuts(origData -> GetEntry(eventIndices[max])))
        outData_.AddEntry(origData -> GetEntry(eventIndices[max]));

    // deincrement
    max--;
  }

  return;
}


void
DataSetGenerator::RandomDrawsWithReplacement(size_t handleIndex_, int nEvents_, 
					          OXSXDataSet& outData_
					     ){

  std::vector<size_t>& eventIndices = fEventIndicies[handleIndex_];
  size_t& max = fMaxs[handleIndex_];
  DataSet* origData = fDataSets.at(handleIndex_);
  if(!(origData->GetNEntries()))
    throw NotFoundError(Formatter() << "DataSetGenerator::RandomDrawsWithReplacement() asked for "
                        << nEvents_ << " but there are no events!");

  if(!eventIndices.size()){
    eventIndices.reserve(origData->GetNEntries());
    for(size_t i = 0; i < origData -> GetNEntries(); i++)
      eventIndices.push_back(i);
    max = eventIndices.size() - 1; // the effective end of the array
  }

  size_t draw  = -1; // the random draw 

  outData_.Reserve(outData_.GetNEntries() + nEvents_);
  for(int i = 0; i < nEvents_; i++){
    if(!(i%10000))
      std::cout << i << "/" << nEvents_ << std::endl;


    if (max==static_cast<size_t>(-999))
      max = eventIndices.size() -1;

    // draw
    draw = Rand::Shoot(max);

    // return 
    if(fCuts.PassesCuts(origData -> GetEntry(eventIndices[draw])))
      outData_.AddEntry(origData -> GetEntry(eventIndices[draw]));
  }

  return;
}


void
DataSetGenerator::AddDataSet(DataSet* data_, double rate_){
    fDataSets.push_back(data_);
    fExpectedRates.push_back(rate_);
    fEventIndicies.push_back(std::vector<size_t>());
    fMaxs.push_back(0);
}


bool
DataSetGenerator::GetBootstrap() const{
    return fBootstrap;
}

void
DataSetGenerator::SetBootstrap(bool b_){
    fBootstrap = b_;
}

void 
DataSetGenerator::ClearDataSets(){
    fExpectedRates.clear();
    fEventIndicies.clear();
    fMaxs.clear();
    fDataSets.clear();
}
