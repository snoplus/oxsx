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
#include <numeric>

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
DataSetGenerator::ExpectedRatesDataSet(std::vector<int>* eventsTaken_){
    if(fExpectedRates.size() != fDataSets.size())
        throw LogicError("Can't generate fake data: need one rate exactly for each data set");

    if(!fDataSets.size())
        throw LogicError("Can't generate fake data: no input data sets");
    
    // if requested note down the events taken
    if(eventsTaken_){
      eventsTaken_->clear();
    }

    OXSXDataSet dataSet;    
    dataSet.SetObservableNames(fDataSets.at(0)->GetObservableNames());
    for(size_t i = 0; i < fDataSets.size(); i++){
        unsigned expectedCounts = round(fExpectedRates.at(i));

	if(fBootstraps.size() == 0)
	  throw NotFoundError("Bootstrap flags not set for dataset");

	if(fBootstraps.at(i))
	  RandomDrawsWithReplacement(i, expectedCounts, dataSet);
	else{

	  if(fSequentialFlags.size() == 0)
	    throw NotFoundError("Sequential flags not set for dataset");

	  if(fSequentialFlags.at(i))
	    SequentialDrawsNoReplacement(i, expectedCounts, dataSet);
	  else
	    RandomDrawsNoReplacement(i, expectedCounts, dataSet);
	}
	if(eventsTaken_)
	  eventsTaken_->push_back(expectedCounts);
    }
            
    return dataSet;
}

OXSXDataSet
DataSetGenerator::PoissonFluctuatedDataSet(std::vector<int>* eventsTaken_){
    if(fExpectedRates.size() != fDataSets.size())
        throw LogicError("Can't generate fake data: need one rate exactly for each data set");
    if(!fDataSets.size())
      throw LogicError("Can't generate fake data: no input data sets");

    // if requested note down the events taken
    if(eventsTaken_){
      eventsTaken_->clear();
    }

    // work out how much memory you'll probably need
    int eventsNeeded = std::accumulate(fExpectedRates.begin(), fExpectedRates.end(), 0.0) * 1.2;
    
    OXSXDataSet dataSet;    
    dataSet.SetObservableNames(fDataSets.at(0)->GetObservableNames());
    dataSet.Reserve(eventsNeeded);

    for(size_t i = 0; i < fDataSets.size(); i++){
        size_t counts = Rand::Poisson(fExpectedRates.at(i));

	if(fBootstraps.size() == 0)
          throw NotFoundError("Bootstrap flags not set for dataset");

        if(fBootstraps.at(i))
            RandomDrawsWithReplacement(i, counts, dataSet);

        else{

	  if(fSequentialFlags.size() == 0)
            throw NotFoundError("Sequential flags not set for dataset");

	  if(fSequentialFlags.at(i))
	    SequentialDrawsNoReplacement(i, counts, dataSet);

	  else
	    RandomDrawsNoReplacement(i, counts, dataSet);
        }
        
        if(eventsTaken_)
            eventsTaken_->push_back(counts);
    }        
    return dataSet;
}


OXSXDataSet
DataSetGenerator::AllValidEvents(std::vector<int>* eventsTaken_){
    // if requested note down the events taken
    if(eventsTaken_){
      eventsTaken_->clear();
    }

    OXSXDataSet dataSet;
    dataSet.SetObservableNames(fDataSets.at(0)->GetObservableNames());
    for(size_t i = 0; i < fDataSets.size(); i++){
      for(size_t j = 0; j < fDataSets.at(i)->GetNEntries(); j++){
        Event event_ = fDataSets.at(i)->GetEntry(j);
        dataSet.AddEntry(event_);
      } // events

      if(eventsTaken_)
	eventsTaken_->push_back(fDataSets.at(i)->GetNEntries());

    } // data sets
    return dataSet;
}


// std::vector<OXSXDataSet*> 
// DataSetGenerator::AllRemainingEvents(std::vector<int>* eventsTaken_){
//   if(eventsTaken_)
//     eventsTaken_->clear();

//   std::vector<OXSXDataSet*> remainders(fDataSets.size(), NULL);
//   int taken;
//   for(size_t iDS = 0; iDS < fDataSets.size(); iDS++){   
//     remainders[iDS] = AllRemainingEvents(iDS, *taken);
//     eventsTaken_->push_back(taken);
//   }

//   return remainders;
// }

void
DataSetGenerator::Reset(){
  fEventIndicies.clear();
  fEventIndicies.resize(fDataSets.size());
  fMaxs = std::vector<size_t>(fDataSets.size(), -999);
}

void
DataSetGenerator::RandomDrawsNoReplacement(size_t handleIndex_, size_t nEvents_,
                                           OXSXDataSet& outData_
                                           ){

  // see http://stackoverflow.com/questions/196017/unique-non-repeating-random-numbers-in-o1#196065
  // credit to Robert Gamble http://stackoverflow.com/users/25222/robert-gamble
  std::vector<size_t>& eventIndices = fEventIndicies[handleIndex_];
  size_t& max = fMaxs[handleIndex_];
  DataSet* origData = fDataSets.at(handleIndex_);

  if(!eventIndices.size()){
    eventIndices.reserve(origData->GetNEntries());
    for(size_t i = 0; i < origData -> GetNEntries(); i++)
      eventIndices.push_back(i);
    max = eventIndices.size() - 1; // the effective end of the array
  }

  if(nEvents_ > origData -> GetNEntries())
    throw NotFoundError(Formatter() << "DataSetGenerator::RandomDrawsNoReplacement() asked for "
                        << nEvents_ << " but only have " << origData -> GetNEntries() 
                        << "events!)");

  size_t cache = -1; // for swapping
  size_t draw  = -1; // the random draw 

  outData_.Reserve(outData_.GetNEntries() + nEvents_);
  for(size_t i = 0; i < nEvents_; i++){
    if(!(i%10000))
      std::cout << i << "/" << nEvents_ << std::endl;

    if (max==static_cast<size_t>(-999))
      max = eventIndices.size() -1;
    // draw
    draw = Rand::Shoot(max);

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
DataSetGenerator::RandomDrawsWithReplacement(size_t handleIndex_, size_t nEvents_, 
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
  for(size_t i = 0; i < nEvents_; i++){
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
DataSetGenerator::AddDataSet(DataSet* data_, double rate_, bool flag_){
    fSequentialFlags.push_back(flag_);
    fDataSets.push_back(data_);
    fExpectedRates.push_back(rate_);
    fEventIndicies.push_back(std::vector<size_t>());
    fMaxs.push_back(0);
}


const std::vector<bool>&
DataSetGenerator::GetBootstrap() const{
    return fBootstraps;
}

void
DataSetGenerator::SetBootstrap(const std::vector<bool>& b_){
    fBootstraps = b_;
}

void 
DataSetGenerator::ClearDataSets(){
    fSequentialFlags.clear();
    fExpectedRates.clear();
    fEventIndicies.clear();
    fMaxs.clear();
    fDataSets.clear();
}

void
DataSetGenerator::SequentialDrawsNoReplacement(size_t handleIndex_, size_t nEvents_, OXSXDataSet& data_){
  std::vector<size_t>& eventIndices = fEventIndicies[handleIndex_];
  size_t& max = fMaxs[handleIndex_];
  DataSet* origData = fDataSets.at(handleIndex_);

  if(nEvents_ > origData->GetNEntries())
      throw NotFoundError(Formatter() << "DataSetGenerator::Not enough events!"
                          << "\n\t (asked for " << nEvents_
                          << " but only have " << max << " (left)"
                          );


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
  size_t cache = -1;
  //data_.Reserve(data_.GetNEntries() + nEvents_);

  size_t oneTenth = nEvents_/10;
  for(size_t i = 0; i < nEvents_; i++){
    if(oneTenth && !(i % oneTenth))
      std::cout << i << "/" << nEvents_ << "   (  " << 10 * i /oneTenth << " %)" << std::endl;
        
    if (max==static_cast<size_t>(-999))
        max = eventIndices.size() -1;

    // draw --> always take the first one 
    draw = 0;
    
    // swap
    cache = eventIndices[draw];
    eventIndices[draw] = eventIndices[max];
    eventIndices[max]  = cache;

    // return 
    if(fCuts.PassesCuts(origData -> GetEntry(eventIndices[draw])))
      data_.AddEntry(origData -> GetEntry(eventIndices[draw]));
    
    // deincrement
    max--;
  }
  
  return;  
}

void
DataSetGenerator::SetSequentialFlags(const std::vector<bool>& flags_){
  fSequentialFlags = flags_;
}


const std::vector<bool>&
DataSetGenerator::GetSequentialFlags() const{
  return fSequentialFlags;
}


OXSXDataSet*
DataSetGenerator::AllRemainingEvents(size_t iDS_, int* eventsTaken_){
  OXSXDataSet* ds = new OXSXDataSet();
  ds -> SetObservableNames(fDataSets.at(iDS_)->GetObservableNames());
  *eventsTaken_ = 0;

  size_t nEvents = fMaxs.at(iDS_);
  if (nEvents==static_cast<size_t>(-999)) 
    nEvents = (fDataSets.at(iDS_) -> GetNEntries()) - 1;
  
  if (nEvents==static_cast<size_t>(-1))
    return ds; //in case all events drawn and not replaced
  
  const std::vector<size_t>& eventIndices = fEventIndicies[iDS_];
  
  size_t oneTenth = nEvents/10;
  
  for(size_t jEV = 0; jEV <= nEvents; jEV++){   
    if(oneTenth && !(jEV % oneTenth))
      std::cout << jEV << "/" << nEvents << "\t" 
		<< "(" <<  10 * jEV/oneTenth
		<< "%)"
		<< std::endl;
    ds -> AddEntry(fDataSets.at(iDS_)->GetEntry(eventIndices.at(jEV)));
  }
  
  *eventsTaken_ = nEvents;
  return ds;
}
