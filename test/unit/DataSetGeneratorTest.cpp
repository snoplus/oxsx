#include <catch.hpp>
#include <OXSXDataSet.h>
#include <DataSetGenerator.h>
#include <iostream>
#include <typeinfo>
#include <numeric>

TEST_CASE("Dataset generation by random draws"){

  //create events in 2 fake datasets
  OXSXDataSet testDataSet1;
  OXSXDataSet testDataSet2;
  size_t eventsInDataSets = 3;
  std::vector<double> eventObs1(4);
  std::vector<double> eventObs2(4);

  for(unsigned i = 0; i < eventsInDataSets; i++){
    for(size_t j = 0; j < eventObs1.size(); j++){
      eventObs1[j] = j+eventsInDataSets*i +1;
      eventObs2[j] = (-1)*eventObs1[j];
    }
    testDataSet1.AddEntry(Event(eventObs1));
    testDataSet2.AddEntry(Event(eventObs2));
  }

  //create dataset generator
  DataSetGenerator gen;
  std::vector<DataSet*> sets;
  sets.push_back(&testDataSet1);
  sets.push_back(&testDataSet2);
  gen.SetDataSets(sets);
  
  std::vector<bool> sequentialFlags(2,0);
  gen.SetSequentialFlags(sequentialFlags);

  std::vector<double> rates;
  rates.push_back(2.0);
  rates.push_back(3.0);
  gen.SetExpectedRates(rates);

  SECTION("Correct number of events drawn without replacement"){
    std::vector<bool> bootstraps;
    bootstraps.push_back(false);
    bootstraps.push_back(false);
    gen.SetBootstrap(bootstraps);
    OXSXDataSet newDataSet = gen.ExpectedRatesDataSet();
    REQUIRE(newDataSet.GetNEntries() == std::accumulate(rates.begin(), rates.end(), 0));
  }
  

  SECTION("Correct number of events left in datasets, no replacement"){
    std::vector<bool> bootstraps;
    bootstraps.push_back(false);
    bootstraps.push_back(false);
    gen.SetBootstrap(bootstraps);
    OXSXDataSet newDataSet = gen.ExpectedRatesDataSet();

    //get number of events remaining in original datasets
    size_t remainingEvents = 0;
    for (size_t i = 0; i<sets.size();i++){
      int countsTaken = 0;
      OXSXDataSet* remainder = gen.AllRemainingEvents(i,&countsTaken);
      remainingEvents+= remainder->GetNEntries();
    }

    REQUIRE(remainingEvents == (2*eventsInDataSets - std::accumulate(rates.begin(), rates.end(), 0)));
  }
  

  SECTION("Correct number of events drawn with replacement"){
    std::vector<bool> bootstraps;
    bootstraps.push_back(true);
    bootstraps.push_back(true);
    gen.SetBootstrap(bootstraps);
    OXSXDataSet newDataSet = gen.ExpectedRatesDataSet();
    REQUIRE(newDataSet.GetNEntries() == std::accumulate(rates.begin(), rates.end(), 0));
  }


  SECTION("Correct number of events left in datasets, with replacement"){
    std::vector<bool> bootstraps;
    bootstraps.push_back(true);
    bootstraps.push_back(true);
    gen.SetBootstrap(bootstraps);
    OXSXDataSet newDataSet = gen.ExpectedRatesDataSet();

    //get number of events remaining in original datasets
    size_t remainingEvents = 0;
    for (size_t i = 0; i<sets.size();i++){
      int countsTaken = 0;
      OXSXDataSet* remainder = gen.AllRemainingEvents(i, &countsTaken);
      remainingEvents+= remainder->GetNEntries();
    }

    REQUIRE(remainingEvents == 2*eventsInDataSets);
  }
}
