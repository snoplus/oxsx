#include <catch2/catch_approx.hpp>
#include <catch2/catch_all.hpp>
#include <Event.h>
#include <EventScale.h>

TEST_CASE("Simple scale on variable number 3"){
    EventScale scaler("scaler");
    std::vector<std::string> observableScaler;
    observableScaler.push_back("obs3");
    scaler.SetOutObservables(observableScaler);
    scaler.SetScale(2);

    // make fake event
    std::vector<std::string> observablesData;
    observablesData.push_back("obs0");
    observablesData.push_back("obs1");
    observablesData.push_back("obs2");
    observablesData.push_back("obs3");

    std::vector<double> fakeObservations;
    fakeObservations.push_back(0);
    fakeObservations.push_back(1);
    fakeObservations.push_back(2);
    fakeObservations.push_back(3);

    Event inData(fakeObservations);
    inData.SetObservableNames(&observablesData);
    
    // apply the scale
    Event outData = scaler(inData);
    std::vector<double> modifiedObs = outData.GetData();

    // Check the values
    std::vector<double> correctVals = fakeObservations;
    correctVals[3] *= 2;

    REQUIRE(modifiedObs == correctVals);
}
