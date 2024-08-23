#include <catch2/catch_all.hpp>
#include <Event.h>
#include <EventShift.h>

TEST_CASE("Simple shift on variable number 2"){
    EventShift shifter("shifter");
    std::vector<std::string> observableShifter;
    observableShifter.push_back("obs2");
    shifter.SetOutObservables(observableShifter);
    shifter.SetShift(1);

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

    // apply the shift
    Event outData = shifter(inData);
    std::vector<double> modifiedObs = outData.GetData();

    // Check the values
    std::vector<double> correctVals = fakeObservations;
    correctVals[2] += 1;

    REQUIRE(modifiedObs == correctVals);
}
