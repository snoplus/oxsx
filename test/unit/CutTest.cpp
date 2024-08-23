#include <catch2/catch_all.hpp>
#include <BoxCut.h>
#include <BoolCut.h>
#include <LineCut.h>
#include <CutCollection.h>

TEST_CASE("1D Box cut"){
    BoxCut cut("boxCut","obs0", 1 , 10);
    cut.SetLowerLimit(1);
    cut.SetUpperLimit(10);
    
    std::vector<std::string> observable;
    observable.push_back("obs0");

    SECTION("in the middle of accepted range"){
        Event ev(std::vector<double>(1, 5));
        ev.SetObservableNames(&observable);
        REQUIRE(cut.PassesCut(ev) == true);
    }

    SECTION("below accepted range"){
        Event ev(std::vector<double>(1, 0));
        ev.SetObservableNames(&observable);
        REQUIRE(cut.PassesCut(ev) == false);
    }

    SECTION("above accepted range"){
        Event ev(std::vector<double>(1, 11));
        ev.SetObservableNames(&observable);
        REQUIRE(cut.PassesCut(ev) == false);
    }

    SECTION("on lower boundary"){
        Event ev(std::vector<double>(1, 1));
        ev.SetObservableNames(&observable);
        REQUIRE(cut.PassesCut(ev) == false);
    }

    SECTION("on upper boundary"){
        Event ev(std::vector<double>(1, 10));
        ev.SetObservableNames(&observable);
        REQUIRE(cut.PassesCut(ev) == false);
    }
}

TEST_CASE("Several Cuts on MultiD data"){
    //  First val must be more than 2
    LineCut lineCut("line", "obs0", 2, "lower");

    // cut if second observable is between 5 and 10
    BoxCut boxCut("box", "obs1", 5, 10);

    // cut if the third value is 1 exactly
    BoolCut boolCut("bool","obs2", 1);

    // Combine the three cuts
    CutCollection combinedCuts;
    combinedCuts.AddCut(lineCut);
    combinedCuts.AddCut(boxCut);
    combinedCuts.AddCut(boolCut);

    // Create a fake event to test
    std::vector<double> observations;
    observations.push_back(1);
    observations.push_back(6);
    observations.push_back(1);
    Event fakeEvent(observations);
    std::vector<std::string> observablesEvent;
    observablesEvent.push_back("obs0");
    observablesEvent.push_back("obs1");
    observablesEvent.push_back("obs2");
    fakeEvent.SetObservableNames(&observablesEvent);

    // Check vals
    REQUIRE(lineCut.PassesCut(fakeEvent) == false);
    REQUIRE(boxCut.PassesCut(fakeEvent)  == true);
    REQUIRE(boolCut.PassesCut(fakeEvent) == true);
    REQUIRE(combinedCuts.PassesCuts(fakeEvent) == false);
}
