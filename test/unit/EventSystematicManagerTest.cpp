#include <catch2/catch_all.hpp>
#include <EventSystematicManager.h>
#include <EventScale.h>
#include <EventShift.h>
#include <Event.h>
#include <iostream>

TEST_CASE("3 observable event, shift and scale")
{
    // Create the systematics
    EventScale scaler("scaler");
    scaler.SetScale(3);

    EventShift shifter("shifter");
    shifter.SetShift(1);

    // Create a fake event
    std::vector<double> data(3, 1);
    Event event(data);

    std::vector<std::string> obsData;
    obsData.push_back("obs0");
    obsData.push_back("obs1");
    obsData.push_back("obs2");
    event.SetObservableNames(&obsData);

    SECTION("Acting of different observables")
    {

        std::vector<std::string> obsScaler;
        obsScaler.push_back("obs0");

        std::vector<std::string> obsShifter;
        obsShifter.push_back("obs1");

        scaler.SetInObservables(obsScaler);
        scaler.SetOutObservables(obsScaler);

        shifter.SetInObservables(obsShifter);
        shifter.SetOutObservables(obsShifter);

        EventSystematicManager sysMan;
        sysMan.Add(&scaler);
        sysMan.Add(&shifter);

        // transform and test
        Event modified = sysMan.ApplySystematics(event);
        REQUIRE(modified.GetData().size() == 3);
        REQUIRE(modified.GetData().at(0) == (1 * 3));
        REQUIRE(modified.GetData().at(1) == (1 + 1));
        REQUIRE(modified.GetData().at(2) == 1);
    }

    SECTION("Acting on same observables")
    {
        std::vector<std::string> obsScalerAndShifter;
        obsScalerAndShifter.push_back("obs0");

        // Act in the order added to sysMan
        scaler.SetInObservables(obsScalerAndShifter);
        shifter.SetInObservables(obsScalerAndShifter);
        scaler.SetOutObservables(obsScalerAndShifter);
        shifter.SetOutObservables(obsScalerAndShifter);

        EventSystematicManager sysMan;
        sysMan.Add(&scaler);
        sysMan.Add(&shifter);

        // transform
        Event modified = sysMan.ApplySystematics(event);
        // Test
        REQUIRE(modified.GetData().size() == 3);
        REQUIRE(modified.GetData().at(0) == (1 * 3 + 1));
        REQUIRE(modified.GetData().at(1) == 1);
        REQUIRE(modified.GetData().at(2) == 1);
    }

    SECTION("Acting on same observables - otherway around")
    {
        // Act in the order added to sysMan
        std::vector<std::string> obsScalerAndShifter;
        obsScalerAndShifter.push_back("obs0");

        // Act in the order added to sysMan
        scaler.SetInObservables(obsScalerAndShifter);
        shifter.SetInObservables(obsScalerAndShifter);
        scaler.SetOutObservables(obsScalerAndShifter);
        shifter.SetOutObservables(obsScalerAndShifter);

        EventSystematicManager sysMan;
        sysMan.Add(&shifter);
        sysMan.Add(&scaler);

        // transform
        Event modified = sysMan.ApplySystematics(event);

        // Test
        REQUIRE(modified.GetData().at(0) == ((1 + 1) * 3));
    }
}
