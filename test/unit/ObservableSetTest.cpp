#include <catch2/catch_all.hpp>
#include <ObsSet.h>
#include <OXSXDataSet.h>

TEST_CASE("Build from one index")
{
    ObsSet drep("obs19");
    REQUIRE(drep.GetNObservables() == 1);
    REQUIRE(drep.GetIndex("obs19") == 0);
}

TEST_CASE("Build from vec")
{
    std::vector<std::string> obsNames;
    obsNames.push_back("obs5");
    obsNames.push_back("obs9");
    obsNames.push_back("obs11");

    ObsSet drep(obsNames);
    REQUIRE(drep.GetNObservables() == 3);
    REQUIRE(drep.GetIndex("obs5") == 0);
    REQUIRE(drep.GetIndex("obs9") == 1);
    REQUIRE(drep.GetIndex("obs11") == 2);

    SECTION("Then copy construct")
    {
        ObsSet drep2(drep);
        REQUIRE(drep2.GetNObservables() == 3);
        REQUIRE(drep2.GetIndex("obs5") == 0);
        REQUIRE(drep2.GetIndex("obs9") == 1);
        REQUIRE(drep2.GetIndex("obs11") == 2);
    }

    SECTION("Testing Inverse Functionality")
    {
        REQUIRE(drep.GetNames().at(2) == "obs11");
        REQUIRE(drep.GetNames().at(0) == "obs5");
        REQUIRE(drep.GetNames().at(1) == "obs9");
    }
}

TEST_CASE("Looking for the relative obsNames of 2d rep in 4d rep")
{
    std::vector<std::string> smallRepObsNames;
    smallRepObsNames.push_back("obs6");
    smallRepObsNames.push_back("obs0");

    std::vector<std::string> bigRepObsNames;
    bigRepObsNames.push_back("obs0");
    bigRepObsNames.push_back("obs3");
    bigRepObsNames.push_back("obs6");

    ObsSet bigRep(bigRepObsNames);
    ObsSet smallRep(smallRepObsNames);

    std::vector<size_t> relativeIndices = smallRep.GetRelativeIndices(bigRep);

    REQUIRE(relativeIndices.size() == 2);
    REQUIRE(relativeIndices.at(0) == 2);
    REQUIRE(relativeIndices.at(1) == 0);
}
