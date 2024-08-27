#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>
#include <ComponentManager.h>
#include <SpectralFitDist.h>
#include <Formatter.hpp>
#include <iostream>

TEST_CASE("Stand alone component manager")
{
    // some fake objects
    AxisCollection axes;
    axes.AddAxis(BinAxis("", 0, 10, 5));
    SpectralFitDist pdf1("test1", axes);
    SpectralFitDist pdf2("test2", axes);

    ComponentManager cmpMan;

    SECTION("initialised correctly")
    {
        REQUIRE(cmpMan.GetTotalParameterCount() == 0);
        REQUIRE(cmpMan.GetParameters() == ParameterDict());
        REQUIRE(cmpMan.GetParameterNames() == std::set<std::string>());
    }

    SECTION("adding a component")
    {
        cmpMan.AddComponent(&pdf1);
        REQUIRE(cmpMan.GetTotalParameterCount() == 5);
    }

    SECTION("adding two")
    {
        cmpMan.AddComponent(&pdf1);
        cmpMan.AddComponent(&pdf2);
        REQUIRE(cmpMan.GetTotalParameterCount() == 5 * 2);

        ParameterDict params1;
        for (int i = 0; i < 5; i++)
            params1[Formatter() << "test1_bin_" << i] = 8;

        ParameterDict params2;
        for (int i = 0; i < 5; i++)
            params2[Formatter() << "test2_bin_" << i] = 8;

        ParameterDict params = params1;
        params.insert(params2.begin(), params2.end());

        cmpMan.SetParameters(params);

        REQUIRE(pdf1.GetParameters() == params1);
        REQUIRE(pdf2.GetParameters() == params2);
    }
    SECTION("clearing")
    {
        REQUIRE(cmpMan.GetTotalParameterCount() == 0);
        REQUIRE(cmpMan.GetParameters() == ParameterDict());
        REQUIRE(cmpMan.GetParameterNames() == std::set<std::string>());
    }
    SECTION("getting parameter by name")
    {
        cmpMan.AddComponent(&pdf1);
        REQUIRE(cmpMan.GetParameter("test1_bin_1") == Catch::Approx(0.));
        pdf1.SetBinContent(0, 10);
        REQUIRE(cmpMan.GetParameter("test1_bin_0") == Catch::Approx(10.));
    }
}
