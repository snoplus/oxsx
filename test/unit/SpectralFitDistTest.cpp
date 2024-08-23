#include <catch2/catch_all.hpp>
#include <Formatter.hpp>
#include <ContainerTools.hpp>
#include <SpectralFitDist.h>

TEST_CASE(""){    
    AxisCollection axes;
    axes.AddAxis(BinAxis("test1", 0, 10, 10));
    axes.AddAxis(BinAxis("test2", 0, 10, 10));

    SpectralFitDist fitED("test", axes);
    
    SECTION("Behaves like a binned pdf"){
        REQUIRE(fitED.GetNDims() == 2);
        REQUIRE(fitED.GetNBins() == 100);
        REQUIRE(fitED.GetBinContents() == std::vector<double>(100, 0));
    }
    SECTION("Bin contents are adjustable"){
        REQUIRE(fitED.GetParameterCount() == 100);
        ParameterDict testP;
        for(size_t i = 0; i < fitED.GetNBins(); i++)
            testP[Formatter() << "test_bin_" << i] = 0;

        REQUIRE(fitED.GetParameters() == testP);


        for(size_t i = 0; i < fitED.GetNBins(); i++)
            testP[Formatter() << "test_bin_" << i] = 0;

        fitED.SetParameters(testP);
        REQUIRE(fitED.GetParameters() == testP);
        REQUIRE(fitED.GetBinContents() == ContainerTools::GetValues(testP, fitED.GetParameterNames()));
    }
    
}
