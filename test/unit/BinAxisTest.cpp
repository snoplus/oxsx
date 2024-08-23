#include <catch2/catch_all.hpp>
#include <BinAxis.h>
#include <iostream>

TEST_CASE("Bin Finding For Uniform Binning 1000 bins 0-100", "[BinAxis]"){
    BinAxis uniform("uniform axis", 0 , 100, 1000);
    
    SECTION("Last bin behaviour"){
        REQUIRE(uniform.FindBin(100) == 999);
        REQUIRE(uniform.FindBin(102) == 999);
        REQUIRE(uniform.GetNBins() == 1000);
    }


    SECTION("First bin behaviour"){
        REQUIRE(uniform.FindBin(0)   == 0);
        REQUIRE(uniform.FindBin(-20) == 0);
        REQUIRE(uniform.GetNBins()   == 1000);
    }

}

TEST_CASE("Bin Finding For Custom Binning 1000 bins 0-1001", "[BinAxis]"){

    std::vector<double> lowEdges(1000);
    std::vector<double> highEdges(1000);

    for(size_t i = 0; i < lowEdges.size(); i++){
        lowEdges[i] = i;
        highEdges[i] = i+1;
    }

    BinAxis custom("custombins", lowEdges, highEdges);
    
    SECTION("Last bin behaviour"){
        REQUIRE(custom.FindBin(1000) == 999);
        REQUIRE(custom.FindBin(1020) == 999);
        REQUIRE(custom.GetNBins() == 1000);
    }


    SECTION("First bin behaviour"){
        REQUIRE(custom.FindBin(0)   == 0);
        REQUIRE(custom.FindBin(-20) == 0);
			 
        REQUIRE(custom.GetNBins()   == 1000);

    }
}
