#include <catch2/catch_all.hpp>
#include <AxisCollection.h>

TEST_CASE("Correct Indexing for 3x3 Uniform Axes","[AxisCollection]"){

    BinAxis axis1("axis1", 0,  20, 100);
    BinAxis axis2("axis2", 10, 30, 219);
    BinAxis axis3("axis3", 0, 1, 2341);

    AxisCollection axes;
    axes.AddAxis(axis1);
    axes.AddAxis(axis2);
    axes.AddAxis(axis3);

    SECTION("Correct Number of Bins"){
        REQUIRE(axes.GetNBins() == (100 * 219 * 2341));
    }

    SECTION("Correct Number of Dimensions"){
        REQUIRE(axes.GetNDimensions() == 3);
    }


    SECTION("Global Bin Conversions for first and last bin"){
        std::vector<size_t> firstBins(3, 0);
        std::vector<size_t> lastBins;
        lastBins.push_back(99);
        lastBins.push_back(218);
        lastBins.push_back(2340);

        REQUIRE(axes.FlattenIndices(firstBins) == 0);
        REQUIRE(axes.FlattenIndices(lastBins)  ==  (2341 * 100 * 219 - 1));
    }

    SECTION("Global Bin Self Conversions"){
        // pick a random bin
        for(unsigned i = 0; i < 100; i++ ){
            size_t bin = rand()%(2341 * 219 * 100);

            std::vector<size_t> indices = axes.UnpackIndices(bin);
            size_t bin2 = axes.FlattenIndices(indices);

            REQUIRE(bin == bin2);
        }
    }

}

TEST_CASE("Correct Indexing when there's just one axis in there", "[AxisCollection]"){
    AxisCollection axes;
    BinAxis axis("test", -20, 55, 9875);
    axes.AddAxis(axis);

    SECTION("Right bin numbers and indices"){
        for(size_t i = 0; i < 100; i++){
            size_t axisBin  = rand() % axis.GetNBins();
            std::vector<size_t> axesBins = axes.UnpackIndices(axisBin);

            REQUIRE(axesBins.size() == 1);
            REQUIRE(axesBins.at(0)  == axisBin);
        }
    }
}


TEST_CASE("Bin Edges for 4x4 binning"){
    BinAxis axis1("axis1", 0,  20, 100);
    BinAxis axis2("axis2", 10, 30, 219);
    BinAxis axis3("axis3", 0, 1, 2341);
    BinAxis axis4("axis4", 0, 2, 2241);

    AxisCollection axes;
    axes.AddAxis(axis1);
    axes.AddAxis(axis2);
    axes.AddAxis(axis3);
    axes.AddAxis(axis4);

    SECTION("Low edge < centre < high edge whole"){
        std::vector<double> centres(4);
        std::vector<double> lows(4);
        std::vector<double> highs(4);

        for(size_t i = 0; i < 100; i++){
            size_t bin = rand() % axes.GetNBins();
            axes.GetBinCentres(bin, centres);
            axes.GetBinLowEdges(bin, lows);
            axes.GetBinHighEdges(bin, highs);
        }
        
        for(size_t i = 0; i < lows.size(); i++){
            REQUIRE(lows.at(i) < highs.at(i));
            REQUIRE(lows.at(i) < centres.at(i));
            REQUIRE(centres.at(i) < highs.at(i));
        }
    }


    SECTION("Low edge < centre < high edge whole"){
        std::vector<double> centres(4);
        std::vector<double> lows(4);
        std::vector<double> highs(4);

        for(size_t i = 0; i < 100; i++){
            size_t bin = rand() % axes.GetNBins();
            axes.GetBinCentres(bin, centres);
            axes.GetBinLowEdges(bin, lows);
            axes.GetBinHighEdges(bin, highs);
        }
        
        for(size_t i = 0; i < lows.size(); i++){
            REQUIRE(lows.at(i) < highs.at(i));
            REQUIRE(lows.at(i) < centres.at(i));
            REQUIRE(centres.at(i) < highs.at(i));
        }
    }    
}

