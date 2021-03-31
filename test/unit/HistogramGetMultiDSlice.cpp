#include <catch.hpp>
#include <Histogram.h>
#include <cstdlib>
#include <iostream>

TEST_CASE("Multi-Dimensional Histogram Slicing"){
    
    // create a 4D histogram from random data
    AxisCollection axes; 
    axes.AddAxis(BinAxis("ax1", 0, 6000, 100));
    axes.AddAxis(BinAxis("ax2", 5, 10, 10)); 
    axes.AddAxis(BinAxis("ax3", 0, 50, 100));
    axes.AddAxis(BinAxis("ax4", 0, 200, 50)); 

    Histogram origHistogram(axes); 
    // fill the histogram bins with values --> probably a horrible method! 
    for (size_t i = 0; i < origHistogram.GetNBins(); i++){
        origHistogram.SetBinContent(i, i);
    }

    // get the original histogram axis 
    BinAxis origAx1 = origHistogram.GetAxes().GetAxis(0);
    BinAxis origAx2 = origHistogram.GetAxes().GetAxis(1);
    BinAxis origAx3 = origHistogram.GetAxes().GetAxis(2);
    BinAxis origAx4 = origHistogram.GetAxes().GetAxis(3);

    SECTION("1D Slice"){
        // slice over 3 axes to yield 1D slice 
        std::map<std::string, size_t> fixedBins; 
        fixedBins["ax1"] = 1;
        fixedBins["ax2"] = 4;
        fixedBins["ax3"] = 10;
        
        Histogram dim1 = origHistogram.GetSlice(fixedBins);
        BinAxis newAx = dim1.GetAxes().GetAxis(0); 
        
        REQUIRE(dim1.GetNDims() == 1); // check it's a 1D slice
        REQUIRE(newAx.GetNBins() == origAx4.GetNBins()); // check sliced out axis is same as 
                                                         // in original histogram   
    }  

    SECTION("2D Slice"){
        // slice over 2 axes to yield 2D slice 
        std::map<std::string, size_t> fixedBins; 
        fixedBins["ax1"] = 1;
        fixedBins["ax2"] = 4;
        
        Histogram dim2 = origHistogram.GetSlice(fixedBins);
        BinAxis newAx1 = dim2.GetAxes().GetAxis(0);
        BinAxis newAx2 = dim2.GetAxes().GetAxis(1);  
        
        REQUIRE(dim2.GetNDims() == 2); // check it's a 2D slice
        REQUIRE(newAx1.GetNBins() == origAx3.GetNBins()); // check sliced out axes are same
        REQUIRE(newAx2.GetNBins() == origAx4.GetNBins()); // as in original histogram   
    }

    SECTION("3D Slice"){
        // slice over 1 axes to yield 3D slice 
        std::map<std::string, size_t> fixedBins; 
        fixedBins["ax1"] = 1;
        
        Histogram dim3 = origHistogram.GetSlice(fixedBins);
        BinAxis newAx1 = dim3.GetAxes().GetAxis(0);
        BinAxis newAx2 = dim3.GetAxes().GetAxis(1);  
        BinAxis newAx3 = dim3.GetAxes().GetAxis(2); 

        REQUIRE(dim3.GetNDims() == 3); // check it's a 3D slice
        REQUIRE(newAx1.GetNBins() == origAx2.GetNBins()); // check sliced out axes are same
        REQUIRE(newAx2.GetNBins() == origAx3.GetNBins()); // as in original histogram   
        REQUIRE(newAx3.GetNBins() == origAx4.GetNBins());
    }
}


