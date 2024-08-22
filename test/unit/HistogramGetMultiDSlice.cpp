#include <catch2/catch_approx.hpp>
#include <catch2/catch_all.hpp>
#include <Histogram.h>
#include <cstdlib>
#include <iostream>

TEST_CASE("Multi-Dimensional Histogram Slicing"){
    
    // create a 4D histogram
    AxisCollection axes; 
    axes.AddAxis(BinAxis("ax1", 0, 10, 3));
    axes.AddAxis(BinAxis("ax2", 0, 10, 3)); 
    axes.AddAxis(BinAxis("ax3", 0, 10, 3));
    axes.AddAxis(BinAxis("ax4", 0, 10, 3)); 

    Histogram origHistogram(axes); 
    // fill the histogram bins with values --> probably a horrible method! 
    for (size_t i = 0; i < origHistogram.GetNBins(); i++){
        origHistogram.SetBinContent(i, i);
        //std::cout << "Filling bin " << i << "with value " << i << std::endl; 
    }

    // get the original histogram axis 
    BinAxis origAx1 = origHistogram.GetAxes().GetAxis(0);
    BinAxis origAx2 = origHistogram.GetAxes().GetAxis(1);
    BinAxis origAx3 = origHistogram.GetAxes().GetAxis(2);
    BinAxis origAx4 = origHistogram.GetAxes().GetAxis(3);

    SECTION("1D Slice"){
        // slice over 3 axes to yield 1D slice 
        std::map<std::string, size_t> fixedBins; 
        fixedBins["ax1"] = 2;
        fixedBins["ax2"] = 1;
        fixedBins["ax3"] = 0;
        Histogram dim1 = origHistogram.GetSlice(fixedBins);
        
        // grab the axis & total bin number in slice for later testing  
        BinAxis newAx = dim1.GetAxes().GetAxis(0); 
        size_t numBins = dim1.GetNBins(); 

        // CODE TO VERIFY BIN CONTENTS OF SLICE IS CORRECT 
        // vector stores the local idx for each bin in slice for each axis 
        std::vector<std::vector<size_t> > localIdx; 

        // loop over the "free index", ie the axis not present in fixedBins
        for(size_t ax4Bin = 0; ax4Bin < 3; ax4Bin++){
            std::vector<size_t> idx;
            // these are the constant bins - where the slice 'pathway' is defined 
            idx.push_back(2); // idx of bin in ax1 
            idx.push_back(1); // idx of bin in ax2 
            idx.push_back(0); // idx of bin in ax3 

            // everything in 4th axis is included - it is the free idx  
            idx.push_back(ax4Bin);
            localIdx.push_back(idx);   
        }
        
        // now sum the bin contents in original histogram corresponding to sliced out bins 
        double sumOrig = 0;
        double sumSlice = dim1.Integral();
        for(size_t i = 0; i < localIdx.size(); i++){
            // convert from local to global (flat) idx 
            size_t globalIdx = origHistogram.FlattenIndices(localIdx.at(i)); 
            std::cout << "globIdx = " << globalIdx << std::endl; 
            // get contents in origHistogram at that flattened idx and sum 
            sumOrig += origHistogram.GetBinContent(globalIdx); 
            std::cout << "Bin content = " << origHistogram.GetBinContent(globalIdx) << std::endl;  
        }
        std::cout << "new sum = " << sumSlice << std::endl; 
        std::cout << "orig sum = " << sumOrig << std::endl;  
        REQUIRE(dim1.GetNDims() == 1); // check it's a 1D slice
        REQUIRE(newAx.GetNBins() == origAx1.GetNBins()); // check sliced out axis is same as in original histogram  
        REQUIRE(sumSlice == sumOrig); // check total contents of slice is equal to sum of bin contents in 
                                      // sliced out original bins
    }  

    SECTION("2D Slice"){
        // slice over 2 axes to yield 2D slice 
        std::map<std::string, size_t> fixedBins; 
        fixedBins["ax1"] = 1;
        fixedBins["ax2"] = 2;
        
        Histogram dim2 = origHistogram.GetSlice(fixedBins);
        BinAxis newAx1 = dim2.GetAxes().GetAxis(0);
        BinAxis newAx2 = dim2.GetAxes().GetAxis(1);  

        // CODE TO VERIFY BIN CONTENTS OF SLICE IS CORRECT 
        // vector stores the local idx for each bin in slice for each axis 
        std::vector<std::vector<size_t> > localIdx; 

        // loop over the "free index", ie the axis not present in fixedBins
        for(size_t ax3Bin = 0; ax3Bin < 3; ax3Bin++){
            for(size_t ax4Bin = 0; ax4Bin < 3; ax4Bin++){
                std::vector<size_t> idx;
                // these are the constant bins - where the slice 'pathway' is defined 
                idx.push_back(1); // idx of bin in ax1 
                idx.push_back(2); // idx of bin in ax2 
                 
                // everything in 3rd & 4th axes are included - they are the free indexes  
                idx.push_back(ax3Bin);
                idx.push_back(ax4Bin);
                localIdx.push_back(idx);   
            }
        }
        
        // now sum the bin contents in original histogram corresponding to sliced out bins 
        double sumOrig = 0;
        double sumSlice = dim2.Integral();
        for(size_t i = 0; i < localIdx.size(); i++){
            // convert from local to global (flat) idx 
            size_t globalIdx = origHistogram.FlattenIndices(localIdx.at(i)); 

            // get contents in origHistogram at that flattened idx and sum 
            sumOrig += origHistogram.GetBinContent(globalIdx); 
        }

        REQUIRE(dim2.GetNDims() == 2); // check it's a 2D slice
        REQUIRE(newAx1.GetNBins() == origAx3.GetNBins()); // check sliced out axes are same
        REQUIRE(newAx2.GetNBins() == origAx4.GetNBins()); // as in original histogram   
        REQUIRE(sumOrig == sumSlice); // verify bin contents are as expected 
    }

    SECTION("3D Slice"){
        // slice over 1 axes to yield 3D slice 
        std::map<std::string, size_t> fixedBins; 
        fixedBins["ax1"] = 1;
        
        Histogram dim3 = origHistogram.GetSlice(fixedBins);
        BinAxis newAx1 = dim3.GetAxes().GetAxis(0);
        BinAxis newAx2 = dim3.GetAxes().GetAxis(1);  
        BinAxis newAx3 = dim3.GetAxes().GetAxis(2); 

        // CODE TO VERIFY BIN CONTENTS OF SLICE IS CORRECT 
        // vector stores the local idx for each bin in slice for each axis 
        std::vector<std::vector<size_t> > localIdx; 

        // loop over the "free index", ie the axis not present in fixedBins
        for(size_t ax4Bin = 0; ax4Bin < 3; ax4Bin++){
            for(size_t ax3Bin = 0; ax3Bin < 3; ax3Bin++){
                for(size_t ax2Bin = 0; ax2Bin < 3; ax2Bin++){
                    std::vector<size_t> idx;
                    // this is the constant bin idx - where the slice 'pathway' is defined 
                    idx.push_back(1); // idx of bin in ax1 
                    
                    // everything in 4th, 3rd & 2nd axes are included - they are the free indexes 
                    idx.push_back(ax2Bin);  
                    idx.push_back(ax3Bin);
                    idx.push_back(ax4Bin);
                    localIdx.push_back(idx);   
                }
            }
        }
        
        // now sum the bin contents in original histogram corresponding to sliced out bins 
        double sumOrig = 0;
        double sumSlice = dim3.Integral();
        for(size_t i = 0; i < localIdx.size(); i++){
            // convert from local to global (flat) idx 
            size_t globalIdx = origHistogram.FlattenIndices(localIdx.at(i)); 

            // get contents in origHistogram at that flattened idx and sum 
            sumOrig += origHistogram.GetBinContent(globalIdx); 
        }

        REQUIRE(dim3.GetNDims() == 3); // check it's a 3D slice
        REQUIRE(newAx1.GetNBins() == origAx2.GetNBins()); // check sliced out axes are same
        REQUIRE(newAx2.GetNBins() == origAx3.GetNBins()); // as in original histogram   
        REQUIRE(newAx3.GetNBins() == origAx4.GetNBins());
        REQUIRE(sumOrig == sumSlice); // verify bin contents of slice 
    }
}


