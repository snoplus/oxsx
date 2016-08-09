#include <catch.hpp>
#include <IO.h>
#include <Histogram.h>
#include <Combinations.hpp>

TEST_CASE("Writing a histogram to disk  and reading back"){
    // make any old histogram
    AxisCollection axes;
    axes.AddAxis(PdfAxis("a", 10, 11, 10));
    axes.AddAxis(PdfAxis("b", 20, 21, 10));    
    axes.AddAxis(PdfAxis("c", 30, 31, 10));
    
    Histogram origHisto(axes);
    origHisto.SetBinContents(SequentialElements<double>(0., axes.GetNBins()));
    
    // save it to disk and read it back
    IO::SaveHistogram(origHisto, "tmp_test_histogram_io.h5");
    Histogram loadedHisto = IO::LoadHistogram("testOrigHistogramWrite.h5");

    SECTION("SAME AXES"){
        for(unsigned i = 0; i < origHisto.GetNDims(); i++){
            const PdfAxis& origAxis   =  origHisto.GetAxes().GetAxis(i);
            const PdfAxis& loadedAxis =  loadedHisto.GetAxes().GetAxis(i);
            REQUIRE(origAxis.GetName() == loadedAxis.GetName());
            REQUIRE(origAxis.GetLatexName() == loadedAxis.GetLatexName());
            REQUIRE(origAxis.GetBinLowEdges() == loadedAxis.GetBinLowEdges());
            REQUIRE(origAxis.GetBinHighEdges() == loadedAxis.GetBinHighEdges());
            REQUIRE(origAxis.GetNBins() == loadedAxis.GetNBins());
        }
    }

    SECTION("SAME BIN CONTENTS"){
        REQUIRE(origHisto.GetBinContents() == loadedHisto.GetBinContents());
    }

    remove("tmp_test_histogram_io.h5");
}

