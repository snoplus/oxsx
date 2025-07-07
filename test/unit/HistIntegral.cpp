#include <catch2/catch_all.hpp>
#include <HistTools.h>
#include <Combinations.hpp>
#include <Formatter.hpp>
#include <BinAxis.h>
#include <Histogram.h>
#include <IO.h>
#include <iostream>
#include <Rand.h>

//TEST_CASE --- SECTION --> TEXT
//REQUIRE --> boolean -- arg->conditions

TEST_CASE("Bounded Integral Test")
{
    // make the axes
    AxisCollection Axes;
    int NBINS = 100;
    BinAxis var0("Variable0",0.0,1.0,NBINS);
    Axes.AddAxis(var0);
    Histogram hist(Axes);

    //IF WE WANT TO VERIFY FOR HIGHER DIMENSIONS
    //BinAxis var1("Variable1",2.0,3.5,100);
    //Axes.AddAxis(var1);
    //Histogram hist2D(Axes);
    
    //filling the histogram with a known function
    for(int i=0; i<= hist.GetNBins(); i++){
      hist.SetBinContent(i, 1.0/(100*100)*i*i ); //function --> x^2
    }

    IO::SaveHistogramROOT(hist,"hTestUnit.root");
    SECTION("1D")
    {
      REQUIRE(hist.Integral(0.2,0.5)/NBINS < 0.045);//does it calculate the result with moderate accuracy? 
      REQUIRE(hist.Integral(0.2,0.5)/NBINS > 0.035);//(Integral(0.2,0.5) of x^2 should be 0.039)
      REQUIRE(hist.Integral(0.0,1.0) == hist.Integral());//consistent with other functions?
      REQUIRE(hist.Integral(0.1,2.0) == hist.Integral(0.1,1.0));//protected to overflow arguments?
      REQUIRE(hist.Integral(-1.0,0.7) == hist.Integral(0.0,0.7));//protected to underflow/negative arguments?
      REQUIRE(hist.Integral(1.0,0.0) == 0);//protected to max < min arguments?
      REQUIRE(hist.Integral(0.7,0.7) == 0);//consistent with math? 
    }
    /*SECTION("2D")
      {
	//hist2D.Integral(2.2,3.3); // does it throw error??
	}*/
}
