#include <catch.hpp>
#include <DistTools.h>
#include <BinnedED.h>
#include <Gaussian.h>
#include <TRandom3.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <math.h>
#include <iostream>
#include <TH1D.h>

TEST_CASE("Writing a 1D pdf to a root histogram", "[DistTools]"){
    BinAxis axis("test", -100, 100, 200);
    AxisCollection axes;
    axes.AddAxis(axis);

    BinnedED binnedPdf("test", axes);

    SECTION("Step Function pdf"){
        // fill a heaviside
        for(size_t i = 0; i < 200; i++){
            if(binnedPdf.GetAxes().GetBinLowEdge(i, 0) >= 0)
                binnedPdf.SetBinContent(i,1);
        }

        TH1D rootPdf = DistTools::ToTH1D(binnedPdf);
        REQUIRE(rootPdf.Integral() == 100);
        
        // root histograms are ordered 1->nbins with 0 and nbins +1 as over/underflow
        // nativex pdfs are ordered 0->nbins-1 with 0 and nbins-1 as underflow

        std::vector<double> firstOneHundred;
        for(size_t i = 1; i < 101; i++)
            firstOneHundred.push_back(rootPdf.GetBinContent(i));

        std::vector<double> secondOneHundred;
        for(size_t i = 101; i < 201; i++)
            secondOneHundred.push_back(rootPdf.GetBinContent(i));

        double underflow = rootPdf.GetBinContent(0);
        double overflow  = rootPdf.GetBinContent(201);


        REQUIRE(firstOneHundred == std::vector<double>(100, 0));
        REQUIRE(secondOneHundred == std::vector<double>(100, 1));
        REQUIRE(underflow == 0);
        REQUIRE(overflow == 0);
    }

}
 

TEST_CASE("Converting a 1D gaussian to a binned pdf", "[DistTools]"){    
    BinAxis axis("test", -100, 100, 2000);                                   
    AxisCollection axes;                                                     
    axes.AddAxis(axis);                                                      
                                                                             
    Gaussian gaus(10, 21.1);                                                 
    gaus.SetCdfCutOff(1E8);                                                  
    BinnedED binnedGaus("bg", DistTools::ToHist(gaus, axes));
    double intBinError  = std::abs(binnedGaus.Integral() - 1);
    REQUIRE(intBinError == Approx(0));
    
    SECTION("Numerical Checks"){             
        binnedGaus.Normalise();
        double meanBinError = std::abs(binnedGaus.Means().at(0) - 10);
        double varBinError  = std::abs(binnedGaus.Variances().at(0) - 21.1 * 21.1);

                
        // integral mean and var right to within 0.3% rms != gaus sigma
        REQUIRE(meanBinError/10 < 0.003);
        REQUIRE(varBinError/21.1/21.1 < 0.003);
        REQUIRE(gaus.Integral() == 1);
    }                                                                        
    
}                                                                           

TEST_CASE("Converting 2D gaussian to binned and marginalise", "[DistTools]"){
    BinAxis axis1("ax1", -1000, 1000, 800);
    BinAxis axis2("ax2", -1000, 1000, 800);

    AxisCollection axes;
    axes.AddAxis(axis1);
    axes.AddAxis(axis2);
    
    std::vector<double> means;
    means.push_back(0);
    means.push_back(10);

    std::vector<double> stDevs;
    stDevs.push_back(20);
    stDevs.push_back(30);

    Gaussian gaus(means, stDevs);
    gaus.SetCdfCutOff(1E8); // max accuracy
    BinnedED binnedGaus("bg", DistTools::ToHist(gaus, axes));

    SECTION("Binning Errors"){
        double integralError = std::abs(binnedGaus.Integral() - 1);        
        double mean1Error = std::abs(binnedGaus.Means().at(0));
        double mean2Error = std::abs(binnedGaus.Means().at(1) - 10);
        double stDev1Error = std::abs(binnedGaus.Variances().at(0) - 20 * 20);
        double stDev2Error = std::abs(binnedGaus.Variances().at(1) - 30 * 30);

        // right dims
        REQUIRE(binnedGaus.GetNDims() == 2);

        // 0.5% error
        REQUIRE(integralError == Approx(0));
        REQUIRE(mean1Error == Approx(0));
        REQUIRE(mean2Error == Approx(0));

        // binning tends to give an error on the varaince 0.1% becuase rms != sigma
        REQUIRE(stDev1Error/20/20 < 0.004);
        REQUIRE(stDev2Error/30/30 < 0.004);        
    }


    SECTION("Marginalise"){
        binnedGaus.Normalise();
        std::vector<std::string> obsNames;
        obsNames.push_back("obs0");
        obsNames.push_back("obs1");
        obsNames.push_back("obs2");
        binnedGaus.SetObservables(obsNames);
	
        BinnedED xProj = binnedGaus.Marginalise("ax1");
        BinnedED yProj = binnedGaus.Marginalise("ax2");

        
        double xMean = xProj.Means().at(0);
        double yMean = yProj.Means().at(0);

        double xVarEr = std::abs(xProj.Variances().at(0) - 20 * 20);
        double yVarEr = std::abs(yProj.Variances().at(0) - 30 * 30);

        REQUIRE(xProj.Integral() == Approx(1));
        REQUIRE(yProj.Integral() == Approx(1));

        REQUIRE(xMean == Approx(0));
        REQUIRE(yMean == Approx(10));

        REQUIRE(xVarEr/20/20 < 0.004);
        REQUIRE(yVarEr/30/30 < 0.004);
        

        //        (DistTools::ToTH1D(xProj)).SaveAs("x_proj.root");
        //        (DistTools::ToTH1D(yProj)).SaveAs("y_proj.root");
    }
}
