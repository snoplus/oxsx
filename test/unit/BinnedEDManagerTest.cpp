#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>
#include <BinnedEDManager.h>
#include <DistTools.h>
#include <Gaussian.h>
#include <iostream>

TEST_CASE("Three pdfs no systematics"){
    Gaussian gaus1(0, 1);
    Gaussian gaus2(2, 3);
    Gaussian gaus3(4, 5);

    AxisCollection axes;
    axes.AddAxis(BinAxis("axis1", -40, 40 , 200));

    BinnedED pdf1("a", DistTools::ToHist(gaus1, axes));
    BinnedED pdf2("b", DistTools::ToHist(gaus2, axes));
    BinnedED pdf3("c", DistTools::ToHist(gaus3, axes));

    std::vector<std::string> observable;
    observable.push_back("obs0");

    pdf1.SetObservables(observable);
    pdf2.SetObservables(observable);    
    pdf3.SetObservables(observable);

    double prob1 = pdf1.GetBinContent(0);
    double prob2 = pdf2.GetBinContent(0);    
    double prob3 = pdf3.GetBinContent(0);

    BinnedEDManager edMan;
    edMan.AddPdf(pdf1);
    edMan.AddPdf(pdf2);
    edMan.AddPdf(pdf3);
    edMan.SetNormalisations(std::vector<double>(3, 1));

    SECTION("Bin Probability Method"){
      REQUIRE(edMan.BinProbability(0) == Catch::Approx(prob1 * prob2 * prob3));
    }
    
    SECTION("Probability Method"){
        Event ev(std::vector<double>(1, -39.5));
        std::vector<std::string> observablesEvent;
        observablesEvent.push_back("obs0");
        observablesEvent.push_back("obs1");
        ev.SetObservableNames(&observablesEvent);
        REQUIRE(edMan.Probability(ev) == Catch::Approx(prob1 * prob2 * prob3));
    }

}
