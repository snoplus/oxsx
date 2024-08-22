#include <catch2/catch_approx.hpp>
#include <catch2/catch_all.hpp>
#include <Gaussian.h>
#include <BinnedED.h>
#include <CompositeED.h>
#include <AnalyticED.h>
#include <iostream>

TEST_CASE("Combining 1D gaussians", "[CompositeED]"){
    Gaussian gausF1(0.5, 0.4);
    Gaussian gausF2(0.5, 0.3);

    AnalyticED gaus1("g1", &gausF1);
    AnalyticED gaus2("g2", &gausF2);

    ObsSet d1("obs0");
    ObsSet d2("obs2");
    
    gaus1.SetObservables(d1);
    gaus2.SetObservables(d2);

    CompositeED compositeED = gaus1 * gaus2;

    // test values
    std::vector<double> vals;
    vals.push_back(1);
    vals.push_back(999);
    vals.push_back(1);
    vals.push_back(1);
    vals.push_back(999);

    std::vector<std::string> observables;
    observables.push_back("obs0");
    observables.push_back("obs1");
    observables.push_back("obs2");
    observables.push_back("obs3");
    observables.push_back("obs4");

    Event ev(vals);
    ev.SetObservableNames(&observables);

    SECTION("Check Dimensions"){
        REQUIRE(compositeED.GetNDims() == 2);
    }

    SECTION("Check Data Flow to internal pdfs "){
        double prob  = compositeED.Probability(ev);
        REQUIRE(prob == Catch::Approx(0.15141173681343614));
        
    }

    SECTION("Check Clone Functionality"){
        EventDistribution* clone = compositeED.Clone();
        REQUIRE(clone -> GetNDims() == 2);
        REQUIRE(clone -> Probability(ev) == Catch::Approx(0.15141173681343614));
    }

    SECTION("Second level of recursion"){
        Gaussian nextF = Gaussian(0.9, 0.8);

        AnalyticED nextED("g3", &nextF);
        nextED.SetObservables(ObsSet("obs3"));
        CompositeED level2 = compositeED * nextED;
        REQUIRE( level2.GetNDims() == 3 );
        REQUIRE( level2.Probability(ev) == Catch::Approx(0.07491808959564718));
                                                                                         
    }
}

TEST_CASE(" Composite of two Binned EDs"){

    // Paired into two couplets
    BinAxis axis1("axis1", -80, 80, 100);
    BinAxis axis2("axis2", -80, 80, 100);

    BinAxis axis3("axis3", -80, 80, 100);
    BinAxis axis4("axis4", -80, 80, 100);

    AxisCollection axes1;
    AxisCollection axes2;
    
    axes1.AddAxis(axis1);
    axes1.AddAxis(axis2);
    
    axes2.AddAxis(axis3);
    axes2.AddAxis(axis4);

    BinnedED pdf1("b1", axes1);
    BinnedED pdf2("b2", axes2);


    // Data, where to look
    std::vector<std::string> indicies1;
    std::vector<std::string> indicies2;

    indicies1.push_back("obs0");
    indicies1.push_back("obs2");
    
    indicies2.push_back("obs1");
    indicies2.push_back("obs3");

    pdf1.SetObservables(indicies1);
    pdf2.SetObservables(indicies2);

    // Now combine
    CompositeED compositeED = pdf1 * pdf2;


}

