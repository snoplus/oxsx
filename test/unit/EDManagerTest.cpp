#include <catch.hpp>
#include <EDManager.h>
#include <Gaussian.h>
#include <AnalyticED.h>
#include <Event.h>
#include <ContainerTools.hpp>

TEST_CASE("Add a couple of analytic pdfs"){
    EDManager pdfMan;
    Gaussian gaus1(0, 1);
    Gaussian gaus2(0, 1);

    AnalyticED pdf1("g1", &gaus1);
    AnalyticED pdf2("g2", &gaus2);
    std::vector<std::string> observable;
    observable.push_back("obs0");
    pdf1.SetObservables(observable);
    pdf2.SetObservables(observable);
    
    SECTION("initialised correctly"){
        REQUIRE(pdfMan.GetNDists() == 0);
        REQUIRE(pdfMan.GetNDims() == 0);
        
    }
    
    SECTION("add pdfs one at a time"){
        pdfMan.AddDist(&pdf1);
        pdfMan.AddDist(&pdf2);
        
        REQUIRE(pdfMan.GetNDims() == 1);
        REQUIRE(pdfMan.GetNDists() == 2);
        REQUIRE(pdfMan.GetNormalisations() == std::vector<double>(2, 0));
    }
    
    SECTION("add dists in one go"){
        std::vector<EventDistribution*> dists;
        dists.push_back(&pdf1);
        dists.push_back(&pdf2);
        pdfMan.AddDists(dists);
        
        REQUIRE(pdfMan.GetNDims() == 1);
        REQUIRE(pdfMan.GetNDists() == 2);
        REQUIRE(pdfMan.GetNormalisations() == std::vector<double>(2, 0));
    }

    SECTION("correct probability"){
        pdfMan.AddDist(&pdf1);
        pdfMan.AddDist(&pdf2);

        Event event(std::vector<double>(1, 0));
        std::vector<std::string> observablesEvent;
        observablesEvent.push_back("obs0");
        observablesEvent.push_back("obs1");
        event.SetObservableNames(&observablesEvent);

        REQUIRE(pdfMan.Probability(event) == 0); // norms are zero

        pdfMan.SetNormalisations(std::vector<double>(2, 1));
        REQUIRE(pdfMan.Probability(event) == Approx(0.7978845607));
        // 0.7978845607 = 2 /sqrt(2 * pi)
    }

    SECTION("works as fit component"){
        pdfMan.AddDist(&pdf1);
        pdfMan.AddDist(&pdf2);

        REQUIRE(pdfMan.GetParameterCount() == 2);
        ParameterDict testPs;
        testPs["g1_norm"] = 0;
        testPs["g2_norm"] = 0;
        REQUIRE(pdfMan.GetParameters() == testPs);
        
        std::set<std::string> expectedNames;
        expectedNames.insert("g1_norm");
        expectedNames.insert("g2_norm");
        REQUIRE(pdfMan.GetParameterNames() == expectedNames);
        
        testPs["g1_norm"] = 10;
        testPs["g2_norm"] = 15;
        pdfMan.SetParameters(testPs);
        REQUIRE(pdfMan.GetParameters()     == testPs);
        // note the line below only works because the normalisations
        // happen to be in alphabetical order.. don't assume this generally
        REQUIRE(pdfMan.GetNormalisations() == ContainerTools::GetValues(testPs));
    }
}
