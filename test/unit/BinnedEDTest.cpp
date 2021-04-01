#include <catch.hpp>
#include <BinnedED.h>
#include <Event.h>

TEST_CASE("Filling a 2x2 PDF"){
    AxisCollection ax;
    ax.AddAxis(BinAxis("axis1", 0, 10 , 100));
    ax.AddAxis(BinAxis("axis2", -12, 34 , 100));

    BinnedED pdf("test", ax);

    SECTION("Intial Binning Correct"){
        REQUIRE(pdf.GetNBins() == 100 * 100);
    }

    SECTION("Filling with weights and a vector of indicies"){
        for(size_t i = 0; i < 100; i++){
            std::vector<double> vals;
            vals.push_back(i);
            vals.push_back(i+1);
            
            pdf.Fill(vals, 0.36);
        }

        REQUIRE(pdf.Integral() == Approx(100 * 0.36));
        

        SECTION("Then Normalise"){
            pdf.Normalise();
            REQUIRE(pdf.Integral() == 1);
        }

        SECTION("Then Empty"){
            pdf.Empty();
            REQUIRE(pdf.Integral() == 0);
        }
    }


    SECTION("Filling from Event with weights"){
        std::vector<std::string> relevantIndicies;
        relevantIndicies.push_back("obs0");
        relevantIndicies.push_back("obs3");
        
        pdf.SetObservables(relevantIndicies);

        for(size_t i = 0; i < 100; i++){
            std::vector<double> vals;
            vals.push_back(i);
            vals.push_back(-1);
            vals.push_back(-1);
            vals.push_back(i + 1);
            Event evData(vals);

            std::vector<std::string> observablesEvent;
            observablesEvent.push_back("obs0");
            observablesEvent.push_back("obs1");
            observablesEvent.push_back("obs2");
            observablesEvent.push_back("obs3");
            evData.SetObservableNames(&observablesEvent);

            pdf.Fill(evData, 0.36);
        }
        
        REQUIRE(pdf.Integral() == Approx(0.36 * 100));

    }


    SECTION("Bin Content initialisation and setting"){
        double initialContent = pdf.GetBinContent(1);

        REQUIRE(initialContent == 0);
        
        pdf.AddBinContent(1, 10.01);
        double newContent = pdf.GetBinContent(1);
        REQUIRE(newContent == 10.01);

        double setContent = 20.9;
        pdf.SetBinContent(1, 20.9);
        newContent = pdf.GetBinContent(1);

        REQUIRE(newContent == 20.9);
        REQUIRE(pdf.Integral() == Approx(20.9));
    }

}


