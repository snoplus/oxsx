#include <catch2/catch_approx.hpp>
#include <catch2/catch_all.hpp>
#include <HistTools.h>
#include <BinAxis.h>
#include <Histogram.h>
#include <Combinations.hpp>
#include <Formatter.hpp>

TEST_CASE("Producing  histograms from a set of axes"){    
    // make the axes
    AxisCollection axes;
    std::stringstream ss;
    for(int i = 0; i < 4; i++)
        axes.AddAxis(BinAxis(Formatter() << i, 0, 1, 1));
    
    std::vector<std::string> dimensionList = axes.GetAxisNames();

    SECTION("2D"){
        typedef std::set<std::pair<std::string, std::string> > CombSet;
        CombSet l2combinations =  Combinations::AllCombsNoDiag(dimensionList);
        
        std::map<std::string, Histogram> hists = HistTools::MakeAllHists(axes, l2combinations);
        int i = 0;
        for(CombSet::iterator it = l2combinations.begin(); it != l2combinations.end(); ++it){
            REQUIRE(hists[it->first + "_" + it->second].GetAxes().GetAxis(0).GetName() == it->first);
            REQUIRE(hists[it->first + "_" + it->second].GetAxes().GetAxis(1).GetName() == it->second);
            i++;
        }
    
    }
    
    SECTION("1D"){
        std::map<std::string, Histogram> hists = HistTools::MakeAllHists(axes, dimensionList);
        for(size_t i = 0; i < hists.size(); i++){
            REQUIRE(hists[dimensionList.at(i)].GetAxes().GetAxis(0).GetName() == dimensionList.at(i));
        }
    }
        
}
    

