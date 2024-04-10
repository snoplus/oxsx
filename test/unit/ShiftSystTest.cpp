#include <catch.hpp>
#include <SystematicManager.h>
#include <Shift.h>

TEST_CASE("Simple shift systematic on 1d PDF"){
    AxisCollection axes;
    axes.AddAxis(BinAxis("axis0", 0, 5, 5));
    std::vector<std::string> observables;
    observables.push_back("obs0");

    BinnedED pdf1("pdf1", axes);
    pdf1.SetBinContent(0, 0);
    pdf1.SetBinContent(1, 10);
    pdf1.SetBinContent(2, 0);
    pdf1.SetBinContent(3, 0);
    pdf1.SetBinContent(4, 0);
    pdf1.SetObservables(observables);

    Shift* shift = new Shift("shift");
    shift->SetShift(2.0);
    shift->SetAxes(axes);
    shift->SetTransformationObs(observables);
    shift->SetDistributionObs(observables);
    shift->Construct();

    SystematicManager man;

    SECTION("With shift == integer number of bins"){

      Systematic* syst1;
      syst1 = shift;
      man.Add(syst1);
      man.AddDist(pdf1,"");
      man.Construct();
      
      std::vector<BinnedED> pdfs;
      pdfs.push_back(pdf1);
      std::vector<BinnedED> OrignalPdfs(pdfs);
      
      man.DistortEDs(OrignalPdfs,pdfs);
      
      std::vector<double> modifiedObs = pdfs.at(0).GetBinContents();
      std::vector<double> correctVals;
      correctVals.push_back(0);
      correctVals.push_back(0);
      correctVals.push_back(0);
      correctVals.push_back(10);
      correctVals.push_back(0);
      
      REQUIRE(modifiedObs == correctVals);
    }
    
    SECTION("With shift == non-integer number of bins"){
      
      shift->SetShift(1.5);
      Systematic* syst1;
      syst1 = shift;
      man.Add(syst1);
      man.AddDist(pdf1,"");
      man.Construct();

      std::vector<BinnedED> pdfs;
      pdfs.push_back(pdf1);
      std::vector<BinnedED> OrignalPdfs(pdfs);

      man.DistortEDs(OrignalPdfs,pdfs);

      std::vector<double> modifiedObs = pdfs.at(0).GetBinContents();
      std::vector<double> correctVals;
      correctVals.push_back(0);
      correctVals.push_back(0);
      correctVals.push_back(5);
      correctVals.push_back(5);
      correctVals.push_back(0);

      REQUIRE(modifiedObs == correctVals);
    }
}
