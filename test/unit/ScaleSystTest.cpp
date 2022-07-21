#include <catch.hpp>
#include <SystematicManager.h>
#include <Scale.h>

TEST_CASE("Simple scale systematic on 1d PDF"){
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

    Scale* scale = new Scale("scale");
    scale->SetScaleFactor(2.0);
    scale->SetAxes(axes);
    scale->SetTransformationObs(observables);
    scale->SetDistributionObs(observables);
    scale->Construct();

    SystematicManager man;
    Systematic* syst1;
    syst1 = scale;
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
