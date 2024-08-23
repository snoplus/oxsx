#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>
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

TEST_CASE("Another simple scale systematic on 2d PDF"){
    AxisCollection axes;
    axes.AddAxis(BinAxis("axis0", 0, 5, 5));
    axes.AddAxis(BinAxis("axis1", 1, 3, 2));
    std::vector<std::string> observables;
    observables.push_back("obs0");
    observables.push_back("obs1");

    BinnedED pdf1("pdf1", axes);
    pdf1.SetBinContent(axes.FlattenIndices({1,0}), 11);
    pdf1.SetBinContent(axes.FlattenIndices({2,1}), 22);
    pdf1.SetObservables(observables);

    Scale* scale = new Scale("scale");
    scale->SetScaleFactor(1.10);
    scale->SetAxes(axes);
    scale->SetTransformationObs({"obs0"});
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
    std::vector<double> correctVals(10,0);
    correctVals[axes.FlattenIndices({1,0})] = 9;
    correctVals[axes.FlattenIndices({2,0})] = 2;
    correctVals[axes.FlattenIndices({2,1})] = 16;
    correctVals[axes.FlattenIndices({3,1})] = 6;
    
    for (size_t i=0; i<modifiedObs.size(); i++) {
        REQUIRE(modifiedObs.at(i) == Catch::Approx(correctVals.at(i)));
    }
}

TEST_CASE("Another simple scale systematic on 3d PDF"){
    AxisCollection axes;
    axes.AddAxis(BinAxis("axis0", 0, 5, 5));
    axes.AddAxis(BinAxis("axis1", 1, 3, 2));
    axes.AddAxis(BinAxis("axis2", 1, 13, 12));
    std::vector<std::string> observables;
    observables.push_back("obs0");
    observables.push_back("obs1");
    observables.push_back("obs2");

    BinnedED pdf1("pdf1", axes);
    pdf1.SetBinContent(axes.FlattenIndices({1,0,5}), 11);
    pdf1.SetBinContent(axes.FlattenIndices({2,1,3}), 22);
    pdf1.SetObservables(observables);

    Scale* scale = new Scale("scale");
    scale->SetScaleFactor(1.10);
    scale->SetAxes(axes);
    scale->SetTransformationObs({"obs0"});
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
    std::vector<double> correctVals(axes.GetNBins(),0);
    correctVals[axes.FlattenIndices({1,0,5})] = 9;
    correctVals[axes.FlattenIndices({2,0,5})] = 2;
    correctVals[axes.FlattenIndices({2,1,3})] = 16;
    correctVals[axes.FlattenIndices({3,1,3})] = 6;
    
    for (size_t i=0; i<modifiedObs.size(); i++) {
        REQUIRE(modifiedObs.at(i) == Catch::Approx(correctVals.at(i)));
    }
}
