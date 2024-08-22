#include <catch2/catch_approx.hpp>
#include <catch2/catch_all.hpp>
#include <SystematicManager.h>
#include <ScaleFunction.h>

double linear_func(const ParameterDict& params, const double& obs_val) {
    /*
    * Simple non-trivial function, linear in x:
    * returns = a*x + b
    */
    return params.at("a") * obs_val + params.at("b");
}

TEST_CASE("Simple ScaleFunction systematic on 1d PDF"){
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

    const std::vector<std::string> scale_param_names {"a", "b"};
    const std::set<std::string> scale_param_names_set {"a", "b"};

    ScaleFunction* scale_func = new ScaleFunction("scale_func_sys");
    scale_func->SetScaleFunction(linear_func, scale_param_names);
    const auto params = ParameterDict({{"a", 2}, {"b", 0}});
    scale_func->SetParameters(params);
    scale_func->SetAxes(axes);
    scale_func->SetTransformationObs(observables);
    scale_func->SetDistributionObs(observables);
    scale_func->Construct();

    SECTION("Check systematic name") {
        REQUIRE( scale_func->GetName() == "scale_func_sys" );
    }

    SECTION("Check observables used") {
        REQUIRE( scale_func->GetDistributionObs() == observables );
        REQUIRE( scale_func->GetTransformationObs() == observables );
    }

    SECTION("Check axes") {
        REQUIRE( scale_func->GetAxes() == axes );
    }

    SECTION("Check argument information") {
        REQUIRE( scale_func->GetParameterCount() == 2 );
        REQUIRE( scale_func->GetParameterNames() == scale_param_names_set );
        REQUIRE( scale_func->GetParameter("a") == 2 );
        REQUIRE( scale_func->GetParameter("b") == 0 );
        REQUIRE( scale_func->GetParameters() == params );
    }

    SystematicManager man;
    Systematic* syst1;
    syst1 = scale_func;
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

TEST_CASE("Another simple ScaleFunction systematic on 2d PDF"){
    AxisCollection axes;
    axes.AddAxis(BinAxis("axis0", 0, 5, 5));
    axes.AddAxis(BinAxis("axis1", 1, 3, 2));
    std::vector<std::string> observables;
    observables.push_back("obs0");
    observables.push_back("obs1");

    BinnedED pdf1("pdf1", axes);
    pdf1.SetBinContent(axes.FlattenIndices({0,0}), 12);
    pdf1.SetBinContent(axes.FlattenIndices({1,1}), 6);
    pdf1.SetObservables(observables);

    const std::vector<std::string> scale_param_names {"a", "b"};
    const std::set<std::string> scale_param_names_set {"a", "b"};

    ScaleFunction* scale_func = new ScaleFunction("scale_func_sys");
    scale_func->SetScaleFunction(linear_func, scale_param_names);
    const auto params = ParameterDict({{"a", 3}, {"b", 1}});
    scale_func->SetParameters(params);
    scale_func->SetAxes(axes);
    scale_func->SetTransformationObs({"obs0"});
    scale_func->SetDistributionObs(observables);
    scale_func->Construct();

    SystematicManager man;
    Systematic* syst1;
    syst1 = scale_func;
    man.Add(syst1);
    man.AddDist(pdf1,"");
    man.Construct();

    std::vector<BinnedED> pdfs;
    pdfs.push_back(pdf1);
    std::vector<BinnedED> OrignalPdfs(pdfs);

    man.DistortEDs(OrignalPdfs,pdfs);

    std::vector<double> modifiedObs = pdfs.at(0).GetBinContents();
    std::vector<double> correctVals(10,0);
    correctVals[axes.FlattenIndices({1,0})] = 4;
    correctVals[axes.FlattenIndices({2,0})] = 4;
    correctVals[axes.FlattenIndices({3,0})] = 4;
    correctVals[axes.FlattenIndices({4,1})] = 2;
    
    for (size_t i=0; i<modifiedObs.size(); i++) {
        REQUIRE(modifiedObs.at(i) == Catch::Approx(correctVals.at(i)));
    }
}
