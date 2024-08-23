#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>
#include <Shape.h>
#include <ParameterDict.h>
#include <Gaussian.h>
#include <DistTools.h>
#include <BinnedNLLH.h>
#include <OXSXDataSet.h>

double linear_func(const ParameterDict& params, const std::vector<double>& obs_vals) {
    /*
    * Simple non-trivial function, linear in x:
    * returns = a*x + b
    */
    return params.at("a") * obs_vals.at(0) + params.at("b");
}

TEST_CASE("Shape") {
    AxisCollection axes;
    const size_t n_bins_x = 50;
    const size_t n_bins_y = 10;
    const size_t tot_n_bins = n_bins_x*n_bins_y;
    axes.AddAxis(BinAxis("x", -10, 10 ,50));
    axes.AddAxis(BinAxis("y", 20, 50, 10));

    ObsSet dist_obs(std::vector<std::string>({"x", "y"}));
    ObsSet trans_obs("x");

    const std::vector<std::string> shape_param_names {"a", "b"};
    const std::set<std::string> shape_param_names_s {"a", "b"};

    Shape shape_sys("shape_sys");
    shape_sys.SetAxes(axes);
    shape_sys.SetDistributionObs(dist_obs);
    shape_sys.SetTransformationObs(trans_obs);
    shape_sys.SetShapeFunction(linear_func, shape_param_names);
    const auto params = ParameterDict({{"a", 1}, {"b", 1000}});
    shape_sys.SetParameters(params);

    shape_sys.Construct();

    SECTION("Check systematic name") {
        REQUIRE( shape_sys.GetName() == "shape_sys" );
    }

    SECTION("Check observables used") {
        REQUIRE( shape_sys.GetDistributionObs() == dist_obs );
        REQUIRE( shape_sys.GetTransformationObs() == trans_obs );
    }

    SECTION("Check axes") {
        REQUIRE( shape_sys.GetAxes() == axes );
    }

    SECTION("Check argument information") {
        REQUIRE( shape_sys.GetParameterCount() == 2 );
        REQUIRE( shape_sys.GetParameterNames() == shape_param_names_s );
        REQUIRE( shape_sys.GetParameter("a") == 1 );
        REQUIRE( shape_sys.GetParameter("b") == 1000 );
        REQUIRE( shape_sys.GetParameters() == params );
    }

    SECTION("Check constructed response matrix is exactly correct") {
        const auto response_matrix = shape_sys.GetResponse();
        REQUIRE( response_matrix.GetNCols() == tot_n_bins );
        REQUIRE( response_matrix.GetNRows() == tot_n_bins );
        for (size_t i = 0; i < tot_n_bins; i++) {
            for (size_t j = 0; j < tot_n_bins; j++) {
                double scale;
                if (i == j) {
                    const auto x = axes.GetBinCentre(i, 0);
                    scale = linear_func(params, {x});
                } else { scale = 0; }
                REQUIRE( response_matrix.GetComponent(i, j) == scale );
            }
        }
    }

    SECTION("Check interface with BinnedNLLH") {
        const Gaussian gauss_sig({0, 30}, {3, 20}, "signal");
        const Gaussian gauss_back({5, 40}, {10, 5}, "background");

        BinnedED pdf_sig("signal", DistTools::ToHist(gauss_sig, axes));
        BinnedED pdf_back("back", DistTools::ToHist(gauss_back, axes));
        pdf_sig.Normalise();
        pdf_back.Normalise();

        const std::vector<std::string> observables {"obs0", "obs1"};
        pdf_sig.SetObservables(observables);
        pdf_back.SetObservables(observables);

        BinnedNLLH lh;
        lh.AddSystematic(&shape_sys, "shape_group");
        lh.AddPdf(pdf_sig, {"shape_group"}, FALSE);
        lh.AddPdf(pdf_back);

        const std::vector<double> test_point {0, 30};
        const size_t central_bin = pdf_sig.FindBin(test_point);
        double shape_norm_factor;
        BinnedED pdf_sig_mod = shape_sys(pdf_sig, &shape_norm_factor);
        double prob_sig = pdf_sig_mod.GetBinContent(central_bin);
        double prob_back = pdf_back.GetBinContent(central_bin);

        double sum_log_prob = -log(prob_back + prob_sig);
        double sum_norm = shape_norm_factor + 1.;

        OXSXDataSet data;
        data.AddEntry(Event(test_point));
        data.SetObservableNames(observables);
        lh.SetDataSet(&data);

        lh.RegisterFitComponents();
        
        lh.SetParameters({{"back", 1}, {"a", 1}, {"b", 1000}});

        REQUIRE(lh.Evaluate() == Catch::Approx(sum_log_prob + sum_norm));

        // Now change params, and check if lh can still handle it!
        lh.SetParameters({{"back", 2}, {"a", 1}, {"b", 500}});
        auto shape_sys_2(shape_sys);
        shape_sys_2.SetParameters({{"a", 1}, {"b", 500}});
        shape_sys_2.Construct();
        pdf_sig_mod = shape_sys_2(pdf_sig, &shape_norm_factor);
        prob_sig = pdf_sig_mod.GetBinContent(central_bin);
        prob_back = pdf_back.GetBinContent(central_bin)*2.;
        sum_log_prob = -log(prob_back + prob_sig);
        sum_norm = shape_norm_factor + 2.;

        REQUIRE(lh.Evaluate() == Catch::Approx(sum_log_prob + sum_norm));
    }
}
