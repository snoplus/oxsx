#include <catch.hpp>
#include <Shape.h>
#include <ParameterDict.h>

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
}