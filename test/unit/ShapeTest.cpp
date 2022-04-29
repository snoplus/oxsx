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
    axes.AddAxis(BinAxis("x", -10, 10 ,50));
    axes.AddAxis(BinAxis("y", 20, 50, 10));

    ObsSet dist_obs(std::vector<std::string>({"x", "y"}));
    ObsSet trans_obs("x");

    const std::vector<std::string> shape_param_names {"a", "b"};

    Shape shape_sys("shape_sys");
    shape_sys.SetAxes(axes);
    shape_sys.SetDistributionObs(dist_obs);
    shape_sys.SetTransformationObs(trans_obs);
    shape_sys.SetShapeFunction(linear_func, shape_param_names);

    shape_sys.SetParameters({{"a", 1}, {"b", 1000}});

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
}