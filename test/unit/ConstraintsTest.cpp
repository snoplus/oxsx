#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>
#include <QuadraticConstraint.h>
#include <BivariateQuadraticConstraint.h>
#include <ShapeInterpConstraint.h>
#include <ConstraintManager.h>
#include <ParameterDict.h>
#include <iostream>

TEST_CASE("ConstraintManager")
{
    // First: create a ConstraintManager, add 1 individual symmetric constraint.
    ConstraintManager c_man;
    c_man.SetConstraint("a", 1, 2);
    ParameterDict params = {{"a", 2}};
    SECTION("1 Symmetric Constraint Evaluation")
    {
        REQUIRE(c_man.Evaluate(params) == Catch::Approx(1. / 8.));
    }
    // Modify constraint to be asymmetric
    c_man.SetConstraint("a", 1, 1, 2);
    SECTION("1 Asymmetric Constraint Evaluation")
    {
        REQUIRE(c_man.Evaluate(params) == Catch::Approx(1. / 2.));
    }
    // Add a second individual constraint
    c_man.SetConstraint("b", 10, 5);
    params["b"] = 0;
    SECTION("2 Symmetric Indiviudal Constraints Evaluation")
    {
        REQUIRE(c_man.Evaluate(params) == Catch::Approx(0.5 + 2.));
    }
    // Now add a pair constraint, with no correlation
    c_man.SetConstraint("c", 0, 2, "d", 2, 10, 0);
    params["c"] = 10;
    params["d"] = 1;
    SECTION("Uncorrelated pair constraint Evaluation")
    {
        REQUIRE(c_man.Evaluate(params) == Catch::Approx(0.5 + 2. + 25. / 2. + 1. / 200.));
    }
    // Modify pair constraint to handle correlation
    c_man.SetConstraint("c", 0, 2, "d", 2, 10, 0.8);
    SECTION("Uncorrelated pair constraint Evaluation")
    {
        REQUIRE(c_man.Evaluate(params) == Catch::Approx(0.5 + 2. + (25. / 2. + 1. / 200. - 0.8 * 5 * (-1. / 10.)) / (1. - 0.8 * 0.8)));
    }
    // Now add a ratio constraint
    c_man.SetConstraint("c", "d", 5, 1);
    SECTION("Ratio constraint Evaluation")
    {
        REQUIRE(c_man.Evaluate(params) == Catch::Approx(0.5 + 2. + (25. / 2. + 1. / 200. - 0.8 * 5 * (-1. / 10.)) / (1. - 0.8 * 0.8) + (25. / 2.)));
    }
    // Now add a shape constraint
    const BinAxis ax("x", -0.5, 9.5, 10);
    AxisCollection axs;
    axs.AddAxis(ax);
    Histogram hist(axs);
    for (size_t idx=0; idx<10; idx++) { hist.Fill(idx, idx*idx); }
    c_man.SetConstraint("x", hist);
    params["x"] = 4.99;
    SECTION("Shape constraint Evaluation")
    {
        REQUIRE(c_man.Evaluate(params) == Catch::Approx(0.5 + 2. + (25. / 2. + 1. / 200. - 0.8 * 5 * (-1. / 10.)) / (1. - 0.8 * 0.8) + (25. / 2.) + 16.0 + (4.99-4.0)*(25.-16.)/(5.-4.)));
    }
}

TEST_CASE("ShapeInterpConstraints")
{
    SECTION("Simple 1D interpolated constraint")
    {
        // Set up data used for interpolation; let's use a simple y=x^2 distribution
        std::vector<double> ys;
        for (size_t idx=0; idx<10; idx++) { ys.push_back(idx*idx); }
        const BinAxis ax("x", -0.5, 9.5, 10);
        AxisCollection axs;
        axs.AddAxis(ax);
        // Create histogram used for constraint
        Histogram hist(axs);
        hist.SetBinContents(ys);
        // Create constraint object
        const ShapeInterpConstraint con(hist);
        // Try at grid points - should get exact value back
        ParameterDict params;
        params["x"] = 0.;
        REQUIRE(con.Evaluate(params) == Catch::Approx(0.));
        params["x"] = 1.;
        REQUIRE(con.Evaluate(params) == Catch::Approx(1.));
        params["x"] = 2.;
        REQUIRE(con.Evaluate(params) == Catch::Approx(4.));
        // Try some actual interpolation!
        params["x"] = 0.5;
        REQUIRE(con.Evaluate(params) == Catch::Approx(0.5));
        params["x"] = 1.2;
        REQUIRE(con.Evaluate(params) == Catch::Approx(1.0 + (1.2-1.0)*(4.-1.)/(2.-1.)));
        params["x"] = 4.99;
        REQUIRE(con.Evaluate(params) == Catch::Approx(16.0 + (4.99-4.0)*(25.-16.)/(5.-4.)));
    }

    SECTION("2D interpolated constraints")
    {
        // Set up data used for interpolation; let's use a z=x^2 + y^3 + 2 distribution
        const BinAxis ax_x("x", -0.5, 9.5, 10);
        const BinAxis ax_y("y", -1.0, 9.0, 5);
        AxisCollection axs;
        axs.AddAxis(ax_x);
        axs.AddAxis(ax_y);
        Histogram hist(axs);
        for (double x=0; x<10; x++)
        {
            for (double y=0; y<10; y+=2) // 0, 2, 4, 6, 8
            {
                hist.Fill(std::vector<double>({x,y}), x*x + y*y*y + 2.0);
            }
        }
        //Create constraint object
        const ShapeInterpConstraint con(hist);
        //
        ParameterDict params = {{"x",0} , {"y",0}};
        REQUIRE(con.Evaluate(params) == Catch::Approx(2.));
        params = {{"x",1} , {"y",1}};
        REQUIRE(con.Evaluate(params) == Catch::Approx(7.));
        params = {{"x",0.5} , {"y",0.5}};
        REQUIRE(con.Evaluate(params) == Catch::Approx(4.5));
    }
}

// Set up analytical function (a-b) / (a+b)
double frac_func(const ParameterDict &params)
{
    double val = (params.at("a") - params.at("b")) / (params.at("a")+params.at("b"));
    return val;
}

// Set up analytical function, here mimicking a specific use case for the antinu 
// analysis. We want a constraint on the fractional difference between the total
// U and Th geonu rates. There are two fixed parameters for each of U and Th
// which store the nominal rates in the PPO and bisMSB phases. A third parameter
// for each then the scales the two parameters together. So the total rate for
// each is scaling*(nom_bismsb + nom_ppo), and we want a constraint on the
// fractional difference between those. 
double scale_frac_func(const ParameterDict &params)
{
    double total_rate1 = params.at("c") * (params.at("d") + params.at("e"));
    double total_rate2 = params.at("f") * (params.at("g") + params.at("h"));

    double val = ( total_rate1 - total_rate2 ) / ( total_rate1 + total_rate2 );
    return val;
}

TEST_CASE("ShapeConstraints")
{
    SECTION("Simple 1D analytical constraint")
    {
        ConstraintManager c_man;
        ParameterDict params = {{"a", 3}, {"b", 1}};
        c_man.SetConstraint(params, frac_func, 1, 0.5);
        REQUIRE(c_man.Evaluate(params) == Catch::Approx(0.5));
        params["c"] = 2;
        params["d"] = 3;
        params["e"] = 1;
        params["f"] = 3;
        params["g"] = 1;
        params["h"] = 1./3.;
        c_man.SetConstraint(params, scale_frac_func, 4./3., 0.5);
        REQUIRE(c_man.Evaluate(params) == Catch::Approx( 0.5 + 2 ));
    }

}
