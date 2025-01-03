#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>
#include <QuadraticConstraint.h>
#include <BivariateQuadraticConstraint.h>
#include <ConstraintManager.h>
#include <ParameterDict.h>
#include <iostream>

TEST_CASE("Constraints")
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
}
