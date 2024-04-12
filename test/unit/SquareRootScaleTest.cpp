#include <catch.hpp>
#include <SquareRootScale.h>
#include <iostream>

TEST_CASE("Square Root Scale", "[SquareRootScale]"){
  double gradient = 3;

  SquareRootScale sqrtscalefunc("sqrt_scale");
  sqrtscalefunc.SetGradient(gradient);

  SECTION("Check paramter storage"){

    REQUIRE(sqrtscalefunc.GetParameter("grad") == 3);

  }

  SECTION("Check Function Values"){

    std::vector<double> test_vals;
    test_vals.push_back(4);

    // Sqrt(4) = 2, multiplied by the gradient, 3, gives 6
    REQUIRE ( sqrtscalefunc(test_vals) == 6);

  }
}
