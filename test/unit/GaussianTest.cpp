#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>
#include <Gaussian.h>
#include <iostream>

TEST_CASE("2D gaussian", "[Gaussian]")
{
    std::vector<double> means;
    means.push_back(0);
    means.push_back(1);

    std::vector<double> stDevs;
    stDevs.push_back(2);
    stDevs.push_back(3);

    Gaussian gaus(means, stDevs);

    // set params right
    REQUIRE(gaus.GetMean(0) == 0);
    REQUIRE(gaus.GetMean(1) == 1);

    REQUIRE(gaus.GetStDev(0) == 2);
    REQUIRE(gaus.GetStDev(1) == 3);

    // organised right
    REQUIRE(gaus.GetMeans().size() == 2);
    REQUIRE(gaus.GetStdDevs().size() == 2);
}

TEST_CASE("1D gaussian", "[Gaussian]")
{
    Gaussian gaus(0, 1);

    SECTION("Check paramter storage")
    {

        REQUIRE(gaus.GetMean(0) == 0);
        REQUIRE(gaus.GetStDev(0) == 1);
    }

    SECTION("Check probability")
    {
        REQUIRE_THAT(gaus(std::vector<double>(1, 1)), Catch::Matchers::WithinAbs(0.24197, 0.0001));
        REQUIRE_THAT(gaus(std::vector<double>(1, 1E8)), Catch::Matchers::WithinAbs(0, 0.0001));
        REQUIRE_THAT(gaus(std::vector<double>(1, -1E8)), Catch::Matchers::WithinAbs(0, 0.0001));
    }

    SECTION("Test CDF")
    {
        REQUIRE_THAT(gaus.Cdf(0, 1), Catch::Matchers::WithinAbs(0.841344, 0.0001));
        REQUIRE_THAT(gaus.Cdf(0, 100), Catch::Matchers::WithinAbs(1, 0.0001));
        REQUIRE_THAT(gaus.Cdf(0, -100), Catch::Matchers::WithinAbs(0, 0.0001));
    }

    SECTION("Test Integral")
    {
        REQUIRE_THAT(gaus.Integral(std::vector<double>(1, -1), std::vector<double>(1, 1)), Catch::Matchers::WithinAbs(0.6827, 0.0001));
        REQUIRE_THAT(gaus.Integral(std::vector<double>(1, -100), std::vector<double>(1, 100)), Catch::Matchers::WithinAbs(1, 0.0001));
        REQUIRE_THAT(gaus.Integral(std::vector<double>(1, 0), std::vector<double>(1, 0)), Catch::Matchers::WithinAbs(0, 0.0001));
    }
}
