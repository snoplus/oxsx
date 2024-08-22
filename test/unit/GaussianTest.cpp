#include <catch2/catch_approx.hpp>
#include <catch2/catch_all.hpp>
#include <Gaussian.h>
#include <iostream>

TEST_CASE("2D gaussian", "[Gaussian]"){
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

TEST_CASE("1D gaussian", "[Gaussian]"){
    Gaussian gaus(0, 1);

    SECTION("Check paramter storage"){

        REQUIRE(gaus.GetMean(0) == 0);
        REQUIRE(gaus.GetStDev(0) == 1);
           
    }

    SECTION("Check probability"){
        REQUIRE(gaus(std::vector<double>(1,1)) == Catch::Approx(0.24197));
        REQUIRE(gaus(std::vector<double> (1,1E8)) == Catch::Approx(0));
        REQUIRE(gaus(std::vector<double> (1,-1E8)) == Catch::Approx(0));
    }


    SECTION("Test CDF"){
        REQUIRE(gaus.Cdf(0, 1) == Catch::Approx(0.841344));
        REQUIRE(gaus.Cdf(0, 100) == Catch::Approx(1));
        REQUIRE(gaus.Cdf(0, -100) == Catch::Approx(0));
    }

    SECTION("Test Integral"){
        REQUIRE(gaus.Integral(std::vector<double>(1,-1),std::vector<double>(1,1)) == Catch::Approx(0.6827));
        REQUIRE(gaus.Integral(std::vector<double>(1, -100), std::vector<double>(1, 100)) == Catch::Approx(1));
        REQUIRE(gaus.Integral(std::vector<double>(1,0), std::vector<double>(1,0)) == Catch::Approx(0));
    }
}

