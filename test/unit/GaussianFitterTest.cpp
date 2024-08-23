#include <catch2/catch_all.hpp>
#include <Gaussian.h>
#include <ParameterDict.h>
#include <iostream>

TEST_CASE("GaussianFitter", "[GaussianFitter]"){
    std::vector<double> means;
    means.push_back(0);
    means.push_back(1);

    std::vector<double> stDevs;
    stDevs.push_back(2);
    stDevs.push_back(3);


    Gaussian gaus(means, stDevs);
    gaus.RenameParameter("means_0","mu0");
    gaus.RenameParameter("means_1","mu1");
    gaus.RenameParameter("stddevs_0","sigma0");
    gaus.RenameParameter("stddevs_1","sigma1");
    SECTION("Check parameter renaming"){

        REQUIRE(gaus.GetMeanNames().at(0) == "mu0");
        REQUIRE(gaus.GetStDevNames().at(1) == "sigma1");
    }

    ParameterDict values;
    values["mu0"]= 100;
    values["mu1"]= 200;
    values["sigma0"]= 300;
    values["sigma1"]= 400;

    gaus.SetParameters(values);

    SECTION("Check setting parameters correctly"){

        REQUIRE(gaus.GetMean(0) == 100);
        REQUIRE(gaus.GetStDev(1) == 400);
    }

}

