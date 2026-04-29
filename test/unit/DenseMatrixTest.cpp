#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>
#include <DenseMatrix.h>
#include <Exceptions.h>

TEST_CASE("DenseMatrix tests")
{
    DenseMatrix A(2, 3);

    SECTION("Confirm Initialisation Size of Matrix A")
    {
        REQUIRE(A.GetNRows() == 2);
        REQUIRE(A.GetNCols() == 3);
        REQUIRE(A.GetComponent(0,0) == 0.);
        REQUIRE(A.GetComponent(1,0) == 0.);
    }

    SECTION("Set/Get entries, 1-by-1")
    {
        A.SetComponent(0, 0, 2.5);
        A.SetComponent(1, 0, 1.2);

        REQUIRE(A.GetComponent(0, 0) == 2.5);
        REQUIRE(A.GetComponent(1, 0) == 1.2);
    }

    SECTION("Reset to zeroes")
    {
        A.SetZeros();
        REQUIRE(A.GetComponent(0,0) == 0.);
        REQUIRE(A.GetComponent(1,0) == 0.);
    }

    SECTION("Confirm that SetToIdentity fails on non-square matrix")
    {
        REQUIRE_THROWS_AS(A.SetToIdentity(), DimensionError);
    }

    DenseMatrix B(2, 2);
    
    SECTION("SetToIdentity using square-matrix")
    {
        B.SetToIdentity();
        REQUIRE(B.GetComponent(0, 0) == 1);
        REQUIRE(B.GetComponent(1, 0) == 0);
        REQUIRE(B.GetComponent(0, 1) == 0);
        REQUIRE(B.GetComponent(1, 1) == 1);
    }

    SECTION("SetComponent and apply to vector")
    {
        B.SetComponent(0, 0, 2);
        B.SetComponent(1, 0, 1);
        const std::vector<double> x {1, 2};
        const std::vector<double> y = B(x);
        const std::vector<double> y_exp = {2.*1+0.*2, 1.*1+0*2};
        REQUIRE_THAT(y.at(0), Catch::Matchers::WithinAbs(y_exp.at(0), 0.00001));
        REQUIRE_THAT(y.at(1), Catch::Matchers::WithinAbs(y_exp.at(1), 0.00001));
    }

    SECTION("Multiply two matrices, *=")
    {
        A.SetComponent(0, 0, 1);
        A.SetComponent(0, 1, 1);
        A.SetComponent(0, 2, 1);
        
        B.SetComponent(0, 0, 2);
        B.SetComponent(1, 0, 1);

        B *= A;
        
        REQUIRE(B.GetNRows() == 2);
        REQUIRE(B.GetNCols() == 3);
        REQUIRE_THAT(B.GetComponent(0, 0), Catch::Matchers::WithinAbs(2.*1+0*1, 0.00001));
        REQUIRE_THAT(B.GetComponent(0, 1), Catch::Matchers::WithinAbs(2.*1+0*1, 0.00001));
        REQUIRE_THAT(B.GetComponent(0, 2), Catch::Matchers::WithinAbs(2.*1+0*1, 0.00001));
        REQUIRE_THAT(B.GetComponent(1, 0), Catch::Matchers::WithinAbs(1.*1+0*0, 0.00001));
        REQUIRE_THAT(B.GetComponent(1, 1), Catch::Matchers::WithinAbs(1.*1+0*0, 0.00001));
        REQUIRE_THAT(B.GetComponent(1, 2), Catch::Matchers::WithinAbs(1.*1+0*0, 0.00001));
    }
}