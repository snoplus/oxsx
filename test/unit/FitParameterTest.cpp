#include <catch2/catch_all.hpp>
#include <ContainerParameter.h>
#include <DoubleParameter.h>
#include <list>

TEST_CASE("Double Parameter"){
    double a = 5;
    DoubleParameter dParam(a);
    
    REQUIRE(dParam.Get() == 5);
    dParam.Set(2);
    REQUIRE(dParam.Get() == 2);
}

TEST_CASE("vector parameter"){
    std::vector<double> testV(100, 0);
    ContainerParameter<std::vector<double> > cntrParam(testV, 50);
    SECTION("no reallocation"){
        REQUIRE(cntrParam.Get() == 0);
        cntrParam.Set(10);
        REQUIRE(cntrParam.Get() == 10);
    }

    SECTION("after reallocation"){
        testV.reserve(1000);
        REQUIRE(cntrParam.Get() == 0);
    }
}

TEST_CASE("list parameter"){
    std::list<double> testL(100, 0);
    ContainerParameter<std::list<double> > cntrParam(testL, 50);
    SECTION("no reallocation"){
        REQUIRE(cntrParam.Get() == 0);
        cntrParam.Set(10);
        REQUIRE(cntrParam.Get() == 10);
    }

    SECTION("after reallocation"){
        testL.resize(1000);
        REQUIRE(cntrParam.Get() == 0);
    }
}
