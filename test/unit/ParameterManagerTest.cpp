#include <catch2/catch_all.hpp>
#include <list>
#include <ParameterManager.h>
#include <DoubleParameter.h>
#include <vector>
#include <set>
#include <iostream>

    
TEST_CASE("Do parameters register correctly?"){
    ParameterManager paramMan;
    std::vector<double> vecParams(10, 0);
    double doubleParam = 0;
    double doubleParam2 = 0;
    std::list<double> listParams(10, 0);
    
    paramMan.Add(new DoubleParameter(doubleParam), "double");
    paramMan.AddDouble(doubleParam2, "double2");

    paramMan.AddContainer(listParams, "list");
    for(size_t i = 0; i < vecParams.size(); i++){
        paramMan.Add(
                     new ContainerParameter<std::vector<double> >(vecParams,i),
                     Formatter() << "container_" << i
                     );
    }

    SECTION("Parameter values"){
        ParameterDict testP;
        testP["double"] = 0;
        testP["double2"] = 0;        
        for(size_t i = 0; i < listParams.size(); i++)
            testP[Formatter() << "list_" << i] = 0;

        for(size_t i = 0; i < vecParams.size(); i++)
            testP[Formatter() << "container_" << i] = 0;

        paramMan.SetParameters(testP);
            
        REQUIRE(paramMan.GetParameterCount() == 22);
        REQUIRE(paramMan.GetParameters() == testP);
    }
    SECTION("Parameter Names"){
        std::set<std::string> expectedNames;
        expectedNames.insert("double");    
        expectedNames.insert("double2");
        for(size_t i = 0; i < listParams.size(); i++){
            std::stringstream ss;
            ss << "list_" << i;
            expectedNames.insert(ss.str());
        }

        for(size_t i = 0; i < vecParams.size(); i++){
            std::stringstream ss;
            ss << "container_" << i;
            expectedNames.insert(ss.str());
        }

        std::set<std::string> names = paramMan.GetParameterNames();
        REQUIRE(names == expectedNames);

    }

    SECTION("after copy"){
        ParameterManager paramManCopy(paramMan);
        REQUIRE(paramManCopy.GetParameterCount() == 0);
    }



    SECTION("setting parameter values"){        
        ParameterDict testP;
        testP["double"] = 5;
        testP["double2"] = 5;

        for(size_t i = 0; i < listParams.size(); i++)
            testP[Formatter() << "list_" << i] = 20;

        for(size_t i = 0; i < vecParams.size(); i++)
            testP[Formatter() << "container_" << i] = 10;


        paramMan.SetParameters(testP);
        REQUIRE(paramMan.GetParameters() == testP);
        REQUIRE(doubleParam == 5);
        REQUIRE(doubleParam2 == 5);

        REQUIRE(listParams == std::list<double> (10, 20));
        REQUIRE(vecParams == std::vector<double> (10, 10));

        // do it again - checks for state dependence
        paramMan.SetParameters(testP);
        REQUIRE(paramMan.GetParameters() == testP);
        REQUIRE(doubleParam == 5);
        REQUIRE(doubleParam2 == 5);
        REQUIRE(listParams == std::list<double> (10, 20));
        REQUIRE(vecParams == std::vector<double> (10, 10));

    
    
    }

}
