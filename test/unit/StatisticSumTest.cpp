#include <catch.hpp>
#include <StatisticSum.h>

class FakeStatistic : public TestStatistic{
public:
    double Evaluate() {
        return fVal;
    }

    size_t GetParameterCount() const {
        return 1;
    }

    void SetParameters(const ParameterDict& p){
        fVal = p.at(fParamName);
    }

    ParameterDict GetParameters() const{
        ParameterDict p; 
        p[fParamName] = fVal;
        return p;
    } 

    std::set<std::string> GetParameterNames() const{
        std::set<std::string> set;
        set.insert(fParamName);
        return set;
    }

    void RegisterFitComponents() {}

    double fVal;
    std::string fParamName;
};

TEST_CASE("Adding test statistics using StatisticSum constr"){
    FakeStatistic s1;
    FakeStatistic s2;
    
    s1.fVal = 1;
    s1.fParamName = "p1";
    s2.fVal = 2;
    s2.fParamName = "p2";

    SECTION("no shared parameters"){
        StatisticSum sum(s1, s2);
        REQUIRE(sum.GetParameterCount() == 2);
        std::set<std::string> expectedNames;
        expectedNames.insert("p1");
        expectedNames.insert("p2");

        ParameterDict expectedVals;
        expectedVals["p1"] = 1;
        expectedVals["p2"] = 2;

        REQUIRE(sum.GetParameters() == expectedVals);
        REQUIRE(sum.GetParameterNames() == expectedNames);

        REQUIRE(sum.Evaluate() == 3); // 2 + 1
    
        ParameterDict setVals;
        setVals["p1"] = 3;
        setVals["p2"] = 4;

        sum.SetParameters(setVals);
        REQUIRE(s1.fVal == 3);
        REQUIRE(s2.fVal == 4);
        
        REQUIRE(sum.Evaluate() == 7); // 3 + 4
    }
    
    SECTION("setting with shared parameters"){
        StatisticSum sum(s1, s2);
        s1.fParamName = "p2";
        std::set<std::string> expectedNames;
        expectedNames.insert("p2");
        REQUIRE(sum.GetParameterNames() == expectedNames);
        REQUIRE(sum.GetParameterCount() == 1);


        ParameterDict p;
        p["p2"] = 10;
        sum.SetParameters(p);
        REQUIRE(sum.Evaluate() == 20); // 10 + 10
        
    }

    SECTION("adding test stats with +"){
        StatisticSum sum = s1 + s2;
        REQUIRE(sum.GetParameterCount() == 2);
        std::set<std::string> expectedNames;
        expectedNames.insert("p1");
        expectedNames.insert("p2");

        ParameterDict expectedVals;
        expectedVals["p1"] = 1;
        expectedVals["p2"] = 2;

        REQUIRE(sum.GetParameters() == expectedVals);
        REQUIRE(sum.GetParameterNames() == expectedNames);

        REQUIRE(sum.Evaluate() == 3); // 2 + 1    


        FakeStatistic s3;
        s3.fVal = 3;
        s3.fParamName = "p3";

        StatisticSum sum2 = sum + s3;
        StatisticSum sum3 = s3 + sum;

        REQUIRE(sum2.Evaluate() == 6);
        REQUIRE(sum3.Evaluate() == 6);
    }

    SECTION("Test stat sum operator"){
        FakeStatistic s3;
        s3.fVal = 3;
        s3.fParamName = "p3";

        std::vector<TestStatistic*> stats;
        stats.push_back(&s1);
        stats.push_back(&s2);
        stats.push_back(&s3);

        StatisticSum sum = Sum(stats);
        REQUIRE(sum.Evaluate() == 6);
    }

}
