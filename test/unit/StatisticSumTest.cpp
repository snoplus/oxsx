#include <catch2/catch_approx.hpp>
#include <catch2/catch_all.hpp>
#include <StatisticSum.h>
#include <BinnedNLLH.h>

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

TEST_CASE("StatisticSum working with multiple BinnedNLLH instances") {
    /*
     * Scenario for test: two datasets taken, one with background only, another with same background & signal.
     * Essentially an on/off comparison measurement.
     * Truth is a flat background, 12 events in the 3 bins of 1-4 MeV; 4 events of signal in the 2-3 MeV bin only.
     */
    // Set up signal & background PDFs
    AxisCollection ax;
    ax.AddAxis(BinAxis("energy", 1, 4, 3));
    const std::vector<std::string> observables = {"energy"};
    
    BinnedED signal("signal", ax);
    BinnedED background("background", ax);
    
    signal.SetObservables(observables);
    background.SetObservables(observables);
    
    signal.Fill(2.5);
    background.Fill(1.5);
    background.Fill(2.5);
    background.Fill(3.5);
    background.Normalise();
    // Create (fake) datasets
    const double n_signal = 4;
    const double n_back = 12;
    
    BinnedED* data_b = dynamic_cast<BinnedED*>(background.Clone());
    data_b->Scale(n_back);
    
    BinnedED* data_sb = dynamic_cast<BinnedED*>(signal.Clone());
    data_sb->Scale(n_signal);
    data_sb->Add(*data_b);
    // Set up log-likelihood test statistics for each dataset
    BinnedNLLH lh_b;
    lh_b.AddPdf(background);
    lh_b.SetDataDist(*data_b);
    
    BinnedNLLH lh_sb;
    lh_sb.AddPdf(background);
    lh_sb.AddPdf(signal);
    lh_sb.SetDataDist(*data_sb);

    SECTION("StatisticSum setup with BinnedNLLH instances") {
        // Set up StatisticSum object
        std::vector<TestStatistic*> stats = {&lh_b, &lh_sb};
        StatisticSum sum_lh = Sum(stats);
        sum_lh.RegisterFitComponents();
        // Set fit components to true values
        ParameterDict params = {{"signal", 4}, {"background", 12}};
        sum_lh.SetParameters(params);
        
        const double lh_eval = sum_lh.Evaluate();
        const double lh_exp_b = n_back - 3.*(n_back/3.)*log(n_back/3.);
        const double lh_exp_sb = n_signal + n_back - 2.*(n_back/3.)*log(n_back/3.) - (n_signal + n_back/3.)*log(n_signal + n_back/3.);
        const double lh_exp_tot = lh_exp_b + lh_exp_sb;
        REQUIRE(lh_eval == lh_exp_tot);
    }
}
