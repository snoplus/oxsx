#include <catch2/catch_all.hpp>
#include <StatisticSum.h>
#include <BinnedNLLH.h>

class FakeStatistic : public TestStatistic
{
public:
    double Evaluate()
    {
        return fVal;
    }

    size_t GetParameterCount() const
    {
        return 1;
    }

    void SetParameters(const ParameterDict &p)
    {
        fVal = p.at(fParamName);
    }

    ParameterDict GetParameters() const
    {
        ParameterDict p;
        p[fParamName] = fVal;
        return p;
    }

    std::set<std::string> GetParameterNames() const
    {
        std::set<std::string> set;
        set.insert(fParamName);
        return set;
    }

    void RegisterFitComponents() {}

    double fVal;
    std::string fParamName;
};

TEST_CASE("Adding test statistics using StatisticSum constr")
{
    FakeStatistic s1;
    FakeStatistic s2;

    s1.fVal = 1;
    s1.fParamName = "p1";
    s2.fVal = 2;
    s2.fParamName = "p2";

    SECTION("no shared parameters")
    {
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

    SECTION("setting with shared parameters")
    {
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

    SECTION("adding test stats with +")
    {
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

    SECTION("Test stat sum operator")
    {
        FakeStatistic s3;
        s3.fVal = 3;
        s3.fParamName = "p3";

        std::vector<TestStatistic *> stats;
        stats.push_back(&s1);
        stats.push_back(&s2);
        stats.push_back(&s3);

        StatisticSum sum = Sum(stats);
        REQUIRE(sum.Evaluate() == 6);
    }
}

TEST_CASE("StatisticSum working with multiple BinnedNLLH instances")
{
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

    BinnedED *data_b = dynamic_cast<BinnedED *>(background.Clone());
    data_b->Scale(n_back);

    BinnedED *data_sb = dynamic_cast<BinnedED *>(signal.Clone());
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

    SECTION("StatisticSum setup with BinnedNLLH instances")
    {
        // Set up StatisticSum object
        std::vector<TestStatistic *> stats = {&lh_b, &lh_sb};
        StatisticSum sum_lh = Sum(stats);
        sum_lh.RegisterFitComponents();
        // Set fit components to true values
        ParameterDict params = {{"signal", 4}, {"background", 12}};
        sum_lh.SetParameters(params);

        double lh_eval = sum_lh.Evaluate();
        const double lh_exp_b = n_back - 3. * (n_back / 3.) * log(n_back / 3.);
        const double lh_exp_sb = n_signal + n_back - 2. * (n_back / 3.) * log(n_back / 3.) - (n_signal + n_back / 3.) * log(n_signal + n_back / 3.);
        const double lh_exp_tot = lh_exp_b + lh_exp_sb;
        REQUIRE(lh_eval == lh_exp_tot);

        double constr_mean_signal = 2;
        double constr_uncert_signal = 0.4;
        double constr_mean_bg = 10;
        double constr_uncert_bg = 0.12;
        double corr_factor = 0.1;
        sum_lh.SetConstraint("signal", constr_mean_signal, constr_uncert_signal, "background", constr_mean_bg, constr_uncert_bg, corr_factor);

        sum_lh.SetParameters(params);
        sum_lh.RegisterFitComponents();

        const double z1 = (n_signal - constr_mean_signal) * (n_signal - constr_mean_signal) / (2 * constr_uncert_signal * constr_uncert_signal);
        const double z2 = (n_back - constr_mean_bg) * (n_back - constr_mean_bg) / (2 * constr_uncert_bg * constr_uncert_bg);
        const double z12 = corr_factor * (n_signal - constr_mean_signal) * (n_back - constr_mean_bg) / (constr_uncert_signal * constr_uncert_bg);
        double lh_constr = (z1 - z12 + z2) / (1. - corr_factor * corr_factor);

        lh_eval = sum_lh.Evaluate();
        REQUIRE(lh_eval == lh_exp_tot + lh_constr);
    }

    SECTION("All Constraints Set in StatisticSum")
    {

        BinnedED signal1("signal1", ax);
        BinnedED background1("background1", ax);
        BinnedED signal2("signal2", ax);
        BinnedED background2("background2", ax);

        signal1.SetObservables(observables);
        background1.SetObservables(observables);
        signal2.SetObservables(observables);
        background2.SetObservables(observables);

        signal1.Fill(2.5);
        signal2.Fill(1.5);
        background1.Fill(1.5);
        background1.Fill(2.5);
        background1.Fill(3.5);
        background2.Fill(1.5);
        background2.Fill(2.5);
        background2.Fill(3.5);
        background1.Normalise();
        background2.Normalise();

        // Create (fake) datasets
        const double n_signal1 = 4;
        const double n_back1 = 12;
        const double n_signal2 = 2;
        const double n_back2 = 6;

        BinnedED *data_b1 = dynamic_cast<BinnedED *>(background1.Clone());
        data_b1->Scale(n_back1);
        BinnedED *data_b2 = dynamic_cast<BinnedED *>(background2.Clone());
        data_b2->Scale(n_back2);

        BinnedED *data_sb1 = dynamic_cast<BinnedED *>(signal1.Clone());
        data_sb1->Scale(n_signal1);
        data_sb1->Add(*data_b1);

        BinnedED *data_sb2 = dynamic_cast<BinnedED *>(signal2.Clone());
        data_sb2->Scale(n_signal2);
        data_sb2->Add(*data_b2);

        // Set up log-likelihood test statistics for each dataset
        BinnedNLLH lh1;
        lh1.AddPdf(background1);
        lh1.AddPdf(signal1);
        lh1.SetDataDist(*data_sb1);

        BinnedNLLH lh2;
        lh2.AddPdf(background2);
        lh2.AddPdf(signal2);
        lh2.SetDataDist(*data_sb2);

        // Set up StatisticSum object
        std::vector<TestStatistic *> stats = {&lh1, &lh2};
        StatisticSum sum_lh = Sum(stats);
        sum_lh.RegisterFitComponents();
        // Set fit components to true values
        ParameterDict params = {{"signal1", 4}, {"background1", 12}, {"signal2", 2}, {"background2", 6}};
        sum_lh.SetParameters(params);

        double constr_mean_bg1 = 10;
        double constr_uncert_bg1 = 0.4;
        double constr_mean_bg2 = 5;
        double constr_uncert_bg2 = 0.12;
        double corr_factor = 0.1;
        double constr_mean_signal1 = 3;
        double constr_uncert_signal1 = 1.0;
        sum_lh.SetConstraint("background1", constr_mean_bg1, constr_uncert_bg1, "background2", constr_mean_bg2, constr_uncert_bg2, corr_factor);
        sum_lh.SetConstraint("signal1", constr_mean_signal1, constr_uncert_signal1);

        sum_lh.RegisterFitComponents();

        double lh_eval = sum_lh.Evaluate();
        const double lh_exp_1 = n_signal1 + n_back1 - 2. * (n_back1 / 3.) * log(n_back1 / 3.) - (n_signal1 + n_back1 / 3.) * log(n_signal1 + n_back1 / 3.);
        const double lh_exp_2 = n_signal2 + n_back2 - 2. * (n_back2 / 3.) * log(n_back2 / 3.) - (n_signal2 + n_back2 / 3.) * log(n_signal2 + n_back2 / 3.);

        const double z1 = (n_back1 - constr_mean_bg1) * (n_back1 - constr_mean_bg1) / (2 * constr_uncert_bg1 * constr_uncert_bg1);
        const double z2 = (n_back2 - constr_mean_bg2) * (n_back2 - constr_mean_bg2) / (2 * constr_uncert_bg2 * constr_uncert_bg2);
        const double z12 = corr_factor * (n_back1 - constr_mean_bg1) * (n_back2 - constr_mean_bg2) / (constr_uncert_bg1 * constr_uncert_bg2);
        double lh_constr = (z1 - z12 + z2) / (1. - corr_factor * corr_factor);

        lh_constr += (n_signal1 - constr_mean_signal1) * (n_signal1 - constr_mean_signal1) / (2 * constr_uncert_signal1 * constr_uncert_signal1);

        const double lh_exp_tot = lh_exp_1 + lh_exp_2 + lh_constr;
        REQUIRE(lh_eval == Catch::Approx(lh_exp_tot).margin(1e-7));
    }

    SECTION("Constraints Set in Inidividual Component TestStatistics")
    {

        BinnedED signal1("signal1", ax);
        BinnedED background1("background1", ax);
        BinnedED signal2("signal2", ax);
        BinnedED background2("background2", ax);

        signal1.SetObservables(observables);
        background1.SetObservables(observables);
        signal2.SetObservables(observables);
        background2.SetObservables(observables);

        signal1.Fill(2.5);
        signal2.Fill(1.5);
        background1.Fill(1.5);
        background1.Fill(2.5);
        background1.Fill(3.5);
        background2.Fill(1.5);
        background2.Fill(2.5);
        background2.Fill(3.5);
        background1.Normalise();
        background2.Normalise();

        // Create (fake) datasets
        const double n_signal1 = 4;
        const double n_back1 = 12;
        const double n_signal2 = 2;
        const double n_back2 = 6;

        BinnedED *data_b1 = dynamic_cast<BinnedED *>(background1.Clone());
        data_b1->Scale(n_back1);
        BinnedED *data_b2 = dynamic_cast<BinnedED *>(background2.Clone());
        data_b2->Scale(n_back2);

        BinnedED *data_sb1 = dynamic_cast<BinnedED *>(signal1.Clone());
        data_sb1->Scale(n_signal1);
        data_sb1->Add(*data_b1);

        BinnedED *data_sb2 = dynamic_cast<BinnedED *>(signal2.Clone());
        data_sb2->Scale(n_signal2);
        data_sb2->Add(*data_b2);

        // Set up log-likelihood test statistics for each dataset
        BinnedNLLH lh1;
        lh1.AddPdf(background1);
        lh1.AddPdf(signal1);
        lh1.SetDataDist(*data_sb1);
        double constr_mean_signal1 = 3;
        double constr_uncert_signal1 = 1.0;
        lh1.SetConstraint("signal1", constr_mean_signal1, constr_uncert_signal1);
        lh1.RegisterFitComponents();

        BinnedNLLH lh2;
        lh2.AddPdf(background2);
        lh2.AddPdf(signal2);
        lh2.SetDataDist(*data_sb2);
        lh2.RegisterFitComponents();

        // Set up StatisticSum object
        std::vector<TestStatistic *> stats = {&lh1, &lh2};
        // Set fit components to true values in individual LLHs
        ParameterDict params1 = {{"signal1", 4}, {"background1", 12}};
        ParameterDict params2 = {{"signal2", 2}, {"background2", 6}};
        lh1.SetParameters(params1);
        lh2.SetParameters(params2);

        StatisticSum sum_lh = Sum(stats);

        double constr_mean_bg1 = 10;
        double constr_uncert_bg1 = 0.4;
        double constr_mean_bg2 = 5;
        double constr_uncert_bg2 = 0.12;
        double corr_factor = 0.1;
        sum_lh.SetConstraint("background1", constr_mean_bg1, constr_uncert_bg1, "background2", constr_mean_bg2, constr_uncert_bg2, corr_factor);

        sum_lh.RegisterFitComponents();

        double lh_eval = sum_lh.Evaluate();
        const double lh_exp_1 = n_signal1 + n_back1 - 2. * (n_back1 / 3.) * log(n_back1 / 3.) - (n_signal1 + n_back1 / 3.) * log(n_signal1 + n_back1 / 3.);
        const double lh_exp_2 = n_signal2 + n_back2 - 2. * (n_back2 / 3.) * log(n_back2 / 3.) - (n_signal2 + n_back2 / 3.) * log(n_signal2 + n_back2 / 3.);

        const double z1 = (n_back1 - constr_mean_bg1) * (n_back1 - constr_mean_bg1) / (2 * constr_uncert_bg1 * constr_uncert_bg1);
        const double z2 = (n_back2 - constr_mean_bg2) * (n_back2 - constr_mean_bg2) / (2 * constr_uncert_bg2 * constr_uncert_bg2);
        const double z12 = corr_factor * (n_back1 - constr_mean_bg1) * (n_back2 - constr_mean_bg2) / (constr_uncert_bg1 * constr_uncert_bg2);
        double lh_constr = (z1 - z12 + z2) / (1. - corr_factor * corr_factor);

        lh_constr += (n_signal1 - constr_mean_signal1) * (n_signal1 - constr_mean_signal1) / (2 * constr_uncert_signal1 * constr_uncert_signal1);

        const double lh_exp_tot = lh_exp_1 + lh_exp_2 + lh_constr;
        REQUIRE(lh_eval == Catch::Approx(lh_exp_tot).margin(1e-7));
    }
}
