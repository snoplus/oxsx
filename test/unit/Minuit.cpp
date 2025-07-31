#include <catch2/catch_all.hpp>
#include <Minuit.h>
#include <TestStatistic.h>

class FakeStatistic : public TestStatistic {
public:
    FakeStatistic() {
        fParamName = "x";
        fVal = 1.0;
    }

    double Evaluate() override { return fVal; }

    size_t GetParameterCount() const override { return 1; }

    void SetParameters(const ParameterDict &p) override {
        fVal = p.at(fParamName);
    }

    ParameterDict GetParameters() const override {
        return {{fParamName, fVal}};
    }

    std::set<std::string> GetParameterNames() const override {
        return {fParamName};
    }

    void RegisterFitComponents() override {}

    std::string fParamName;
    double fVal;
};

TEST_CASE("Minuit configuration methods") {
    Minuit min;

    SECTION("Set and get method") {
        min.SetMethod("Simplex");
        REQUIRE(min.GetMethod() == "Simplex");
    }

    SECTION("Set and get tolerance") {
        min.SetTolerance(0.01);
        REQUIRE(min.GetTolerance() == Catch::Approx(0.01));
    }

    SECTION("Set and get strategy") {
        min.SetStrategy(2);
        REQUIRE(min.GetStrategy() == 2);
    }

    SECTION("Set and get max calls") {
        min.SetMaxCalls(500);
        REQUIRE(min.GetMaxCalls() == 500);
    }

    SECTION("Set and get upper contour edge") {
        min.SetUpperContourEdge(1.23);
        REQUIRE(min.GetUpperContourEdge() == Catch::Approx(1.23));
    }

    SECTION("Fix and release parameters") {
        min.Fix("x");
        min.Release("x"); // Just ensure no exceptions or crashes
    }

    SECTION("Set and get minima/maxima") {
        ParameterDict minima = {{"x", -1.}};
        ParameterDict maxima = {{"x", 5.}};
        min.SetMinima(minima);
        min.SetMaxima(maxima);
        REQUIRE(min.GetMinima() == minima);
        REQUIRE(min.GetMaxima() == maxima);
    }

    SECTION("Set and get initial values/errors") {
        ParameterDict vals = {{"x", 2.}};
        ParameterDict errs = {{"x", 0.5}};
        min.SetInitialValues(vals);
        min.SetInitialErrors(errs);
        // No getters, so just checking that it doesn't crash
    }

    SECTION("Set and get maximising flag") {
        min.SetMaximising(true);
        REQUIRE(min.GetMaximising() == true);
    }
}

TEST_CASE("Minuit optimise returns FitResult") {
    Minuit min;
    FakeStatistic stat;

    ParameterDict vals = {{"x", 1.}};
    ParameterDict errs = {{"x", 0.1}};
    ParameterDict mins = {{"x", 0.}};
    ParameterDict maxs = {{"x", 5.}};

    min.SetInitialValues(vals);
    min.SetInitialErrors(errs);
    min.SetMinima(mins);
    min.SetMaxima(maxs);
    min.SetMaxCalls(100);
    min.SetStrategy(1);
    min.SetTolerance(0.01);
    min.SetMethod("Simplex");

    const FitResult &result = min.Optimise(&stat);

    REQUIRE(result.GetValid());
    REQUIRE(result.GetBestFit().count("x") == 1);
}

TEST_CASE("Minuit fit result getter returns last fit") {
    Minuit min;
    FakeStatistic stat;

    ParameterDict vals = {{"x", 2.}};
    ParameterDict errs = {{"x", 0.2}};
    ParameterDict mins = {{"x", 1.}};
    ParameterDict maxs = {{"x", 3.}};

    min.SetInitialValues(vals);
    min.SetInitialErrors(errs);
    min.SetMinima(mins);
    min.SetMaxima(maxs);
    min.SetMaxCalls(100);
    min.SetMethod("Simplex");

    min.Optimise(&stat);
    FitResult result = min.GetFitResult();
    REQUIRE(result.GetBestFit().count("x") == 1);
}
