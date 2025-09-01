#include <catch2/catch_all.hpp>
#include <Minuit.h>
#include <TestStatistic.h>

class DummyStatistic : public TestStatistic
{
public:
    DummyStatistic()
    {
        fParamName = "x";
        fVal = 1.0;
    }

    double Evaluate() override { return fVal; }

    size_t GetParameterCount() const override { return 1; }

    void SetParameters(const ParameterDict &p) override
    {
        fVal = p.at(fParamName);
    }

    ParameterDict GetParameters() const override
    {
        return {{fParamName, fVal}};
    }

    std::set<std::string> GetParameterNames() const override
    {
        return {fParamName};
    }

    void RegisterFitComponents() override {}

    std::string fParamName;
    double fVal;
};

class DummyQuadStatistic : public TestStatistic
{
public:
    DummyQuadStatistic()
    {
        fParamNames = {"a", "b"};
        fParamValues = {1, 1};
    }

    // NLL = (a-2)^2 + (b-3)^2
    double Evaluate() override
    {
        return (fParamValues.at(0) - 2.0) * (fParamValues.at(0) - 2.0) + (fParamValues.at(1) - 3.0) * (fParamValues.at(1) - 3.0);
    }

    size_t GetParameterCount() const override { return 2; }

    void SetParameters(const ParameterDict &p) override
    {
        fParamValues.at(0) = p.at("a");
        fParamValues.at(1) = p.at("b");
    }

    ParameterDict GetParameters() const override
    {
        return {{"a", fParamValues.at(0)}, {"b", fParamValues.at(1)}};
    }

    std::set<std::string> GetParameterNames() const override
    {
        return fParamNames;
    }

    void RegisterFitComponents() override {}

private:
    std::vector<double> fParamValues;
    std::set<std::string> fParamNames;
};


TEST_CASE("Minuit configuration methods")
{
    Minuit min;

    SECTION("Set and get method")
    {
        min.SetMethod("Simplex");
        REQUIRE(min.GetMethod() == "Simplex");
    }

    SECTION("Set and get tolerance")
    {
        min.SetTolerance(0.01);
        REQUIRE(min.GetTolerance() == Catch::Approx(0.01));
    }

    SECTION("Set and get strategy")
    {
        min.SetStrategy(2);
        REQUIRE(min.GetStrategy() == 2);
    }

    SECTION("Set and get max calls")
    {
        min.SetMaxCalls(500);
        REQUIRE(min.GetMaxCalls() == 500);
    }

    SECTION("Set and get upper contour edge")
    {
        min.SetUpperContourEdge(1.23);
        REQUIRE(min.GetUpperContourEdge() == Catch::Approx(1.23));
    }

    SECTION("Fix and release parameters")
    {
        min.Fix("x");
        min.Release("x"); // Just ensure no exceptions or crashes
    }

    SECTION("Set and get minima/maxima")
    {
        ParameterDict minima = {{"x", -1.}};
        ParameterDict maxima = {{"x", 5.}};
        min.SetMinima(minima);
        min.SetMaxima(maxima);
        REQUIRE(min.GetMinima() == minima);
        REQUIRE(min.GetMaxima() == maxima);
    }

    SECTION("Set and get initial values/errors")
    {
        ParameterDict vals = {{"x", 2.}};
        ParameterDict errs = {{"x", 0.5}};
        min.SetInitialValues(vals);
        min.SetInitialErrors(errs);
        // No getters, so just checking that it doesn't crash
    }

    SECTION("Set and get maximising flag")
    {
        min.SetMaximising(true);
        REQUIRE(min.GetMaximising() == true);
    }
}

TEST_CASE("Minuit optimise returns FitResult")
{
    Minuit min;
    DummyStatistic stat;

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

TEST_CASE("Minuit fit result getter returns last fit")
{
    Minuit min;
    DummyStatistic stat;

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

TEST_CASE("Fixing Minuit parameters returns correct values")
{
    Minuit min;
    DummyQuadStatistic stat;

    // Start at a = 10; b = 10
    // Best point should be a = 2; b = 3 if both are free
    ParameterDict vals = {{"a", 10.0}, {"b", 10.0}};
    ParameterDict errs = {{"a", 1.0}, {"b", 1.0}};
    ParameterDict mins = {{"a", -100.0}, {"b", -100.0}};
    ParameterDict maxs = {{"a", 100.0}, {"b", 100.0}};

    min.SetInitialValues(vals);
    min.SetInitialErrors(errs);
    min.SetMinima(mins);
    min.SetMaxima(maxs);
    min.SetMaxCalls(10000);
    min.SetMethod("Migrad");

    min.Optimise(&stat);
    FitResult result = min.GetFitResult();

    ParameterDict best = result.GetBestFit();

    REQUIRE(best.at("a") == Catch::Approx(2.0).margin(1e-3));
    REQUIRE(best.at("b") == Catch::Approx(3.0).margin(1e-3));

    min.Fix("b");

    min.Optimise(&stat);
    result = min.GetFitResult();

    best = result.GetBestFit();

    REQUIRE(best.at("a") == Catch::Approx(2.0).margin(1e-3));
    REQUIRE(best.at("b") == Catch::Approx(10.0));

    min.Release("b");
    min.Fix("a");

    min.Optimise(&stat);
    result = min.GetFitResult();

    best = result.GetBestFit();

    REQUIRE(best.at("a") == Catch::Approx(10.0));
    REQUIRE(best.at("b") == Catch::Approx(3.0).margin(1e-3));
}
