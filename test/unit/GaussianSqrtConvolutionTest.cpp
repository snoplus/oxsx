#include <catch2/catch_all.hpp>
#include <SystematicManager.h>
#include <GaussianSqrtConvolution.h>
#include <JumpPDF.h>
#include <gsl/gsl_cdf.h>

TEST_CASE("Simple GaussianSqrtConvolution systematic on 1d PDF")
{
  // First - build base 1D BinnedED object for applying systematic to
  AxisCollection axes;
  axes.AddAxis(BinAxis("axis0", 0, 4, 4));
  std::vector<std::string> observables;
  observables.push_back("obs0");

  BinnedED pdf1("pdf1", axes);
  pdf1.SetBinContent(1, 10);
  pdf1.SetObservables(observables);

  // Now build the GaussianSqrtConvolution systematic
  double sigma = 0.5;

  GaussianSqrtConvolution conv("smear_sys");
  conv.RenameSigma("sigma");
  conv.SetSigma(sigma);
  conv.SetAxes(axes);
  conv.SetTransformationObs(observables);
  conv.SetDistributionObs(observables);

  SECTION("Check GetSigmaName() method")
  {
    REQUIRE(conv.GetSigmaName() == "sigma");
  }

  SECTION("Check FitComponent interface for GaussianSqrtConvolution, 1D")
  {
    REQUIRE(conv.GetParameterCount() == 1); // one param: Gaussian's sigma
    REQUIRE(conv.GetParameter("sigma") == sigma);
    REQUIRE(conv.GetName() == "smear_sys");

    sigma = 0.4;
    conv.SetParameter("sigma", sigma);
    REQUIRE(conv.GetParameter("sigma") == sigma);
  }

  // Add systematic & BinnedED objects to a systematic manager; test smearing
  SystematicManager man;

  SECTION("Test 1D GaussianSqrtConvolution smearing")
  {
    man.Add(&conv);
    man.AddDist(pdf1, "");
    man.Construct();

    std::vector<BinnedED> pdfs = {pdf1};
    std::vector<BinnedED> OrignalPdfs(pdfs);

    man.DistortEDs(OrignalPdfs, pdfs);

    const double mu = 1.5; // centre of bin which has data in it
    const double s = sigma*std::sqrt(mu);
    std::vector<double> modifiedObs = pdfs.at(0).GetBinContents();
    std::vector<double> correctVals = {
        10. * (gsl_cdf_gaussian_P(1 - mu, s) - gsl_cdf_gaussian_P(0 - mu, s)),
        10. * (gsl_cdf_gaussian_P(2 - mu, s) - gsl_cdf_gaussian_P(1 - mu, s)),
        10. * (gsl_cdf_gaussian_P(3 - mu, s) - gsl_cdf_gaussian_P(2 - mu, s)),
        10. * (gsl_cdf_gaussian_P(4 - mu, s) - gsl_cdf_gaussian_P(3 - mu, s)),
    };

    REQUIRE(modifiedObs == correctVals);
  }
}

TEST_CASE("Simple GaussianSqrtConvolution systematic on 1d PDF, small CDF cutoff")
{
  // First - build base 1D BinnedED object for applying systematic to
  AxisCollection axes;
  axes.AddAxis(BinAxis("axis0", 0, 4, 4));
  std::vector<std::string> observables;
  observables.push_back("obs0");

  BinnedED pdf1("pdf1", axes);
  pdf1.SetBinContent(1, 10);
  pdf1.SetObservables(observables);

  // Now build the GaussianSqrtConvolution systematic
  double sigma = 0.5;
  const double cdfcutoff = 1.; // comedically small cutoff, just to test 
  GaussianSqrtConvolution conv("smear_sys", cdfcutoff);
  conv.RenameSigma("sigma");
  conv.SetSigma(sigma);
  conv.SetAxes(axes);
  conv.SetTransformationObs(observables);
  conv.SetDistributionObs(observables);

  SECTION("Check GetSigmaName() method, with custom cutoff in place")
  {
    REQUIRE(conv.GetSigmaName() == "sigma");
  }

  SECTION("Check FitComponent interface for GaussianSqrtConvolution, 1D, with custom cutoff in place")
  {
    REQUIRE(conv.GetParameterCount() == 1); // one param: Gaussian's sigma
    REQUIRE(conv.GetParameter("sigma") == sigma);
    REQUIRE(conv.GetName() == "smear_sys");

    sigma = 0.4;
    conv.SetParameter("sigma", sigma);
    REQUIRE(conv.GetParameter("sigma") == sigma);
  }

  // Add systematic & BinnedED objects to a systematic manager; test smearing
  SystematicManager man;

  SECTION("Test 1D GaussianSqrtConvolution smearing, with custom cutoff in place")
  {
    man.Add(&conv);
    man.AddDist(pdf1, "");
    man.Construct();

    std::vector<BinnedED> pdfs = {pdf1};
    std::vector<BinnedED> OrignalPdfs(pdfs);

    man.DistortEDs(OrignalPdfs, pdfs);

    const double mu = 1.5; // centre of bin which has data in it
    const double s = sigma*std::sqrt(mu);
    std::vector<double> modifiedObs = pdfs.at(0).GetBinContents();
    // values different now because of the CDF
    std::vector<double> correctVals = {
        10. * (gsl_cdf_gaussian_P(1 - mu, s) - 0.),
        10. * (gsl_cdf_gaussian_P(2 - mu, s) - gsl_cdf_gaussian_P(1 - mu, s)),
        10. * (1. - gsl_cdf_gaussian_P(2 - mu, s)),
        10. * (0. - 0.),
    };

    REQUIRE(modifiedObs == correctVals);
  }
}

TEST_CASE("Simple GaussianSqrtConvolution systematic on 2d PDF")
{
  // First - build base 1D BinnedED object for applying systematic to
  AxisCollection axes;
  axes.AddAxis(BinAxis("axis0", 0, 4, 4));
  axes.AddAxis(BinAxis("axis1", 0, 2, 2));
  const std::vector<std::string> observables = {"obs0", "obs1"};
  const std::vector<std::string> obs_smear = {"obs0"};

  BinnedED pdf1("pdf1", axes);
  pdf1.SetBinContent(axes.FlattenIndices({1, 0}), 10);
  pdf1.SetBinContent(axes.FlattenIndices({2, 1}), 20);
  pdf1.SetObservables(observables);

  // Now build the GaussianSqrtConvolution systematic and related objects
  double sigma = 0.5;

  GaussianSqrtConvolution conv("smear_sys");
  conv.RenameSigma("sigma");
  conv.SetSigma(sigma);
  conv.SetAxes(axes);
  conv.SetTransformationObs(obs_smear);
  conv.SetDistributionObs(observables);

  SECTION("Check FitComponent interface for GaussianSqrtConvolution, 2D")
  {
    REQUIRE(conv.GetParameterCount() == 1); // one param: Gaussian's sigma
    REQUIRE(conv.GetParameter("sigma") == sigma);
    REQUIRE(conv.GetName() == "smear_sys");

    sigma = 0.4;
    conv.SetParameter("sigma", sigma);
    REQUIRE(conv.GetParameter("sigma") == sigma);
  }

  // Add systematic & BinnedED objects to a systematic manager; test smearing
  SystematicManager man;

  SECTION("Test 2D GaussianSqrtConvolution smearing")
  {
    man.Add(&conv);
    man.AddDist(pdf1, "");
    man.Construct();

    std::vector<BinnedED> pdfs = {pdf1};
    std::vector<BinnedED> OrignalPdfs(pdfs);

    man.DistortEDs(OrignalPdfs, pdfs);

    const double mu1 = 1.5; // centres of bins which has data in it
    const double mu2 = 2.5; // centres of bins which has data in it
    const double s1 = sigma*std::sqrt(mu1);
    const double s2 = sigma*std::sqrt(mu2);
    std::vector<double> modifiedObs = pdfs.at(0).GetBinContents();
    std::vector<double> correctVals(pdf1.GetNBins(), 0);
    correctVals[axes.FlattenIndices({0, 0})] = 10. * (gsl_cdf_gaussian_P(1 - mu1, s1) - gsl_cdf_gaussian_P(0 - mu1, s1));
    correctVals[axes.FlattenIndices({1, 0})] = 10. * (gsl_cdf_gaussian_P(2 - mu1, s1) - gsl_cdf_gaussian_P(1 - mu1, s1));
    correctVals[axes.FlattenIndices({2, 0})] = 10. * (gsl_cdf_gaussian_P(3 - mu1, s1) - gsl_cdf_gaussian_P(2 - mu1, s1));
    correctVals[axes.FlattenIndices({3, 0})] = 10. * (gsl_cdf_gaussian_P(4 - mu1, s1) - gsl_cdf_gaussian_P(3 - mu1, s1));
    correctVals[axes.FlattenIndices({0, 1})] = 20. * (gsl_cdf_gaussian_P(1 - mu2, s2) - gsl_cdf_gaussian_P(0 - mu2, s2));
    correctVals[axes.FlattenIndices({1, 1})] = 20. * (gsl_cdf_gaussian_P(2 - mu2, s2) - gsl_cdf_gaussian_P(1 - mu2, s2));
    correctVals[axes.FlattenIndices({2, 1})] = 20. * (gsl_cdf_gaussian_P(3 - mu2, s2) - gsl_cdf_gaussian_P(2 - mu2, s2));
    correctVals[axes.FlattenIndices({3, 1})] = 20. * (gsl_cdf_gaussian_P(4 - mu2, s2) - gsl_cdf_gaussian_P(3 - mu2, s2));

    REQUIRE(modifiedObs == correctVals);
  }
}

TEST_CASE("Simple GaussianSqrtConvolution systematic on 2d PDF, alt axis ordering")
{
  // First - build base 1D BinnedED object for applying systematic to
  AxisCollection axes;
  axes.AddAxis(BinAxis("axis1", 0, 2, 2));
  axes.AddAxis(BinAxis("axis0", 0, 4, 4));
  const std::vector<std::string> observables = {"obs1", "obs0"};
  const std::vector<std::string> obs_smear = {"obs0"};

  BinnedED pdf1("pdf1", axes);
  pdf1.SetBinContent(axes.FlattenIndices({0, 1}), 10);
  pdf1.SetBinContent(axes.FlattenIndices({1, 2}), 20);
  pdf1.SetObservables(observables);

  // Now build the GaussianSqrtConvolution systematic and related objects
  double sigma = 0.5;

  GaussianSqrtConvolution conv("smear_sys");
  conv.RenameSigma("sigma");
  conv.SetSigma(sigma);
  conv.SetAxes(axes);
  conv.SetTransformationObs(obs_smear);
  conv.SetDistributionObs(observables);

  SECTION("Check FitComponent interface for GaussianSqrtConvolution, 2D")
  {
    REQUIRE(conv.GetParameterCount() == 1); // one params: Gaussian's sigma
    REQUIRE(conv.GetParameter("sigma") == sigma);
    REQUIRE(conv.GetName() == "smear_sys");

    sigma = 0.4;
    conv.SetParameter("sigma", sigma);
    REQUIRE(conv.GetParameter("sigma") == sigma);
  }

  // Add systematic & BinnedED objects to a systematic manager; test smearing
  SystematicManager man;

  SECTION("Test 2D GaussianSqrtConvolution smearing")
  {
    man.Add(&conv);
    man.AddDist(pdf1, "");
    man.Construct();

    std::vector<BinnedED> pdfs = {pdf1};
    std::vector<BinnedED> OrignalPdfs(pdfs);

    man.DistortEDs(OrignalPdfs, pdfs);

    const double mu1 = 1.5; // centres of bins which has data in it
    const double mu2 = 2.5; // centres of bins which has data in it
    const double s1 = sigma*std::sqrt(mu1);
    const double s2 = sigma*std::sqrt(mu2);
    std::vector<double> modifiedObs = pdfs.at(0).GetBinContents();
    std::vector<double> correctVals(pdf1.GetNBins(), 0);
    correctVals[axes.FlattenIndices({0, 0})] = 10. * (gsl_cdf_gaussian_P(1 - mu1, s1) - gsl_cdf_gaussian_P(0 - mu1, s1));
    correctVals[axes.FlattenIndices({0, 1})] = 10. * (gsl_cdf_gaussian_P(2 - mu1, s1) - gsl_cdf_gaussian_P(1 - mu1, s1));
    correctVals[axes.FlattenIndices({0, 2})] = 10. * (gsl_cdf_gaussian_P(3 - mu1, s1) - gsl_cdf_gaussian_P(2 - mu1, s1));
    correctVals[axes.FlattenIndices({0, 3})] = 10. * (gsl_cdf_gaussian_P(4 - mu1, s1) - gsl_cdf_gaussian_P(3 - mu1, s1));
    correctVals[axes.FlattenIndices({1, 0})] = 20. * (gsl_cdf_gaussian_P(1 - mu2, s2) - gsl_cdf_gaussian_P(0 - mu2, s2));
    correctVals[axes.FlattenIndices({1, 1})] = 20. * (gsl_cdf_gaussian_P(2 - mu2, s2) - gsl_cdf_gaussian_P(1 - mu2, s2));
    correctVals[axes.FlattenIndices({1, 2})] = 20. * (gsl_cdf_gaussian_P(3 - mu2, s2) - gsl_cdf_gaussian_P(2 - mu2, s2));
    correctVals[axes.FlattenIndices({1, 3})] = 20. * (gsl_cdf_gaussian_P(4 - mu2, s2) - gsl_cdf_gaussian_P(3 - mu2, s2));

    REQUIRE(modifiedObs == correctVals);
  }
}

TEST_CASE("GaussianSqrtConvolution systematic on 1d PDF, variable binning")
{
  // First - build base 1D BinnedED object for applying systematic to
  AxisCollection axes;
  const std::vector<double> lowEdges {0., 1., 3., 4.};
  const std::vector<double> highEdges {1., 3., 4., 4.5};
  axes.AddAxis(BinAxis("axis0", lowEdges, highEdges));
  std::vector<std::string> observables;
  observables.push_back("obs0");

  BinnedED pdf1("pdf1", axes);
  pdf1.SetBinContent(1, 10);
  pdf1.SetObservables(observables);

  // Now build the GaussianSqrtConvolution systematic
  double sigma = 0.5;

  GaussianSqrtConvolution conv("smear_sys");
  conv.RenameSigma("sigma");
  conv.SetSigma(sigma);
  conv.SetAxes(axes);
  conv.SetTransformationObs(observables);
  conv.SetDistributionObs(observables);

  SECTION("Check GetSigmaName() method")
  {
    REQUIRE(conv.GetSigmaName() == "sigma");
  }

  SECTION("Check FitComponent interface for GaussianSqrtConvolution, 1D")
  {
    REQUIRE(conv.GetParameterCount() == 1); // one param: Gaussian's sigma
    REQUIRE(conv.GetParameter("sigma") == sigma);
    REQUIRE(conv.GetName() == "smear_sys");

    sigma = 1.;
    conv.SetParameter("sigma", sigma);
    REQUIRE(conv.GetParameter("sigma") == sigma);
  }

  // Add systematic & BinnedED objects to a systematic manager; test smearing
  SystematicManager man;

  SECTION("Test 1D GaussianSqrtConvolution smearing")
  {
    man.Add(&conv);
    man.AddDist(pdf1, "");
    man.Construct();

    std::vector<BinnedED> pdfs = {pdf1};
    std::vector<BinnedED> OrignalPdfs(pdfs);

    man.DistortEDs(OrignalPdfs, pdfs);

    const double mu = 2.; // centre of bin which has data in it
    const double s = sigma*std::sqrt(mu);
    std::vector<double> modifiedObs = pdfs.at(0).GetBinContents();
    std::vector<double> correctVals = {
        10. * (gsl_cdf_gaussian_P(1 - mu, s) - gsl_cdf_gaussian_P(0 - mu, s)),
        10. * (gsl_cdf_gaussian_P(3 - mu, s) - gsl_cdf_gaussian_P(1 - mu, s)),
        10. * (gsl_cdf_gaussian_P(4 - mu, s) - gsl_cdf_gaussian_P(3 - mu, s)),
        10. * (gsl_cdf_gaussian_P(4.5 - mu, s) - gsl_cdf_gaussian_P(4 - mu, s)),
    };

    REQUIRE(modifiedObs == correctVals);
  }
}