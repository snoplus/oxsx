#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>
#include <BinnedNLLH.h>
#include <Gaussian.h>
#include <DistTools.h>
#include <OXSXDataSet.h>
#include <Shift.h>
#include <iostream>

TEST_CASE("Binned NLLH, 3 rates no systematics")
{
  Gaussian gaus1(1, 5);
  Gaussian gaus2(1, 8);
  Gaussian gaus3(1, 10);

  AxisCollection axes;
  axes.AddAxis(BinAxis("axis1", -40, 40, 200));

  BinnedED pdf1("a", DistTools::ToHist(gaus1, axes));
  BinnedED pdf2("b", DistTools::ToHist(gaus2, axes));
  BinnedED pdf3("c", DistTools::ToHist(gaus3, axes));

  size_t centralBin = pdf1.FindBin(std::vector<double>(1, 0));
  double prob1 = pdf1.GetBinContent(centralBin);
  double prob2 = pdf2.GetBinContent(centralBin);
  double prob3 = pdf3.GetBinContent(centralBin);

  std::vector<std::string> observable;
  observable.push_back("obs0");
  pdf1.SetObservables(observable);
  pdf2.SetObservables(observable);
  pdf3.SetObservables(observable);

  BinnedNLLH lh;
  lh.AddPdf(pdf1);
  lh.AddPdf(pdf2);
  lh.AddPdf(pdf3);

  OXSXDataSet data;
  data.AddEntry(Event(std::vector<double>(1, 0)));
  data.SetObservableNames(observable);

  lh.SetDataSet(&data);

  lh.RegisterFitComponents();
  SECTION("Correct Probability")
  {
    double sumLogProb = -log(prob1 + prob2 + prob3);
    double sumNorm = 3;

    ParameterDict params;
    params["a"] = 1;
    params["b"] = 1;
    params["c"] = 1;
    lh.SetParameters(params);
    REQUIRE(lh.GetParameters() == params);
    REQUIRE_THAT(lh.Evaluate(), Catch::Matchers::WithinAbs(sumNorm + sumLogProb, 0.0001));
  }
  SECTION("Correct Probability with constraint")
  {
    lh.SetConstraint("a", 3, 1);

    double sumLogProb = -log(prob1 + prob2 + prob3);
    double sumNorm = 3;
    double constraint = 2;

    ParameterDict params;
    params["a"] = 1;
    params["b"] = 1;
    params["c"] = 1;
    lh.SetParameters(params);
    REQUIRE_THAT(lh.Evaluate(), Catch::Matchers::WithinAbs(sumNorm + sumLogProb + constraint, 0.0001));
  }
  SECTION("Correct Probability with asymmetric constraint")
  {
    lh.SetConstraint("b", 5, 1, 2);

    double sumLogProb = -log(prob1 + prob2 + prob3);
    double sumNorm = 3;
    double constraint = 8;

    ParameterDict params;
    params["a"] = 1;
    params["b"] = 1;
    params["c"] = 1;
    lh.SetParameters(params);
    REQUIRE_THAT(lh.Evaluate(), Catch::Matchers::WithinAbs(sumNorm + sumLogProb + constraint, 0.0001));
  }
  SECTION("Correct Probability with asymmetric constraint 2")
  {
    lh.SetConstraint("b", -3, 1, 2);

    double sumLogProb = -log(prob1 + prob2 + prob3);
    double sumNorm = 3;
    double constraint = 2;

    ParameterDict params;
    params["a"] = 1;
    params["b"] = 1;
    params["c"] = 1;
    lh.SetParameters(params);
    REQUIRE_THAT(lh.Evaluate(), Catch::Matchers::WithinAbs(sumNorm + sumLogProb + constraint, 0.0001));
  }

  std::vector<int> genRates(3, pow(10, 6));
  BinnedNLLH lh2;
  lh2.SetBarlowBeeston(true);
  lh2.AddPdf(pdf1, genRates.at(0));
  lh2.AddPdf(pdf2, genRates.at(1));
  lh2.AddPdf(pdf3, genRates.at(2));
  lh2.SetDataSet(&data);
  lh2.RegisterFitComponents();

  SECTION("Correct Probability with Barlow Beeston")
  {
    double betaPen = 0;
    double sumLogProb = 0;
    double sumNorm = 0;

    for (unsigned int i = 0; i < pdf1.GetNBins(); i++)
    {

      double binprob1 = pdf1.GetBinContent(i);
      double binprob2 = pdf2.GetBinContent(i);
      double binprob3 = pdf3.GetBinContent(i);
      double binprob = binprob1 + binprob2 + binprob3;

      double dat = 0;
      if (i == centralBin)
        dat = 1;

      double sig1 = sqrt(1 / (float)genRates.at(0));
      double sig2 = sqrt(1 / (float)genRates.at(1));
      double sig3 = sqrt(1 / (float)genRates.at(2));
      double sig = sqrt(sig1 * sig1 + sig2 * sig2 + sig3 * sig3) / binprob;
      double beta = (-(binprob * sig * sig - 1) + sqrt((binprob * sig * sig - 1) * (binprob * sig * sig - 1) + 4 * dat * sig * sig)) / 2;

      betaPen += (beta - 1) * (beta - 1) / (2 * sig * sig);
      sumLogProb += -dat * log(beta * (binprob));
      sumNorm += beta * (binprob);
    }

    ParameterDict params;
    params["a"] = 1;
    params["b"] = 1;
    params["c"] = 1;
    lh2.SetParameters(params);
    REQUIRE_THAT(lh2.Evaluate(), Catch::Matchers::WithinAbs(sumNorm + sumLogProb + betaPen, 0.0001));
  }
  SECTION("Correct Probability with Barlow Beeston and constraint")
  {
    lh2.SetConstraint("a", 3, 1);

    double betaPen = 0;
    double sumLogProb = 0;
    double sumNorm = 0;
    double constraint = 2;

    for (unsigned int i = 0; i < pdf1.GetNBins(); i++)
    {

      double binprob1 = pdf1.GetBinContent(i);
      double binprob2 = pdf2.GetBinContent(i);
      double binprob3 = pdf3.GetBinContent(i);
      double binprob = binprob1 + binprob2 + binprob3;

      double dat = 0;
      if (i == centralBin)
        dat = 1;

      double sig1 = sqrt(1 / (float)genRates.at(0));
      double sig2 = sqrt(1 / (float)genRates.at(1));
      double sig3 = sqrt(1 / (float)genRates.at(2));
      double sig = sqrt(sig1 * sig1 + sig2 * sig2 + sig3 * sig3) / binprob;
      double beta = (-(binprob * sig * sig - 1) + sqrt((binprob * sig * sig - 1) * (binprob * sig * sig - 1) + 4 * dat * sig * sig)) / 2;

      betaPen += (beta - 1) * (beta - 1) / (2 * sig * sig);
      sumLogProb += -dat * log(beta * (binprob));
      sumNorm += beta * (binprob);
    }

    ParameterDict params;
    params["a"] = 1;
    params["b"] = 1;
    params["c"] = 1;
    lh2.SetParameters(params);
    REQUIRE_THAT(lh2.Evaluate(), Catch::Matchers::WithinAbs(sumNorm + sumLogProb + betaPen + constraint, 0.0001));
  }
  SECTION("Correct probability with Barlow Beeston and asymmetric constraint")
  {
    lh2.SetConstraint("b", 5, 1, 2);

    double betaPen = 0;
    double sumLogProb = 0;
    double sumNorm = 0;
    double constraint = 8;

    for (unsigned int i = 0; i < pdf1.GetNBins(); i++)
    {

      double binprob1 = pdf1.GetBinContent(i);
      double binprob2 = pdf2.GetBinContent(i);
      double binprob3 = pdf3.GetBinContent(i);
      double binprob = binprob1 + binprob2 + binprob3;

      double dat = 0;
      if (i == centralBin)
        dat = 1;

      double sig1 = sqrt(1 / (float)genRates.at(0));
      double sig2 = sqrt(1 / (float)genRates.at(1));
      double sig3 = sqrt(1 / (float)genRates.at(2));
      double sig = sqrt(sig1 * sig1 + sig2 * sig2 + sig3 * sig3) / binprob;
      double beta = (-(binprob * sig * sig - 1) + sqrt((binprob * sig * sig - 1) * (binprob * sig * sig - 1) + 4 * dat * sig * sig)) / 2;

      betaPen += (beta - 1) * (beta - 1) / (2 * sig * sig);
      sumLogProb += -dat * log(beta * (binprob));
      sumNorm += beta * (binprob);
    }

    ParameterDict params;
    params["a"] = 1;
    params["b"] = 1;
    params["c"] = 1;
    lh2.SetParameters(params);
    REQUIRE_THAT(lh2.Evaluate(), Catch::Matchers::WithinAbs(sumNorm + sumLogProb + betaPen + constraint, 0.0001));
  }
  SECTION("Correct Probability with Barlow Beeston and constraint 2")
  {
    lh2.SetConstraint("b", -3, 1, 2);

    double betaPen = 0;
    double sumLogProb = 0;
    double sumNorm = 0;
    double constraint = 2;

    for (unsigned int i = 0; i < pdf1.GetNBins(); i++)
    {

      double binprob1 = pdf1.GetBinContent(i);
      double binprob2 = pdf2.GetBinContent(i);
      double binprob3 = pdf3.GetBinContent(i);
      double binprob = binprob1 + binprob2 + binprob3;

      double dat = 0;
      if (i == centralBin)
        dat = 1;

      double sig1 = sqrt(1 / (float)genRates.at(0));
      double sig2 = sqrt(1 / (float)genRates.at(1));
      double sig3 = sqrt(1 / (float)genRates.at(2));
      double sig = sqrt(sig1 * sig1 + sig2 * sig2 + sig3 * sig3) / binprob;
      double beta = (-(binprob * sig * sig - 1) + sqrt((binprob * sig * sig - 1) * (binprob * sig * sig - 1) + 4 * dat * sig * sig)) / 2;

      betaPen += (beta - 1) * (beta - 1) / (2 * sig * sig);
      sumLogProb += -dat * log(beta * (binprob));
      sumNorm += beta * (binprob);
    }

    ParameterDict params;
    params["a"] = 1;
    params["b"] = 1;
    params["c"] = 1;
    lh2.SetParameters(params);
    REQUIRE_THAT(lh2.Evaluate(), Catch::Matchers::WithinAbs(sumNorm + sumLogProb + betaPen + constraint, 0.0001));
  }
}

TEST_CASE("Binned NLLH, Grouping Systematics")
{

  SECTION("Two PDFs, two systematics in \"\" Group")
  {

    // Make BinnedEDs
    AxisCollection ax;
    // Two bins, each unit width
    ax.AddAxis(BinAxis("xaxis", 0, 2, 2));
    BinnedED pdf1("pdf1", ax);
    BinnedED pdf2("pdf2", ax);
    BinnedED data("data", ax);

    // Some simple bin contents
    double bin0_pdf1 = 4.;
    double bin1_pdf1 = 6.;
    double bin0_pdf2 = 3.;
    double bin1_pdf2 = 2.;
    double bin0_data = 3.;
    double bin1_data = 4.;
    pdf1.SetBinContent(0, bin0_pdf1);
    pdf1.SetBinContent(1, bin1_pdf1);
    pdf2.SetBinContent(0, bin0_pdf2);
    pdf2.SetBinContent(1, bin1_pdf2);
    data.SetBinContent(0, bin0_data);
    data.SetBinContent(1, bin1_data);

    // Now make 2 simple shift systematics
    Shift *shift1 = new Shift("shift1");
    shift1->RenameParameter("shift", "xshift1");
    shift1->SetShift(0.0);
    ObsSet obs = {"xaxis"};
    shift1->SetAxes(ax);
    shift1->SetTransformationObs(obs);
    shift1->SetDistributionObs(obs);
    shift1->Construct();

    Shift *shift2 = new Shift("shift2");
    shift2->RenameParameter("shift", "xshift2");
    shift2->SetShift(0.0);
    shift2->SetAxes(ax);
    shift2->SetTransformationObs(obs);
    shift2->SetDistributionObs(obs);
    shift2->Construct();

    // Make a LLH object
    BinnedNLLH lh1;
    lh1.SetDataDist(data);
    lh1.AddSystematic(shift1, "");
    lh1.AddSystematic(shift2, "");
    std::vector<std::string> group1vec = {""};
    std::vector<std::string> group2vec = {""};
    lh1.AddPdf(pdf1, group1vec, 100);
    lh1.AddPdf(pdf2, group2vec, 100);
    lh1.RegisterFitComponents();

    // And set some parameter values
    double shiftval1 = 0.5;
    double shiftval2 = 0.1;
    double pdfnorm1 = 1.0;
    double pdfnorm2 = 1.0;
    ParameterDict parameterValues;
    parameterValues["pdf1"] = pdfnorm1;
    parameterValues["pdf2"] = pdfnorm2;
    parameterValues["xshift1"] = shiftval1;
    parameterValues["xshift2"] = shiftval2;

    lh1.SetParameters(parameterValues);
    double llh = lh1.Evaluate();

    // Now let's calculate the llh by hand, applying the systematics ourselves
    // Bin 0 loses events being shifted into Bin 1:
    double bin0_pdf1_shift1 = bin0_pdf1 - shiftval1 * bin0_pdf1;
    // Bin 1 loses events being shifted into the overflow, but gains from Bin 0:
    double bin1_pdf1_shift1 = bin1_pdf1 + (shiftval1 * bin0_pdf1) - (shiftval1 * bin1_pdf1);
    // And we can repeat for the second shift:
    double bin0_pdf1_shift2 = bin0_pdf1_shift1 - shiftval2 * bin0_pdf1_shift1;
    double bin1_pdf1_shift2 = bin1_pdf1_shift1 + (shiftval2 * bin0_pdf1_shift1) - (shiftval2 * bin1_pdf1_shift1);

    // And we have the same for PDF 2:
    double bin0_pdf2_shift1 = bin0_pdf2 - shiftval1 * bin0_pdf2;
    double bin1_pdf2_shift1 = bin1_pdf2 + (shiftval1 * bin0_pdf2) - (shiftval1 * bin1_pdf2);
    double bin0_pdf2_shift2 = bin0_pdf2_shift1 - shiftval2 * bin0_pdf2_shift1;
    double bin1_pdf2_shift2 = bin1_pdf2_shift1 + (shiftval2 * bin0_pdf2_shift1) - (shiftval2 * bin1_pdf2_shift1);

    // Now we can sum the total MC bin contents:
    double bin0_totmc = (pdfnorm1 * bin0_pdf1_shift2) + (pdfnorm2 * bin0_pdf2_shift2);
    double bin1_totmc = (pdfnorm1 * bin1_pdf1_shift2) + (pdfnorm2 * bin1_pdf2_shift2);

    // And calculate the LLH! Summing -data*log(mc) + mc for each bin
    double calcllh = (-bin0_data * log(bin0_totmc) + bin0_totmc) + (-bin1_data * log(bin1_totmc) + bin1_totmc);

    REQUIRE_THAT(llh, Catch::Matchers::WithinAbs(calcllh, 0.0001));
    // And compare to hand calculated value just to be sure:
    REQUIRE_THAT(llh, Catch::Matchers::WithinAbs(-1.03259, 0.0001));
  }

  SECTION("Two PDFs, one systematic in \"\" Group, one in separate group")
  {

    // Make BinnedEDs
    AxisCollection ax;
    // Two bins, each unit width
    ax.AddAxis(BinAxis("xaxis", 0, 2, 2));
    BinnedED pdf1("pdf1", ax);
    BinnedED pdf2("pdf2", ax);
    BinnedED data("data", ax);

    // Some simple bin contents
    double bin0_pdf1 = 4.;
    double bin1_pdf1 = 6.;
    double bin0_pdf2 = 3.;
    double bin1_pdf2 = 2.;
    double bin0_data = 3.;
    double bin1_data = 4.;
    pdf1.SetBinContent(0, bin0_pdf1);
    pdf1.SetBinContent(1, bin1_pdf1);
    pdf2.SetBinContent(0, bin0_pdf2);
    pdf2.SetBinContent(1, bin1_pdf2);
    data.SetBinContent(0, bin0_data);
    data.SetBinContent(1, bin1_data);

    // Now make 2 simple shift systematics
    Shift *shift1 = new Shift("shift1");
    shift1->RenameParameter("shift", "xshift1");
    shift1->SetShift(0.0);
    ObsSet obs = {"xaxis"};
    shift1->SetAxes(ax);
    shift1->SetTransformationObs(obs);
    shift1->SetDistributionObs(obs);
    shift1->Construct();

    Shift *shift2 = new Shift("shift2");
    shift2->RenameParameter("shift", "xshift2");
    shift2->SetShift(0.0);
    shift2->SetAxes(ax);
    shift2->SetTransformationObs(obs);
    shift2->SetDistributionObs(obs);
    shift2->Construct();

    // Make a LLH object
    BinnedNLLH lh1;
    lh1.SetDataDist(data);
    lh1.AddSystematic(shift1, "");
    lh1.AddSystematic(shift2, "syst2");
    std::vector<std::string> group1vec = {""};
    std::vector<std::string> group2vec = {"syst2"};
    lh1.AddPdf(pdf1, group1vec, 100);
    lh1.AddPdf(pdf2, group2vec, 100);
    lh1.RegisterFitComponents();

    // And set some parameter values
    double shiftval1 = 0.5;
    double shiftval2 = 0.1;
    double pdfnorm1 = 1.0;
    double pdfnorm2 = 1.0;
    ParameterDict parameterValues;
    parameterValues["pdf1"] = pdfnorm1;
    parameterValues["pdf2"] = pdfnorm2;
    parameterValues["xshift1"] = shiftval1;
    parameterValues["xshift2"] = shiftval2;

    lh1.SetParameters(parameterValues);
    double llh = lh1.Evaluate();

    // Now let's calculate the llh by hand, applying the systematics ourselves
    // Bin 0 loses events being shifted into Bin 1:
    double bin0_pdf1_shift1 = bin0_pdf1 - shiftval1 * bin0_pdf1;
    // Bin 1 loses events being shifted into the overflow, but gains from Bin 0:
    double bin1_pdf1_shift1 = bin1_pdf1 + (shiftval1 * bin0_pdf1) - (shiftval1 * bin1_pdf1);
    // But the second shift doesn''t apply to PDF 1

    // And we have the same for PDF 2:
    double bin0_pdf2_shift1 = bin0_pdf2 - shiftval1 * bin0_pdf2;
    double bin1_pdf2_shift1 = bin1_pdf2 + (shiftval1 * bin0_pdf2) - (shiftval1 * bin1_pdf2);
    double bin0_pdf2_shift2 = bin0_pdf2_shift1 - shiftval2 * bin0_pdf2_shift1;
    double bin1_pdf2_shift2 = bin1_pdf2_shift1 + (shiftval2 * bin0_pdf2_shift1) - (shiftval2 * bin1_pdf2_shift1);

    // Now we can sum the total MC bin contents:
    double bin0_totmc = (pdfnorm1 * bin0_pdf1_shift1) + (pdfnorm2 * bin0_pdf2_shift2);
    double bin1_totmc = (pdfnorm1 * bin1_pdf1_shift1) + (pdfnorm2 * bin1_pdf2_shift2);

    // And calculate the LLH! Summing -data*log(mc) + mc for each bin
    double calcllh = (-bin0_data * log(bin0_totmc) + bin0_totmc) + (-bin1_data * log(bin1_totmc) + bin1_totmc);

    REQUIRE_THAT(llh, Catch::Matchers::WithinAbs(calcllh, 0.0001));
    REQUIRE_THAT(llh, Catch::Matchers::WithinAbs(-0.88280, 0.0001));
  }

  SECTION("Two PDFs, two systematics in separate groups")
  {

    // Make BinnedEDs
    AxisCollection ax;
    // Two bins, each unit width
    ax.AddAxis(BinAxis("xaxis", 0, 2, 2));
    BinnedED pdf1("pdf1", ax);
    BinnedED pdf2("pdf2", ax);
    BinnedED data("data", ax);

    // Some simple bin contents
    double bin0_pdf1 = 4.;
    double bin1_pdf1 = 6.;
    double bin0_pdf2 = 3.;
    double bin1_pdf2 = 2.;
    double bin0_data = 3.;
    double bin1_data = 4.;
    pdf1.SetBinContent(0, bin0_pdf1);
    pdf1.SetBinContent(1, bin1_pdf1);
    pdf2.SetBinContent(0, bin0_pdf2);
    pdf2.SetBinContent(1, bin1_pdf2);
    data.SetBinContent(0, bin0_data);
    data.SetBinContent(1, bin1_data);

    // Now make 2 simple shift systematics
    Shift *shift1 = new Shift("shift1");
    shift1->RenameParameter("shift", "xshift1");
    shift1->SetShift(0.0);
    ObsSet obs = {"xaxis"};
    shift1->SetAxes(ax);
    shift1->SetTransformationObs(obs);
    shift1->SetDistributionObs(obs);
    shift1->Construct();

    Shift *shift2 = new Shift("shift2");
    shift2->RenameParameter("shift", "xshift2");
    shift2->SetShift(0.0);
    shift2->SetAxes(ax);
    shift2->SetTransformationObs(obs);
    shift2->SetDistributionObs(obs);
    shift2->Construct();

    // Make a LLH object
    BinnedNLLH lh1;
    lh1.SetDataDist(data);
    lh1.AddSystematic(shift1, "syst1");
    lh1.AddSystematic(shift2, "syst2");
    std::vector<std::string> group1vec = {"syst1"};
    std::vector<std::string> group2vec = {"syst2"};
    lh1.AddPdf(pdf1, group1vec, 100);
    lh1.AddPdf(pdf2, group2vec, 100);
    lh1.RegisterFitComponents();

    // And set some parameter values
    double shiftval1 = 0.5;
    double shiftval2 = 0.1;
    double pdfnorm1 = 1.0;
    double pdfnorm2 = 1.0;
    ParameterDict parameterValues;
    parameterValues["pdf1"] = pdfnorm1;
    parameterValues["pdf2"] = pdfnorm2;
    parameterValues["xshift1"] = shiftval1;
    parameterValues["xshift2"] = shiftval2;

    lh1.SetParameters(parameterValues);
    double llh = lh1.Evaluate();

    // Now let's calculate the llh by hand, applying the systematics ourselves
    // Bin 0 loses events being shifted into Bin 1:
    double bin0_pdf1_shift1 = bin0_pdf1 - shiftval1 * bin0_pdf1;
    // Bin 1 loses events being shifted into the overflow, but gains from Bin 0:
    double bin1_pdf1_shift1 = bin1_pdf1 + (shiftval1 * bin0_pdf1) - (shiftval1 * bin1_pdf1);
    // But the second shift doesn''t apply to PDF 1

    // And we have the same for PDF 2 and shift2:
    double bin0_pdf2_shift2 = bin0_pdf2 - shiftval2 * bin0_pdf2;
    double bin1_pdf2_shift2 = bin1_pdf2 + (shiftval2 * bin0_pdf2) - (shiftval2 * bin1_pdf2);

    // Now we can sum the total MC bin contents:
    double bin0_totmc = (pdfnorm1 * bin0_pdf1_shift1) + (pdfnorm2 * bin0_pdf2_shift2);
    double bin1_totmc = (pdfnorm1 * bin1_pdf1_shift1) + (pdfnorm2 * bin1_pdf2_shift2);

    // And calculate the LLH! Summing -data*log(mc) + mc for each bin
    double calcllh = (-bin0_data * log(bin0_totmc) + bin0_totmc) + (-bin1_data * log(bin1_totmc) + bin1_totmc);

    REQUIRE_THAT(llh, Catch::Matchers::WithinAbs(calcllh, 0.0001));
    REQUIRE_THAT(llh, Catch::Matchers::WithinAbs(-0.68307, 0.0001));
  }

  SECTION("Two PDFs, one systematic in \"\" Group, one in separate group which contains both PDFs")
  {

    // Make BinnedEDs
    AxisCollection ax;
    // Two bins, each unit width
    ax.AddAxis(BinAxis("xaxis", 0, 2, 2));
    BinnedED pdf1("pdf1", ax);
    BinnedED pdf2("pdf2", ax);
    BinnedED data("data", ax);

    // Some simple bin contents
    double bin0_pdf1 = 4.;
    double bin1_pdf1 = 6.;
    double bin0_pdf2 = 3.;
    double bin1_pdf2 = 2.;
    double bin0_data = 3.;
    double bin1_data = 4.;
    pdf1.SetBinContent(0, bin0_pdf1);
    pdf1.SetBinContent(1, bin1_pdf1);
    pdf2.SetBinContent(0, bin0_pdf2);
    pdf2.SetBinContent(1, bin1_pdf2);
    data.SetBinContent(0, bin0_data);
    data.SetBinContent(1, bin1_data);

    // Now make 2 simple shift systematics
    Shift *shift1 = new Shift("shift1");
    shift1->RenameParameter("shift", "xshift1");
    shift1->SetShift(0.0);
    ObsSet obs = {"xaxis"};
    shift1->SetAxes(ax);
    shift1->SetTransformationObs(obs);
    shift1->SetDistributionObs(obs);
    shift1->Construct();

    Shift *shift2 = new Shift("shift2");
    shift2->RenameParameter("shift", "xshift2");
    shift2->SetShift(0.0);
    shift2->SetAxes(ax);
    shift2->SetTransformationObs(obs);
    shift2->SetDistributionObs(obs);
    shift2->Construct();

    // Make a LLH object
    BinnedNLLH lh1;
    lh1.SetDataDist(data);
    lh1.AddSystematic(shift1, "");
    lh1.AddSystematic(shift2, "syst2");
    std::vector<std::string> group1vec = {"syst2"};
    std::vector<std::string> group2vec = {"syst2"};
    lh1.AddPdf(pdf1, group1vec, 100);
    lh1.AddPdf(pdf2, group2vec, 100);
    lh1.RegisterFitComponents();

    // And set some parameter values
    double shiftval1 = 0.5;
    double shiftval2 = 0.1;
    double pdfnorm1 = 1.0;
    double pdfnorm2 = 1.0;
    ParameterDict parameterValues;
    parameterValues["pdf1"] = pdfnorm1;
    parameterValues["pdf2"] = pdfnorm2;
    parameterValues["xshift1"] = shiftval1;
    parameterValues["xshift2"] = shiftval2;

    lh1.SetParameters(parameterValues);
    double llh = lh1.Evaluate();

    // Now let's calculate the llh by hand, applying the systematics ourselves
    // Bin 0 loses events being shifted into Bin 1:
    double bin0_pdf1_shift1 = bin0_pdf1 - shiftval1 * bin0_pdf1;
    // Bin 1 loses events being shifted into the overflow, but gains from Bin 0:
    double bin1_pdf1_shift1 = bin1_pdf1 + (shiftval1 * bin0_pdf1) - (shiftval1 * bin1_pdf1);
    // And we can repeat for the second shift:
    double bin0_pdf1_shift2 = bin0_pdf1_shift1 - shiftval2 * bin0_pdf1_shift1;
    double bin1_pdf1_shift2 = bin1_pdf1_shift1 + (shiftval2 * bin0_pdf1_shift1) - (shiftval2 * bin1_pdf1_shift1);

    // And we have the same for PDF 2 and shift2:
    double bin0_pdf2_shift1 = bin0_pdf2 - shiftval1 * bin0_pdf2;
    double bin1_pdf2_shift1 = bin1_pdf2 + (shiftval1 * bin0_pdf2) - (shiftval1 * bin1_pdf2);
    double bin0_pdf2_shift2 = bin0_pdf2_shift1 - shiftval2 * bin0_pdf2_shift1;
    double bin1_pdf2_shift2 = bin1_pdf2_shift1 + (shiftval2 * bin0_pdf2_shift1) - (shiftval2 * bin1_pdf2_shift1);

    // Now we can sum the total MC bin contents:
    double bin0_totmc = (pdfnorm1 * bin0_pdf1_shift2) + (pdfnorm2 * bin0_pdf2_shift2);
    double bin1_totmc = (pdfnorm1 * bin1_pdf1_shift2) + (pdfnorm2 * bin1_pdf2_shift2);

    // And calculate the LLH! Summing -data*log(mc) + mc for each bin
    double calcllh = (-bin0_data * log(bin0_totmc) + bin0_totmc) + (-bin1_data * log(bin1_totmc) + bin1_totmc);

    REQUIRE_THAT(llh, Catch::Matchers::WithinAbs(calcllh, 0.0001));
    REQUIRE_THAT(llh, Catch::Matchers::WithinAbs(-1.03259, 0.0001));
  }

  SECTION("Three PDFs, two systematics with separate groups")
  {

    // Make BinnedEDs
    AxisCollection ax;
    // Two bins, each unit width
    ax.AddAxis(BinAxis("xaxis", 0, 2, 2));
    BinnedED pdf1("pdf1", ax);
    BinnedED pdf2("pdf2", ax);
    BinnedED pdf3("pdf3", ax);
    BinnedED data("data", ax);

    // Some simple bin contents
    double bin0_pdf1 = 4.;
    double bin1_pdf1 = 6.;
    double bin0_pdf2 = 3.;
    double bin1_pdf2 = 2.;
    double bin0_pdf3 = 2.;
    double bin1_pdf3 = 4.;
    double bin0_data = 3.;
    double bin1_data = 4.;
    pdf1.SetBinContent(0, bin0_pdf1);
    pdf1.SetBinContent(1, bin1_pdf1);
    pdf2.SetBinContent(0, bin0_pdf2);
    pdf2.SetBinContent(1, bin1_pdf2);
    pdf3.SetBinContent(0, bin0_pdf3);
    pdf3.SetBinContent(1, bin1_pdf3);
    data.SetBinContent(0, bin0_data);
    data.SetBinContent(1, bin1_data);

    // Now make 2 simple shift systematics
    Shift *shift1 = new Shift("shift1");
    shift1->RenameParameter("shift", "xshift1");
    shift1->SetShift(0.0);
    ObsSet obs = {"xaxis"};
    shift1->SetAxes(ax);
    shift1->SetTransformationObs(obs);
    shift1->SetDistributionObs(obs);
    shift1->Construct();

    Shift *shift2 = new Shift("shift2");
    shift2->RenameParameter("shift", "xshift2");
    shift2->SetShift(0.0);
    shift2->SetAxes(ax);
    shift2->SetTransformationObs(obs);
    shift2->SetDistributionObs(obs);
    shift2->Construct();

    // Make a LLH object
    BinnedNLLH lh1;
    lh1.SetDataDist(data);
    lh1.AddSystematic(shift1, "syst1");
    lh1.AddSystematic(shift2, "syst2");
    std::vector<std::string> group1vec = {"syst1"};
    std::vector<std::string> group2vec = {"syst2"};
    std::vector<std::string> group3vec = {"syst2"};
    lh1.AddPdf(pdf1, group1vec, 100);
    lh1.AddPdf(pdf2, group2vec, 100);
    lh1.AddPdf(pdf3, group3vec, 100);
    lh1.RegisterFitComponents();

    // And set some parameter values
    double shiftval1 = 0.5;
    double shiftval2 = 0.1;
    double pdfnorm1 = 1.0;
    double pdfnorm2 = 1.0;
    double pdfnorm3 = 1.0;
    ParameterDict parameterValues;
    parameterValues["pdf1"] = pdfnorm1;
    parameterValues["pdf2"] = pdfnorm2;
    parameterValues["pdf3"] = pdfnorm3;
    parameterValues["xshift1"] = shiftval1;
    parameterValues["xshift2"] = shiftval2;

    lh1.SetParameters(parameterValues);
    double llh = lh1.Evaluate();

    // Now let's calculate the llh by hand, applying the systematics ourselves
    // Bin 0 loses events being shifted into Bin 1:
    double bin0_pdf1_shift1 = bin0_pdf1 - shiftval1 * bin0_pdf1;
    // Bin 1 loses events being shifted into the overflow, but gains from Bin 0:
    double bin1_pdf1_shift1 = bin1_pdf1 + (shiftval1 * bin0_pdf1) - (shiftval1 * bin1_pdf1);
    // But the second shift doesn't apply to PDF 1

    // And we have shift2 for PDFs 2 and 3:
    double bin0_pdf2_shift2 = bin0_pdf2 - shiftval2 * bin0_pdf2;
    double bin1_pdf2_shift2 = bin1_pdf2 + (shiftval2 * bin0_pdf2) - (shiftval2 * bin1_pdf2);
    double bin0_pdf3_shift2 = bin0_pdf3 - shiftval2 * bin0_pdf3;
    double bin1_pdf3_shift2 = bin1_pdf3 + (shiftval2 * bin0_pdf3) - (shiftval2 * bin1_pdf3);

    // Now we can sum the total MC bin contents:
    double bin0_totmc = (pdfnorm1 * bin0_pdf1_shift1) + (pdfnorm2 * bin0_pdf2_shift2) + (pdfnorm3 * bin0_pdf3_shift2);
    double bin1_totmc = (pdfnorm1 * bin1_pdf1_shift1) + (pdfnorm2 * bin1_pdf2_shift2) + (pdfnorm3 * bin1_pdf3_shift2);

    // And calculate the LLH! Summing -data*log(mc) + mc for each bin
    double calcllh = (-bin0_data * log(bin0_totmc) + bin0_totmc) + (-bin1_data * log(bin1_totmc) + bin1_totmc);

    REQUIRE_THAT(llh, Catch::Matchers::WithinAbs(calcllh, 0.0001));
    REQUIRE_THAT(llh, Catch::Matchers::WithinAbs(2.22954, 0.0001));
  }

  SECTION("Three PDFs, two systematics in different groups, one in \"\" group")
  {

    // Make BinnedEDs
    AxisCollection ax;
    // Two bins, each unit width
    ax.AddAxis(BinAxis("xaxis", 0, 2, 2));
    BinnedED pdf1("pdf1", ax);
    BinnedED pdf2("pdf2", ax);
    BinnedED pdf3("pdf3", ax);
    BinnedED data("data", ax);

    // Some simple bin contents
    double bin0_pdf1 = 4.;
    double bin1_pdf1 = 6.;
    double bin0_pdf2 = 3.;
    double bin1_pdf2 = 2.;
    double bin0_pdf3 = 2.;
    double bin1_pdf3 = 4.;
    double bin0_data = 3.;
    double bin1_data = 4.;
    pdf1.SetBinContent(0, bin0_pdf1);
    pdf1.SetBinContent(1, bin1_pdf1);
    pdf2.SetBinContent(0, bin0_pdf2);
    pdf2.SetBinContent(1, bin1_pdf2);
    pdf3.SetBinContent(0, bin0_pdf3);
    pdf3.SetBinContent(1, bin1_pdf3);
    data.SetBinContent(0, bin0_data);
    data.SetBinContent(1, bin1_data);

    // Now make 2 simple shift systematics
    Shift *shift1 = new Shift("shift1");
    shift1->RenameParameter("shift", "xshift1");
    shift1->SetShift(0.0);
    ObsSet obs = {"xaxis"};
    shift1->SetAxes(ax);
    shift1->SetTransformationObs(obs);
    shift1->SetDistributionObs(obs);
    shift1->Construct();

    Shift *shift2 = new Shift("shift2");
    shift2->RenameParameter("shift", "xshift2");
    shift2->SetShift(0.0);
    shift2->SetAxes(ax);
    shift2->SetTransformationObs(obs);
    shift2->SetDistributionObs(obs);
    shift2->Construct();

    Shift *shift3 = new Shift("shift3");
    shift3->RenameParameter("shift", "xshift3");
    shift3->SetShift(0.0);
    shift3->SetAxes(ax);
    shift3->SetTransformationObs(obs);
    shift3->SetDistributionObs(obs);
    shift3->Construct();

    // Make a LLH object
    BinnedNLLH lh1;
    lh1.SetDataDist(data);
    lh1.AddSystematic(shift1, "");
    lh1.AddSystematic(shift2, "syst2");
    lh1.AddSystematic(shift3, "syst3");
    std::vector<std::string> group1vec = {""};
    std::vector<std::string> group2vec = {"syst2"};
    std::vector<std::string> group3vec = {"syst3"};
    lh1.AddPdf(pdf1, group1vec, 100);
    lh1.AddPdf(pdf2, group2vec, 100);
    lh1.AddPdf(pdf3, group3vec, 100);
    lh1.RegisterFitComponents();

    // And set some parameter values
    double shiftval1 = 0.5;
    double shiftval2 = 0.1;
    double shiftval3 = -0.3;
    double pdfnorm1 = 1.0;
    double pdfnorm2 = 1.0;
    double pdfnorm3 = 1.0;
    ParameterDict parameterValues;
    parameterValues["pdf1"] = pdfnorm1;
    parameterValues["pdf2"] = pdfnorm2;
    parameterValues["pdf3"] = pdfnorm3;
    parameterValues["xshift1"] = shiftval1;
    parameterValues["xshift2"] = shiftval2;
    parameterValues["xshift3"] = shiftval3;

    lh1.SetParameters(parameterValues);
    double llh = lh1.Evaluate();

    // Now let's calculate the llh by hand, applying the systematics ourselves
    // Bin 0 loses events being shifted into Bin 1:
    double bin0_pdf1_shift1 = bin0_pdf1 - shiftval1 * bin0_pdf1;
    // Bin 1 loses events being shifted into the overflow, but gains from Bin 0:
    double bin1_pdf1_shift1 = bin1_pdf1 + (shiftval1 * bin0_pdf1) - (shiftval1 * bin1_pdf1);
    // But the second and third shifts don't apply to PDF 1

    // And we have shift2 for PDF 2:
    double bin0_pdf2_shift1 = bin0_pdf2 - shiftval1 * bin0_pdf2;
    double bin1_pdf2_shift1 = bin1_pdf2 + (shiftval1 * bin0_pdf2) - (shiftval1 * bin1_pdf2);
    double bin0_pdf2_shift2 = bin0_pdf2_shift1 - shiftval2 * bin0_pdf2_shift1;
    double bin1_pdf2_shift2 = bin1_pdf2_shift1 + (shiftval2 * bin0_pdf2_shift1) - (shiftval2 * bin1_pdf2_shift1);

    // And we have shift3 for PDFs 3:
    double bin0_pdf3_shift1 = bin0_pdf3 - shiftval1 * bin0_pdf3;
    double bin1_pdf3_shift1 = bin1_pdf3 + (shiftval1 * bin0_pdf3) - (shiftval1 * bin1_pdf3);
    // As shift 3 has a negative value, the shift goes the other way
    double bin0_pdf3_shift3 = bin0_pdf3_shift1 - (shiftval3 * bin1_pdf3_shift1) + (shiftval3 * bin0_pdf3_shift1);
    double bin1_pdf3_shift3 = bin1_pdf3_shift1 + shiftval3 * bin1_pdf3_shift1;

    // Now we can sum the total MC bin contents:
    double bin0_totmc = (pdfnorm1 * bin0_pdf1_shift1) + (pdfnorm2 * bin0_pdf2_shift2) + (pdfnorm3 * bin0_pdf3_shift3);
    double bin1_totmc = (pdfnorm1 * bin1_pdf1_shift1) + (pdfnorm2 * bin1_pdf2_shift2) + (pdfnorm3 * bin1_pdf3_shift3);

    // And calculate the LLH! Summing -data*log(mc) + mc for each bin
    double calcllh = (-bin0_data * log(bin0_totmc) + bin0_totmc) + (-bin1_data * log(bin1_totmc) + bin1_totmc);

    REQUIRE_THAT(llh, Catch::Matchers::WithinAbs(calcllh, 0.0001));
    REQUIRE_THAT(llh, Catch::Matchers::WithinAbs(0.64667, 0.0001));
  }

  SECTION("Three PDFs, two systematics in \"\" group, one in a separate group")
  {

    // Make BinnedEDs
    AxisCollection ax;
    // Two bins, each unit width
    ax.AddAxis(BinAxis("xaxis", 0, 2, 2));
    BinnedED pdf1("pdf1", ax);
    BinnedED pdf2("pdf2", ax);
    BinnedED pdf3("pdf3", ax);
    BinnedED data("data", ax);

    // Some simple bin contents
    double bin0_pdf1 = 4.;
    double bin1_pdf1 = 6.;
    double bin0_pdf2 = 3.;
    double bin1_pdf2 = 2.;
    double bin0_pdf3 = 2.;
    double bin1_pdf3 = 4.;
    double bin0_data = 3.;
    double bin1_data = 4.;
    pdf1.SetBinContent(0, bin0_pdf1);
    pdf1.SetBinContent(1, bin1_pdf1);
    pdf2.SetBinContent(0, bin0_pdf2);
    pdf2.SetBinContent(1, bin1_pdf2);
    pdf3.SetBinContent(0, bin0_pdf3);
    pdf3.SetBinContent(1, bin1_pdf3);
    data.SetBinContent(0, bin0_data);
    data.SetBinContent(1, bin1_data);

    // Now make 2 simple shift systematics
    Shift *shift1 = new Shift("shift1");
    shift1->RenameParameter("shift", "xshift1");
    shift1->SetShift(0.0);
    ObsSet obs = {"xaxis"};
    shift1->SetAxes(ax);
    shift1->SetTransformationObs(obs);
    shift1->SetDistributionObs(obs);
    shift1->Construct();

    Shift *shift2 = new Shift("shift2");
    shift2->RenameParameter("shift", "xshift2");
    shift2->SetShift(0.0);
    shift2->SetAxes(ax);
    shift2->SetTransformationObs(obs);
    shift2->SetDistributionObs(obs);
    shift2->Construct();

    Shift *shift3 = new Shift("shift3");
    shift3->RenameParameter("shift", "xshift3");
    shift3->SetShift(0.0);
    shift3->SetAxes(ax);
    shift3->SetTransformationObs(obs);
    shift3->SetDistributionObs(obs);
    shift3->Construct();

    // Make a LLH object
    BinnedNLLH lh1;
    lh1.SetDataDist(data);
    lh1.AddSystematic(shift1, "");
    lh1.AddSystematic(shift2, "");
    lh1.AddSystematic(shift3, "syst3");
    std::vector<std::string> group1vec = {""};
    std::vector<std::string> group2vec = {""};
    std::vector<std::string> group3vec = {"syst3"};
    lh1.AddPdf(pdf1, group1vec, 100);
    lh1.AddPdf(pdf2, group2vec, 100);
    lh1.AddPdf(pdf3, group3vec, 100);
    lh1.RegisterFitComponents();

    // And set some parameter values
    double shiftval1 = 0.5;
    double shiftval2 = 0.1;
    double shiftval3 = -0.3;
    double pdfnorm1 = 1.0;
    double pdfnorm2 = 1.0;
    double pdfnorm3 = 1.0;
    ParameterDict parameterValues;
    parameterValues["pdf1"] = pdfnorm1;
    parameterValues["pdf2"] = pdfnorm2;
    parameterValues["pdf3"] = pdfnorm3;
    parameterValues["xshift1"] = shiftval1;
    parameterValues["xshift2"] = shiftval2;
    parameterValues["xshift3"] = shiftval3;

    lh1.SetParameters(parameterValues);
    double llh = lh1.Evaluate();

    // Now let's calculate the llh by hand, applying the systematics ourselves
    // Bin 0 loses events being shifted into Bin 1:
    double bin0_pdf1_shift1 = bin0_pdf1 - shiftval1 * bin0_pdf1;
    // Bin 1 loses events being shifted into the overflow, but gains from Bin 0:
    double bin1_pdf1_shift1 = bin1_pdf1 + (shiftval1 * bin0_pdf1) - (shiftval1 * bin1_pdf1);
    // And do the same for shift 2:
    double bin0_pdf1_shift2 = bin0_pdf1_shift1 - shiftval2 * bin0_pdf1_shift1;
    double bin1_pdf1_shift2 = bin1_pdf1_shift1 + (shiftval2 * bin0_pdf1_shift1) - (shiftval2 * bin1_pdf1_shift1);

    // And similarly for PDF 2:
    double bin0_pdf2_shift1 = bin0_pdf2 - shiftval1 * bin0_pdf2;
    double bin1_pdf2_shift1 = bin1_pdf2 + (shiftval1 * bin0_pdf2) - (shiftval1 * bin1_pdf2);
    double bin0_pdf2_shift2 = bin0_pdf2_shift1 - shiftval2 * bin0_pdf2_shift1;
    double bin1_pdf2_shift2 = bin1_pdf2_shift1 + (shiftval2 * bin0_pdf2_shift1) - (shiftval2 * bin1_pdf2_shift1);

    // And we have shift3 for PDFs 3:
    double bin0_pdf3_shift1 = bin0_pdf3 - shiftval1 * bin0_pdf3;
    double bin1_pdf3_shift1 = bin1_pdf3 + (shiftval1 * bin0_pdf3) - (shiftval1 * bin1_pdf3);
    double bin0_pdf3_shift2 = bin0_pdf3_shift1 - shiftval2 * bin0_pdf3_shift1;
    double bin1_pdf3_shift2 = bin1_pdf3_shift1 + (shiftval2 * bin0_pdf3_shift1) - (shiftval2 * bin1_pdf3_shift1);
    // As shift 3 has a negative value, the shift goes the other way
    double bin0_pdf3_shift3 = bin0_pdf3_shift2 - (shiftval3 * bin1_pdf3_shift2) + (shiftval3 * bin0_pdf3_shift2);
    double bin1_pdf3_shift3 = bin1_pdf3_shift2 + shiftval3 * bin1_pdf3_shift2;

    // Now we can sum the total MC bin contents:
    double bin0_totmc = (pdfnorm1 * bin0_pdf1_shift2) + (pdfnorm2 * bin0_pdf2_shift2) + (pdfnorm3 * bin0_pdf3_shift3);
    double bin1_totmc = (pdfnorm1 * bin1_pdf1_shift2) + (pdfnorm2 * bin1_pdf2_shift2) + (pdfnorm3 * bin1_pdf3_shift3);

    // And calculate the LLH! Summing -data*log(mc) + mc for each bin
    double calcllh = (-bin0_data * log(bin0_totmc) + bin0_totmc) + (-bin1_data * log(bin1_totmc) + bin1_totmc);

    REQUIRE_THAT(llh, Catch::Matchers::WithinAbs(calcllh, 0.0001));
    REQUIRE_THAT(llh, Catch::Matchers::WithinAbs(0.27339, 0.0001));
  }

  SECTION("Three PDFs, three systematics all in separate groups")
  {

    // Make BinnedEDs
    AxisCollection ax;
    // Two bins, each unit width
    ax.AddAxis(BinAxis("xaxis", 0, 2, 2));
    BinnedED pdf1("pdf1", ax);
    BinnedED pdf2("pdf2", ax);
    BinnedED pdf3("pdf3", ax);
    BinnedED data("data", ax);

    // Some simple bin contents
    double bin0_pdf1 = 4.;
    double bin1_pdf1 = 6.;
    double bin0_pdf2 = 3.;
    double bin1_pdf2 = 2.;
    double bin0_pdf3 = 2.;
    double bin1_pdf3 = 4.;
    double bin0_data = 3.;
    double bin1_data = 4.;
    pdf1.SetBinContent(0, bin0_pdf1);
    pdf1.SetBinContent(1, bin1_pdf1);
    pdf2.SetBinContent(0, bin0_pdf2);
    pdf2.SetBinContent(1, bin1_pdf2);
    pdf3.SetBinContent(0, bin0_pdf3);
    pdf3.SetBinContent(1, bin1_pdf3);
    data.SetBinContent(0, bin0_data);
    data.SetBinContent(1, bin1_data);

    // Now make 2 simple shift systematics
    Shift *shift1 = new Shift("shift1");
    shift1->RenameParameter("shift", "xshift1");
    shift1->SetShift(0.0);
    ObsSet obs = {"xaxis"};
    shift1->SetAxes(ax);
    shift1->SetTransformationObs(obs);
    shift1->SetDistributionObs(obs);
    shift1->Construct();

    Shift *shift2 = new Shift("shift2");
    shift2->RenameParameter("shift", "xshift2");
    shift2->SetShift(0.0);
    shift2->SetAxes(ax);
    shift2->SetTransformationObs(obs);
    shift2->SetDistributionObs(obs);
    shift2->Construct();

    Shift *shift3 = new Shift("shift3");
    shift3->RenameParameter("shift", "xshift3");
    shift3->SetShift(0.0);
    shift3->SetAxes(ax);
    shift3->SetTransformationObs(obs);
    shift3->SetDistributionObs(obs);
    shift3->Construct();

    // Make a LLH object
    BinnedNLLH lh1;
    lh1.SetDataDist(data);
    lh1.AddSystematic(shift1, "syst1");
    lh1.AddSystematic(shift2, "syst2");
    lh1.AddSystematic(shift3, "syst3");
    std::vector<std::string> group1vec = {"syst1"};
    std::vector<std::string> group2vec = {"syst2"};
    std::vector<std::string> group3vec = {"syst3"};
    lh1.AddPdf(pdf1, group1vec, 100);
    lh1.AddPdf(pdf2, group2vec, 100);
    lh1.AddPdf(pdf3, group3vec, 100);
    lh1.RegisterFitComponents();

    // And set some parameter values
    double shiftval1 = 0.5;
    double shiftval2 = 0.1;
    double shiftval3 = -0.3;
    double pdfnorm1 = 1.0;
    double pdfnorm2 = 1.0;
    double pdfnorm3 = 1.0;
    ParameterDict parameterValues;
    parameterValues["pdf1"] = pdfnorm1;
    parameterValues["pdf2"] = pdfnorm2;
    parameterValues["pdf3"] = pdfnorm3;
    parameterValues["xshift1"] = shiftval1;
    parameterValues["xshift2"] = shiftval2;
    parameterValues["xshift3"] = shiftval3;

    lh1.SetParameters(parameterValues);
    double llh = lh1.Evaluate();

    // Now let's calculate the llh by hand, applying the systematics ourselves
    // Bin 0 loses events being shifted into Bin 1:
    double bin0_pdf1_shift1 = bin0_pdf1 - shiftval1 * bin0_pdf1;
    // Bin 1 loses events being shifted into the overflow, but gains from Bin 0:
    double bin1_pdf1_shift1 = bin1_pdf1 + (shiftval1 * bin0_pdf1) - (shiftval1 * bin1_pdf1);

    // And similarly for PDF 2:
    double bin0_pdf2_shift2 = bin0_pdf2 - shiftval2 * bin0_pdf2;
    double bin1_pdf2_shift2 = bin1_pdf2 + (shiftval2 * bin0_pdf2) - (shiftval2 * bin1_pdf2);

    // And we have shift3 for PDFs 3:
    double bin0_pdf3_shift3 = bin0_pdf3 - (shiftval3 * bin1_pdf3) + (shiftval3 * bin0_pdf3);
    double bin1_pdf3_shift3 = bin1_pdf3 + shiftval3 * bin1_pdf3;

    // Now we can sum the total MC bin contents:
    double bin0_totmc = (pdfnorm1 * bin0_pdf1_shift1) + (pdfnorm2 * bin0_pdf2_shift2) + (pdfnorm3 * bin0_pdf3_shift3);
    double bin1_totmc = (pdfnorm1 * bin1_pdf1_shift1) + (pdfnorm2 * bin1_pdf2_shift2) + (pdfnorm3 * bin1_pdf3_shift3);

    // And calculate the LLH! Summing -data*log(mc) + mc for each bin
    double calcllh = (-bin0_data * log(bin0_totmc) + bin0_totmc) + (-bin1_data * log(bin1_totmc) + bin1_totmc);

    REQUIRE_THAT(llh, Catch::Matchers::WithinAbs(calcllh, 0.0001));
    REQUIRE_THAT(llh, Catch::Matchers::WithinAbs(2.06624, 0.0001));
  }

  SECTION("Three PDFs, one systematic in the \"\" group, one applies to one PDF, one applies to two PDFs")
  {

    // Make BinnedEDs
    AxisCollection ax;
    // Two bins, each unit width
    ax.AddAxis(BinAxis("xaxis", 0, 2, 2));
    BinnedED pdf1("pdf1", ax);
    BinnedED pdf2("pdf2", ax);
    BinnedED pdf3("pdf3", ax);
    BinnedED data("data", ax);

    // Some simple bin contents
    double bin0_pdf1 = 4.;
    double bin1_pdf1 = 6.;
    double bin0_pdf2 = 3.;
    double bin1_pdf2 = 2.;
    double bin0_pdf3 = 2.;
    double bin1_pdf3 = 4.;
    double bin0_data = 3.;
    double bin1_data = 4.;
    pdf1.SetBinContent(0, bin0_pdf1);
    pdf1.SetBinContent(1, bin1_pdf1);
    pdf2.SetBinContent(0, bin0_pdf2);
    pdf2.SetBinContent(1, bin1_pdf2);
    pdf3.SetBinContent(0, bin0_pdf3);
    pdf3.SetBinContent(1, bin1_pdf3);
    data.SetBinContent(0, bin0_data);
    data.SetBinContent(1, bin1_data);

    // Now make 2 simple shift systematics
    Shift *shift1 = new Shift("shift1");
    shift1->RenameParameter("shift", "xshift1");
    shift1->SetShift(0.0);
    ObsSet obs = {"xaxis"};
    shift1->SetAxes(ax);
    shift1->SetTransformationObs(obs);
    shift1->SetDistributionObs(obs);
    shift1->Construct();

    Shift *shift2 = new Shift("shift2");
    shift2->RenameParameter("shift", "xshift2");
    shift2->SetShift(0.0);
    shift2->SetAxes(ax);
    shift2->SetTransformationObs(obs);
    shift2->SetDistributionObs(obs);
    shift2->Construct();

    Shift *shift3 = new Shift("shift3");
    shift3->RenameParameter("shift", "xshift3");
    shift3->SetShift(0.0);
    shift3->SetAxes(ax);
    shift3->SetTransformationObs(obs);
    shift3->SetDistributionObs(obs);
    shift3->Construct();

    // Make a LLH object
    BinnedNLLH lh1;
    lh1.SetDataDist(data);
    lh1.AddSystematic(shift1, "");
    lh1.AddSystematic(shift2, "syst2");
    lh1.AddSystematic(shift3, "syst3");
    std::vector<std::string> group1vec = {""};
    std::vector<std::string> group2vec = {"syst2", "syst3"};
    std::vector<std::string> group3vec = {"syst3"};
    lh1.AddPdf(pdf1, group1vec, 100);
    lh1.AddPdf(pdf2, group2vec, 100);
    lh1.AddPdf(pdf3, group3vec, 100);
    lh1.RegisterFitComponents();

    // And set some parameter values
    double shiftval1 = 0.5;
    double shiftval2 = 0.1;
    double shiftval3 = -0.3;
    double pdfnorm1 = 1.0;
    double pdfnorm2 = 1.0;
    double pdfnorm3 = 1.0;
    ParameterDict parameterValues;
    parameterValues["pdf1"] = pdfnorm1;
    parameterValues["pdf2"] = pdfnorm2;
    parameterValues["pdf3"] = pdfnorm3;
    parameterValues["xshift1"] = shiftval1;
    parameterValues["xshift2"] = shiftval2;
    parameterValues["xshift3"] = shiftval3;

    lh1.SetParameters(parameterValues);
    double llh = lh1.Evaluate();

    // Now let's calculate the llh by hand, applying the systematics ourselves
    // Bin 0 loses events being shifted into Bin 1:
    double bin0_pdf1_shift1 = bin0_pdf1 - shiftval1 * bin0_pdf1;
    // Bin 1 loses events being shifted into the overflow, but gains from Bin 0:
    double bin1_pdf1_shift1 = bin1_pdf1 + (shiftval1 * bin0_pdf1) - (shiftval1 * bin1_pdf1);

    // And shifts 1, 2 and 3 for PDF 2:
    double bin0_pdf2_shift1 = bin0_pdf2 - shiftval1 * bin0_pdf2;
    double bin1_pdf2_shift1 = bin1_pdf2 + (shiftval1 * bin0_pdf2) - (shiftval1 * bin1_pdf2);
    double bin0_pdf2_shift2 = bin0_pdf2_shift1 - shiftval2 * bin0_pdf2_shift1;
    double bin1_pdf2_shift2 = bin1_pdf2_shift1 + (shiftval2 * bin0_pdf2_shift1) - (shiftval2 * bin1_pdf2_shift1);
    // As shift 3 has a negative value, the shift goes the other way
    double bin0_pdf2_shift3 = bin0_pdf2_shift2 - (shiftval3 * bin1_pdf2_shift2) + (shiftval3 * bin0_pdf2_shift2);
    double bin1_pdf2_shift3 = bin1_pdf2_shift2 + shiftval3 * bin1_pdf2_shift2;

    // And we have shift 1 and 3 for PDFs 3:
    double bin0_pdf3_shift1 = bin0_pdf3 - shiftval1 * bin0_pdf3;
    double bin1_pdf3_shift1 = bin1_pdf3 + (shiftval1 * bin0_pdf3) - (shiftval1 * bin1_pdf3);
    // As shift 3 has a negative value, the shift goes the other way
    double bin0_pdf3_shift3 = bin0_pdf3_shift1 - (shiftval3 * bin1_pdf3_shift1) + (shiftval3 * bin0_pdf3_shift1);
    double bin1_pdf3_shift3 = bin1_pdf3_shift1 + shiftval3 * bin1_pdf3_shift1;

    // Now we can sum the total MC bin contents:
    double bin0_totmc = (pdfnorm1 * bin0_pdf1_shift1) + (pdfnorm2 * bin0_pdf2_shift3) + (pdfnorm3 * bin0_pdf3_shift3);
    double bin1_totmc = (pdfnorm1 * bin1_pdf1_shift1) + (pdfnorm2 * bin1_pdf2_shift3) + (pdfnorm3 * bin1_pdf3_shift3);

    // And calculate the LLH! Summing -data*log(mc) + mc for each bin
    double calcllh = (-bin0_data * log(bin0_totmc) + bin0_totmc) + (-bin1_data * log(bin1_totmc) + bin1_totmc);

    REQUIRE_THAT(llh, Catch::Matchers::WithinAbs(calcllh, 0.0001));
    REQUIRE_THAT(llh, Catch::Matchers::WithinAbs(0.371851, 0.0001));
  }
}
