#include <catch2/catch_approx.hpp>
#include <catch2/catch_all.hpp>
#include <BinnedNLLH.h>
#include <Gaussian.h>
#include <DistTools.h>
#include <OXSXDataSet.h>
#include <iostream>

TEST_CASE("Binned NLLH, 3 rates no systematics"){
    Gaussian gaus1(1, 5);
    Gaussian gaus2(1, 8);
    Gaussian gaus3(1, 10);

    AxisCollection axes;
    axes.AddAxis(BinAxis("axis1", -40, 40 , 200));

    BinnedED pdf1("a", DistTools::ToHist(gaus1, axes));
    BinnedED pdf2("b", DistTools::ToHist(gaus2, axes));
    BinnedED pdf3("c", DistTools::ToHist(gaus3, axes));


    size_t centralBin = pdf1.FindBin(std::vector<double>(1,0));
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
    SECTION("Correct Probability"){
        double sumLogProb = -log(prob1 + prob2 + prob3);
        double sumNorm    = 3;
        
        ParameterDict params;
        params["a"] = 1;
        params["b"] = 1;
        params["c"] = 1;
        lh.SetParameters(params);
        REQUIRE(lh.GetParameters() == params);
        REQUIRE(lh.Evaluate() == Catch::Approx(sumNorm + sumLogProb));
    }
    SECTION("Correct Probability with constraint"){
        lh.SetConstraint("a", 3, 1);

        double sumLogProb = -log(prob1 + prob2 + prob3);
        double sumNorm    = 3;
        double constraint = 2;

        ParameterDict params;
        params["a"] = 1;
        params["b"] = 1;
        params["c"] = 1;
        lh.SetParameters(params);
        REQUIRE(lh.Evaluate() == Catch::Approx(sumNorm + sumLogProb + constraint));
    }
    SECTION("Correct Probability with asymmetric constraint"){
        lh.SetConstraint("b", 5, 1, 2);

        double sumLogProb = -log(prob1 + prob2 + prob3);
        double sumNorm    = 3;
        double constraint = 8;

        ParameterDict params;
        params["a"] = 1;
        params["b"] = 1;
        params["c"] = 1;
        lh.SetParameters(params);
        REQUIRE(lh.Evaluate() == Catch::Approx(sumNorm + sumLogProb + constraint));
    }
    SECTION("Correct Probability with asymmetric constraint 2"){
        lh.SetConstraint("b", -3, 1, 2);

        double sumLogProb = -log(prob1 + prob2 + prob3);
        double sumNorm    = 3;
        double constraint = 2;

        ParameterDict params;
        params["a"] = 1;
        params["b"] = 1;
        params["c"] = 1;
        lh.SetParameters(params);
        REQUIRE(lh.Evaluate() == Catch::Approx(sumNorm + sumLogProb + constraint));
    }

    std::vector<int> genRates(3, pow(10,6));
    BinnedNLLH lh2;
    lh2.SetBarlowBeeston(true);
    lh2.AddPdf(pdf1, genRates.at(0));
    lh2.AddPdf(pdf2, genRates.at(1));
    lh2.AddPdf(pdf3, genRates.at(2));
    lh2.SetDataSet(&data);
    lh2.RegisterFitComponents();

    SECTION("Correct Probability with Barlow Beeston"){
      double betaPen    = 0;
      double sumLogProb = 0;
      double sumNorm    = 0;

      for(unsigned int i=0; i<pdf1.GetNBins(); i++){

	double binprob1 = pdf1.GetBinContent(i);
	double binprob2 = pdf2.GetBinContent(i);
	double binprob3 = pdf3.GetBinContent(i);	
	double binprob  = binprob1 + binprob2 + binprob3;

	double dat=0;
	if(i==centralBin)
	  dat = 1;

	double sig1     = sqrt(1/(float)genRates.at(0));
	double sig2     = sqrt(1/(float)genRates.at(1));
	double sig3     = sqrt(1/(float)genRates.at(2));	
	double sig      = sqrt(sig1*sig1 + sig2*sig2 + sig3*sig3)/binprob;
	double beta     = (-(binprob*sig*sig - 1) + sqrt((binprob*sig*sig - 1)*(binprob*sig*sig - 1) + 4*dat*sig*sig))/2;

	betaPen         += (beta-1)*(beta-1)/(2*sig*sig);
	sumLogProb      += -dat*log(beta*(binprob));
	sumNorm         += beta*(binprob);
      }

      ParameterDict params;
      params["a"] = 1;
      params["b"] = 1;
      params["c"] = 1;
      lh2.SetParameters(params);
      REQUIRE(lh2.Evaluate() == Catch::Approx(sumNorm + sumLogProb + betaPen));
    }
    SECTION("Correct Probability with Barlow Beeston and constraint"){
      lh2.SetConstraint("a", 3, 1);
    
      double betaPen    = 0;
      double sumLogProb = 0;
      double sumNorm    = 0;
      double constraint = 2;

      for(unsigned int i=0; i<pdf1.GetNBins(); i++){

	double binprob1 = pdf1.GetBinContent(i);
        double binprob2 = pdf2.GetBinContent(i);
        double binprob3 = pdf3.GetBinContent(i);
        double binprob  = binprob1 + binprob2 + binprob3;

        double dat=0;
        if(i==centralBin)
          dat = 1;

        double sig1     = sqrt(1/(float)genRates.at(0));
        double sig2     = sqrt(1/(float)genRates.at(1));
        double sig3     = sqrt(1/(float)genRates.at(2));
        double sig      = sqrt(sig1*sig1 + sig2*sig2 + sig3*sig3)/binprob;
        double beta     = (-(binprob*sig*sig - 1) + sqrt((binprob*sig*sig - 1)*(binprob*sig*sig - 1) + 4*dat*sig*sig))/2;

        betaPen         += (beta-1)*(beta-1)/(2*sig*sig);
        sumLogProb      += -dat*log(beta*(binprob));
        sumNorm         += beta*(binprob);
      }

      ParameterDict params;
      params["a"] = 1;
      params["b"] = 1;
      params["c"] = 1;
      lh2.SetParameters(params);
      REQUIRE(lh2.Evaluate() == Catch::Approx(sumNorm + sumLogProb + betaPen + constraint));
    }
    SECTION("Correct probability with Barlow Beeston and asymmetric constraint"){
      lh2.SetConstraint("b", 5, 1, 2);

      double betaPen    = 0;
      double sumLogProb = 0;
      double sumNorm    = 0;
      double constraint = 8;

      for(unsigned int i=0; i<pdf1.GetNBins(); i++){

        double binprob1 = pdf1.GetBinContent(i);
        double binprob2 = pdf2.GetBinContent(i);
        double binprob3 = pdf3.GetBinContent(i);
        double binprob  = binprob1 + binprob2 + binprob3;

        double dat=0;
        if(i==centralBin)
          dat = 1;

        double sig1     = sqrt(1/(float)genRates.at(0));
        double sig2     = sqrt(1/(float)genRates.at(1));
        double sig3     = sqrt(1/(float)genRates.at(2));
        double sig      = sqrt(sig1*sig1 + sig2*sig2 + sig3*sig3)/binprob;
        double beta     = (-(binprob*sig*sig - 1) + sqrt((binprob*sig*sig - 1)*(binprob*sig*sig - 1) + 4*dat*sig*sig))/2;

        betaPen         += (beta-1)*(beta-1)/(2*sig*sig);
        sumLogProb      += -dat*log(beta*(binprob));
        sumNorm         += beta*(binprob);
      }

      ParameterDict params;
      params["a"] = 1;
      params["b"] = 1;
      params["c"] = 1;
      lh2.SetParameters(params);
      REQUIRE(lh2.Evaluate() == Catch::Approx(sumNorm + sumLogProb + betaPen + constraint));
    }
    SECTION("Correct Probability with Barlow Beeston and constraint 2"){
      lh2.SetConstraint("b", -3, 1, 2);

      double betaPen    = 0;
      double sumLogProb = 0;
      double sumNorm    = 0;
      double constraint = 2;

      for(unsigned int i=0; i<pdf1.GetNBins(); i++){

	double binprob1 = pdf1.GetBinContent(i);
	double binprob2 = pdf2.GetBinContent(i);
	double binprob3 = pdf3.GetBinContent(i);
	double binprob  = binprob1 + binprob2 + binprob3;

	double dat=0;
	if(i==centralBin)
          dat = 1;

	double sig1     = sqrt(1/(float)genRates.at(0));
	double sig2     = sqrt(1/(float)genRates.at(1));
	double sig3     = sqrt(1/(float)genRates.at(2));
	double sig      = sqrt(sig1*sig1 + sig2*sig2 + sig3*sig3)/binprob;
	double beta     = (-(binprob*sig*sig - 1) + sqrt((binprob*sig*sig - 1)*(binprob*sig*sig - 1) + 4*dat*sig*sig))/2;

	betaPen         += (beta-1)*(beta-1)/(2*sig*sig);
	sumLogProb      += -dat*log(beta*(binprob));
	sumNorm         += beta*(binprob);
      }

      ParameterDict params;
      params["a"] = 1;
      params["b"] = 1;
      params["c"] = 1;
      lh2.SetParameters(params);
      REQUIRE(lh2.Evaluate() == Catch::Approx(sumNorm + sumLogProb + betaPen + constraint));
    }
}
