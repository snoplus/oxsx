//////////////////////////////////////////////////////////////
//   Simple example how to fit a fake data distibution      //
//   with a PDF and one global systematic in one dimension  //
//////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <math.h>
#include <Rand.h>
#include <fstream>

#include <TH1D.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TPaveStats.h>

#include <BinnedED.h>
#include <BinnedEDGenerator.h>
#include <SystematicManager.h>
#include <BinnedNLLH.h>
#include <FitResult.h>
#include <Minuit.h>
#include <DistTools.h>
#include <Minuit.h>
#include <Convolution.h>
#include <Scale.h>
#include <BoolCut.h>
#include <BoxCut.h>
#include <Gaussian.h>
#include <ParameterDict.h>
#include <ContainerTools.hpp>

TH1D *diffHist(TH1D *h1, TH1D *h2);

int main()
{

    Rand::SetSeed(0);
    AxisCollection axes;
    axes.AddAxis(BinAxis("axis1", 0, 50, 200));

    /////////////////////
    // 1. Set up a PDF //
    /////////////////////
    Gaussian gausPdf(15, 0.5);
    BinnedED pdf1("a_mc", DistTools::ToHist(gausPdf, axes));
    pdf1.SetObservables({"axis1"});
    pdf1.Scale(100000);

    pdf1.Normalise();
    std::vector<BinnedED> mcPdfs;
    mcPdfs.push_back(pdf1);

    // Padding MC PDFs with small numbers
    std::cout << "Padding MC Now" << std::endl;
    for (size_t i = 0; i < mcPdfs.size(); i++)
    {
        mcPdfs.at(i).AddPadding();
    }

    /////////////////////////
    // 2. Set up Fake Data //
    /////////////////////////
    Gaussian gausData(20, 1.0);
    BinnedED pdf3("a_data", DistTools::ToHist(gausData, axes));
    pdf3.SetObservables({"axis1"});

    std::vector<BinnedED> dataPdfs;
    dataPdfs.push_back(pdf3);

    // Padding data distributions with small numbers
    std::cout << "Padding Data Now" << std::endl;
    for (size_t i = 0; i < dataPdfs.size(); i++)
    {
        dataPdfs.at(i).AddPadding();
    }

    BinnedEDGenerator dataGen;
    dataGen.SetPdfs(dataPdfs);
    std::vector<double> rates(1, 100000);
    dataGen.SetRates(rates);

    // BinnedED fakeData= dataGen.ExpectedRatesED();
    BinnedED fakeData = dataGen.PoissonFluctuatedED();
    std::cout << "Fake data set up\n";

    ///////////////////////
    // 3. The systematic //
    ///////////////////////
    ObsSet obsSet("axis1");

    Convolution *conv_a = new Convolution("conv_a");
    Gaussian *gaus_a = new Gaussian(0, 1, "gaus_a");
    gaus_a->RenameParameter("means_0", "gaus_a_1");
    gaus_a->RenameParameter("stddevs_0", "gaus_a_2");

    conv_a->SetFunction(gaus_a);

    conv_a->SetAxes(axes);
    conv_a->SetTransformationObs(obsSet);
    conv_a->SetDistributionObs(obsSet);
    conv_a->Construct();

    std::cout << "Convolution constructed\n";
    ////////////////////////////////////
    // 4. Setting optimisation limits //
    ////////////////////////////////////
    ParameterDict minima;
    minima["a_mc"] = 10;
    minima["gaus_a_1"] = -15;
    minima["gaus_a_2"] = 0;

    ParameterDict maxima;
    maxima["a_mc"] = 200000;
    maxima["gaus_a_1"] = 15;
    maxima["gaus_a_2"] = 1;

    ParameterDict initialval;
    initialval["a_mc"] = 90000;
    initialval["gaus_a_1"] = 4.;
    initialval["gaus_a_2"] = 1.;

    ParameterDict initialerr;
    initialerr["a_mc"] = 0.1 * initialval["a_mc"];
    initialerr["gaus_a_1"] = 0.1 * initialval["gaus_a_1"];
    initialerr["gaus_a_2"] = 0.1 * initialval["gaus_a_2"];

    std::cout << "Limits set\n";
    //////////////////////////////////////////
    // 5. Setting up likelihood, including  //
    //    PDF, fake data and systematic     //
    //////////////////////////////////////////
    int BuffLow = 20;
    int BuffHigh = 20;

    BinnedNLLH lh;
    lh.SetBufferAsOverflow(false);
    lh.SetBuffer("axis1", BuffLow, BuffHigh);
    // Initialise with the data set
    lh.SetDataDist(fakeData);
    // Add systematics to the global group "", will be applied to all distributions.
    lh.AddSystematic(conv_a);
    // Associate EDs to lh, systematics will be applied
    lh.AddPdf(mcPdfs.at(0));

    std::cout << "BinnedNLLH set up\n";
    ////////////
    // 6. Fit //
    ////////////
    Minuit min;
    min.SetMethod("Simplex");
    min.SetMaxCalls(10000000);
    min.SetMinima(minima);
    min.SetMaxima(maxima);
    min.SetInitialValues(initialval);
    min.SetInitialErrors(initialerr);

    std::cout << "About to Fit" << std::endl;
    FitResult result = min.Optimise(&lh);
    result.SetPrintPrecision(4);
    result.Print();
    ParameterDict bestResult = result.GetBestFit();

    ////////////////////
    // 7. Plot Result //
    ////////////////////
    BinnedED BiHolder = mcPdfs.at(0);
    BinnedED BiResult;

    Convolution BiSmearer("aConv");
    BiSmearer.SetFunction(new Gaussian(bestResult.at("gaus_a_1"), bestResult.at("gaus_a_2")));
    BiSmearer.SetAxes(axes);
    BiSmearer.SetTransformationObs(obsSet);
    BiSmearer.SetDistributionObs(obsSet);
    BiSmearer.Construct();

    BiResult = BiSmearer(BiHolder);
    BiResult.Scale(bestResult.at("a_mc"));

    TH1D fakeDataHist;
    TH1D BiFit;
    BiFit.Sumw2();
    fakeDataHist = DistTools::ToTH1D(fakeData);
    BiFit = DistTools::ToTH1D(BiResult);

    TH1D FullFit("FullFit", "", BiFit.GetNbinsX(), BiFit.GetXaxis()->GetXmin(), BiFit.GetXaxis()->GetXmax());
    FullFit.Sumw2();

    FullFit.Add(&BiFit);

    TLegend *leg = new TLegend(0.8, 0.7, 1, 0.9);
    leg->AddEntry(&fakeDataHist, "FakeData", "lf");
    leg->AddEntry(&BiFit, "a fit result", "lf");

    TCanvas *diff = new TCanvas("diff", "", 800, 800);
    diff->cd();

    // -------------- Top panel
    gStyle->SetOptStat(kFALSE);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
    pad1->Draw();
    pad1->SetLogy(1);
    pad1->cd();
    pad1->SetGrid(kTRUE);
    pad1->SetBottomMargin(0.00);
    gPad->RedrawAxis();

    fakeDataHist.SetTitle("Independent Systematic Fit");
    fakeDataHist.GetYaxis()->SetTitle(Form("Counts"));
    fakeDataHist.GetYaxis()->SetTitleOffset(1.5);
    fakeDataHist.SetFillColorAlpha(kGreen, 0.5);

    fakeDataHist.Draw();
    FullFit.SetFillColorAlpha(kRed, 0.5);
    BiFit.SetLineColor(kRed);
    BiFit.SetLineWidth(3);
    BiFit.Draw("same e");

    TH1D hist1 = DistTools::ToTH1D(pdf1);
    hist1.Scale(4000);

    hist1.SetFillColorAlpha(kRed, 0.5);
    hist1.SetLineWidth(2);
    hist1.Draw("same");

    leg->AddEntry(&hist1, "a pdf", "lf");
    leg->Draw();

    TPaveText pt(0.7, 0.2, 1.0, 0.6, "NDC");
    pt.AddText(Form("a norm = %.2f", bestResult["a_mc"]));
    pt.AddText(Form("a conv mean = %.2f", bestResult["gaus_a_1"]));
    pt.AddText(Form("a conv RMS= %.2f", bestResult["gaus_a_2"]));
    pt.SetFillColor(kWhite);
    pt.SetShadowColor(kWhite);
    pt.Draw();
    diff->cd();

    // -------------- Bottom panel
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
    pad2->SetTopMargin(0.00);
    pad2->Draw();
    pad2->cd();
    pad2->SetBottomMargin(0.3);
    pad2->SetGrid(kTRUE);
    gStyle->SetOptStat(kFALSE);

    TH1D *fracDiff = diffHist(&fakeDataHist, &FullFit);
    fracDiff->SetLineWidth(2);

    fracDiff->GetXaxis()->SetTitle("axis1");
    fracDiff->GetYaxis()->SetTitle("Fit / Fake Data");
    fracDiff->GetYaxis()->SetTitleOffset(0.5);
    fracDiff->GetXaxis()->SetLabelSize(0.1);
    fracDiff->GetXaxis()->SetTitleSize(0.1);
    fracDiff->GetYaxis()->SetLabelSize(0.1);
    fracDiff->GetYaxis()->SetTitleSize(0.1);

    fracDiff->SetMaximum(2);
    fracDiff->SetMinimum(0);
    fracDiff->Draw();

    diff->Print("simpleSystematicFitExample.png");

    return 0;
}

TH1D *diffHist(TH1D *h1, TH1D *h2)
{
    double minBin = h1->GetXaxis()->GetXmin();
    double maxBin = h1->GetXaxis()->GetXmax();
    double numOfBins = h1->GetNbinsX();

    TH1D *rhist = new TH1D("rhist", "", numOfBins, minBin, maxBin);
    for (double i = 0; i < numOfBins; i++)
    {
        double h1cont = h1->GetBinContent(i);
        double h2cont = h2->GetBinContent(i);
        double weight;
        double error;
        if (h1cont != 0 && h2cont != 0)
        {
            weight = h2cont / h1cont;
            error = weight * sqrt(1 / sqrt(h2cont) + 1 / sqrt(h1cont));
        }
        else
        {
            // Don't show that bin on the plot
            weight = 10000000;
            error = 0;
        }
        rhist->SetBinContent(i, weight);
        rhist->SetBinError(i, error);
    }
    return rhist;
}
