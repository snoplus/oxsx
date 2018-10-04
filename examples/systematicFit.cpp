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

TH1D* diffHist(TH1D * h1,TH1D * h2);
void padPDFs(std::vector<BinnedED>& binnedEDList);

void function1();
void function2();
void function3();

int main(int argc, char *argv[]){
    // Using groups of systematics.
    function1();
    // Using global of systematics.
    function2();
    // Using both global systematics with groups of systematics.
    function3();
    return 0;
}

void
function1(){

    Rand::SetSeed(0);
    Gaussian gaus1(10, 0.65);
    Gaussian gaus2(15, 0.5);
    Gaussian gaus3(20, 0.5);
    Gaussian gaus4(30, 1.0);

    AxisCollection axes;
    axes.AddAxis(BinAxis("axis1", 0, 50 ,200));

    // Making PDFs. 
    BinnedED pdf3("a_data", DistTools::ToHist(gaus1, axes));
    BinnedED pdf1("a_mc", DistTools::ToHist(gaus2, axes));
    // Making fake data.
    BinnedED pdf4("b_data", DistTools::ToHist(gaus4, axes));
    BinnedED pdf2("b_mc", DistTools::ToHist(gaus3, axes));

    pdf1.SetObservables(0);
    pdf2.SetObservables(0);    
    pdf3.SetObservables(0);
    pdf4.SetObservables(0);

    pdf1.Scale(100000);
    pdf2.Scale(100000);    

    std::vector<BinnedED> dataPdfs;
    dataPdfs.push_back(pdf3);
    dataPdfs.push_back(pdf4);
    padPDFs(dataPdfs);

    BinnedEDGenerator dataGen;
    dataGen.SetPdfs(dataPdfs);
    std::vector<double> rates(2,100000);
    dataGen.SetRates(rates);

    // BinnedED fakeData= dataGen.ExpectedRatesED();
    // Fluctuate data.
    BinnedED fakeData= dataGen.PoissonFluctuatedED();

    pdf1.Normalise();
    pdf2.Normalise();    
    std::vector<BinnedED> mcPdfs;
    mcPdfs.push_back(pdf1);
    mcPdfs.push_back(pdf2);

    padPDFs(mcPdfs);

    ObsSet  obsSet(0);

    Convolution* conv_a = new Convolution("conv_a");
    Gaussian* gaus_a = new Gaussian(0,1,"gaus_a"); 
    gaus_a->RenameParameter("means_0","gaus_a_1");
    gaus_a->RenameParameter("stddevs_0","gaus_a_2");

    conv_a->SetFunction(gaus_a);

    conv_a->SetAxes(axes);
    conv_a->SetTransformationObs(obsSet);
    conv_a->SetDistributionObs(obsSet);
    conv_a->Construct();

    Convolution* conv_b = new Convolution("conv_b");
    Gaussian * gaus_b = new Gaussian(0,1,"gaus_b"); 
    gaus_b->RenameParameter("means_0","gaus_b_1");
    gaus_b->RenameParameter("stddevs_0","gaus_b_2");
    conv_b->SetFunction(gaus_b);
    conv_b->SetAxes(axes);
    conv_b->SetTransformationObs(obsSet);
    conv_b->SetDistributionObs(obsSet);
    conv_b->Construct();

    // Setting optimisation limits
    ParameterDict minima;
    minima["a_mc_norm"] = 10; 
    minima["b_mc_norm"] = 10; 
    minima["gaus_a_1" ] = -15;
    minima["gaus_a_2" ] = 0;  
    minima["gaus_b_1" ] = -15;
    minima["gaus_b_2" ] = 0;  

    ParameterDict maxima;
    maxima["a_mc_norm"] = 200000;
    maxima["b_mc_norm"] = 200000;
    maxima["gaus_a_1" ] = 15;    
    maxima["gaus_a_2" ] = 1;     
    maxima["gaus_b_1" ] = 16.;   
    maxima["gaus_b_2" ] = 1;     


    // ParameterDict initialval;
    // Rand rand;
    // initialval["a_mc_norm"] = rand.UniformRange(minima["a_mc_norm"],maxima["a_mc_norm"]); 
    // initialval["b_mc_norm"] = rand.UniformRange(minima["b_mc_norm"],maxima["b_mc_norm"]); 
    // initialval["gaus_a_1" ] = rand.UniformRange(minima["gaus_a_1" ],maxima["gaus_a_1" ]); 
    // initialval["gaus_a_2" ] = rand.UniformRange(minima["gaus_a_2" ],maxima["gaus_a_2" ]); 
    // initialval["gaus_b_1" ] = rand.UniformRange(minima["gaus_b_1" ],maxima["gaus_b_1" ]); 
    // initialval["gaus_b_2" ] = rand.UniformRange(minima["gaus_b_2" ],maxima["gaus_b_2" ]); 

    ParameterDict initialval;
    initialval["a_mc_norm"] = 90000; 
    initialval["b_mc_norm"] = 90000; 
    initialval["gaus_a_1"]  = -4.;   
    initialval["gaus_a_2"]  = 1.;    
    initialval["gaus_b_1"]  = 9.;    
    initialval["gaus_b_2"]  = 1.;    

    ParameterDict initialerr;
    initialerr["a_mc_norm"] = 0.1*initialval["a_mc_norm"];
    initialerr["b_mc_norm"] = 0.1*initialval["b_mc_norm"];
    initialerr["gaus_a_1" ] = 0.1*initialval["gaus_a_1"]; 
    initialerr["gaus_a_2" ] = 0.1*initialval["gaus_a_2"]; 
    initialerr["gaus_b_1" ] = 0.1*initialval["gaus_b_1"]; 
    initialerr["gaus_b_2" ] = 0.1*initialval["gaus_b_2"]; 

    //Setting up likelihood.
    int BuffLow  = 20;
    int BuffHigh = 20;

    BinnedNLLH lh; 
    lh.SetBufferAsOverflow(false);
    lh.SetBuffer(0,BuffLow,BuffHigh);
    lh.SetDataDist(fakeData); // initialise with the data set
    
    // Add systematics to groups. The order in which the systematic is added to
    // the group is the same as the order in which they are applied.
    lh.AddSystematic(conv_a,"aGroup");
    lh.AddSystematic(conv_b,"bGroup");

    // Associate EDs with a vector of groups. The order of the vector of groups
    // is the same as the order in which they are applied.
    lh.AddDist(mcPdfs.at(0),std::vector<std::string>(1,"aGroup"));
    lh.AddDist(mcPdfs.at(1),std::vector<std::string>(1,"bGroup"));

    Minuit min;
    min.SetMethod("Simplex");
    min.SetMaxCalls(10000000);
    min.SetMinima(minima);
    min.SetMaxima(maxima);
    min.SetInitialValues(initialval);
    min.SetInitialErrors(initialerr);

    std::cout << "function1 : About to Fit" << std::endl;
    FitResult result = min.Optimise(&lh);
    result.SetPrintPrecision(4);
    result.Print();
    ParameterDict bestResult = result.GetBestFit();

    // Plot Result
    {
        BinnedED BiHolder = mcPdfs.at(0);
        BinnedED PoHolder = mcPdfs.at(1);
        BinnedED BiResult;
        BinnedED PoResult;

        Convolution BiSmearer("aConv");
        BiSmearer.SetFunction(new Gaussian(bestResult.at("gaus_a_1"),bestResult.at("gaus_a_2")));
        BiSmearer.SetAxes(axes);
        BiSmearer.SetTransformationObs(obsSet);
        BiSmearer.SetDistributionObs(obsSet);
        BiSmearer.Construct();

        Convolution PoSmearer("bConv");
        PoSmearer.SetFunction(new Gaussian(bestResult.at("gaus_b_1"),bestResult.at("gaus_b_2")));
        PoSmearer.SetAxes(axes);
        PoSmearer.SetTransformationObs(obsSet);
        PoSmearer.SetDistributionObs(obsSet);
        PoSmearer.Construct();

        PoResult = PoSmearer( PoHolder );
        BiResult = BiSmearer( BiHolder );

        BiResult.Scale(bestResult.at("a_mc_norm"));
        PoResult.Scale(bestResult.at("b_mc_norm"));

        TH1D fakeDataHist;
        TH1D BiFit;
        TH1D PoFit;
        BiFit.Sumw2();
        PoFit.Sumw2();
        fakeDataHist = DistTools::ToTH1D(fakeData);
        BiFit= DistTools::ToTH1D(BiResult);
        PoFit= DistTools::ToTH1D(PoResult);

        TH1D FullFit("FullFit","",BiFit.GetNbinsX(),BiFit.GetXaxis()->GetXmin(),BiFit.GetXaxis()->GetXmax());
        FullFit.Sumw2();

        FullFit.Add(&BiFit,&PoFit);

        TLegend* leg =new TLegend(0.8,0.7,1,0.9); 
        leg->AddEntry(&fakeDataHist,"FakeData","lf"); 
        leg->AddEntry(&BiFit,"a fit result","lf"); 
        leg->AddEntry(&PoFit,"b fit result","lf"); 

        TCanvas* diff = new TCanvas("diff","",800,800);
        diff->cd(); 
        // -------------- Top panel 
        gStyle->SetOptStat(kFALSE);  
        TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1); 
        pad1->Draw(); 
        pad1->SetLogy(1); 
        pad1->cd(); 
        pad1->SetGrid(kTRUE); 
        pad1->SetBottomMargin(0.00); 
        gPad->RedrawAxis();  

        fakeDataHist.SetTitle("Independent Systematic Fit");
        fakeDataHist.GetYaxis()->SetTitle(Form("Counts"));
        fakeDataHist.GetYaxis()->SetTitleOffset(1.5); 
        fakeDataHist.SetFillColorAlpha(kGreen,0.5); 

        fakeDataHist.Draw(); 
        FullFit.SetFillColorAlpha(kRed,0.5); 
        BiFit.SetLineColor(kRed);
        BiFit.SetLineWidth(3);
        BiFit.Draw("same e"); 
        PoFit.SetLineColor(kBlue);
        PoFit.SetLineWidth(3);
        PoFit.Draw("same e"); 

        TH1D hist1 =  DistTools::ToTH1D(pdf1);
        TH1D hist2 =  DistTools::ToTH1D(pdf2);
        hist1.Scale(4000);
        hist2.Scale(4000);

        hist1.SetFillColorAlpha(kRed,0.5); 
        hist1.SetLineWidth(2);
        hist1.Draw("same");
        hist2.SetFillColorAlpha(kBlue,0.5); 
        hist2.SetLineWidth(2);
        hist2.Draw("same");

        leg->AddEntry(&hist1,"a pdf","lf"); 
        leg->AddEntry(&hist2,"b pdf","lf"); 
        leg->Draw(); 

        TPaveText pt(0.7,0.2,1.0,0.6,"NDC");
        pt.AddText(Form("a norm = %.2f" ,bestResult["a_mc_norm"]));
        pt.AddText(Form("b norm = %.2f",bestResult["b_mc_norm"]));
        pt.AddText(Form("a conv mean = %.2f",bestResult["gaus_a_1"]));
        pt.AddText(Form("a conv RMS= %.2f",bestResult["gaus_a_2"]));
        pt.AddText(Form("b conv mean = %.2f",bestResult["gaus_b_1"]));
        pt.AddText(Form("b conv RMS = %.2f",bestResult["gaus_b_2"]));
        pt.SetFillColor(kWhite);
        pt.SetShadowColor(kWhite);
        pt.Draw();
        diff->cd(); 

        // -------------- Bottom panel 
        TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.3); 
        pad2->SetTopMargin(0.00); 
        pad2->Draw(); 
        pad2->cd(); 
        pad2->SetBottomMargin(0.3); 
        pad2->SetGrid(kTRUE); 
        gStyle->SetOptStat(kFALSE);  

        TH1D * fracDiff= diffHist(&fakeDataHist,&FullFit); 
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

        diff->Print("systematicFitExample_withGroups.png"); 
        delete diff;
    }
}

void
function2(){

    Rand::SetSeed(0);
    AxisCollection axes;
    axes.AddAxis(BinAxis("axis1", 0, 50 ,200));
    ObsSet  obsSet(0);

    // Setting up artificial scale with a scaleFactor = 1.5.
    Scale* scale = new Scale("scale");
    scale->SetScaleFactor(1.5);
    scale->SetAxes(axes);
    scale->SetTransformationObs(obsSet);
    scale->SetDistributionObs(obsSet);
    scale->Construct();

    Gaussian gaus1(10, 0.65);
    Gaussian gaus2(25, 0.5);

    BinnedED pdf1("a_mc", DistTools::ToHist(gaus1, axes));
    BinnedED pdf2("b_mc", DistTools::ToHist(gaus2, axes));
    BinnedED pdf3 = scale->operator()(pdf1);
    BinnedED pdf4 = scale->operator()(pdf2);

    pdf1.SetObservables(0);
    pdf2.SetObservables(0);    
    pdf3.SetObservables(0);
    pdf4.SetObservables(0);

    pdf1.Scale(100000);
    pdf2.Scale(100000);    

    std::vector<BinnedED> dataPdfs;
    dataPdfs.push_back(pdf3);
    dataPdfs.push_back(pdf4);
    padPDFs(dataPdfs);

    BinnedEDGenerator dataGen;
    dataGen.SetPdfs(dataPdfs);
    std::vector<double> rates(2,100000);
    dataGen.SetRates(rates);

    // BinnedED fakeData= dataGen.ExpectedRatesED();
    // Fluctuate data.
    BinnedED fakeData= dataGen.PoissonFluctuatedED();

    pdf1.Normalise();
    pdf2.Normalise();    
    std::vector<BinnedED> mcPdfs;
    mcPdfs.push_back(pdf1);
    mcPdfs.push_back(pdf2);

    padPDFs(mcPdfs);

    // Setting optimisation limits
    ParameterDict minima;
    minima["a_mc_norm"] = 10; 
    minima["b_mc_norm"] = 10; 
    minima["scaleFactor" ] = 0.1;

    ParameterDict maxima;
    maxima["a_mc_norm"] = 200000;
    maxima["b_mc_norm"] = 200000;
    maxima["scaleFactor" ] = 2.0;    

    // ParameterDict initialval;
    // Rand rand;
    // initialval["a_mc_norm"] = rand.UniformRange(minima["a_mc_norm"],maxima["a_mc_norm"]); 
    // initialval["b_mc_norm"] = rand.UniformRange(minima["b_mc_norm"],maxima["b_mc_norm"]); 
    // initialval["gaus_a_1" ] = rand.UniformRange(minima["gaus_a_1" ],maxima["gaus_a_1" ]); 
    // initialval["gaus_a_2" ] = rand.UniformRange(minima["gaus_a_2" ],maxima["gaus_a_2" ]); 
    // initialval["gaus_b_1" ] = rand.UniformRange(minima["gaus_b_1" ],maxima["gaus_b_1" ]); 
    // initialval["gaus_b_2" ] = rand.UniformRange(minima["gaus_b_2" ],maxima["gaus_b_2" ]); 

    ParameterDict initialval;
    initialval["a_mc_norm"] = 90000; 
    initialval["b_mc_norm"] = 90000; 
    initialval["scaleFactor"]  = 1;   

    ParameterDict initialerr;
    initialerr["a_mc_norm"] = 0.1*initialval["a_mc_norm"];
    initialerr["b_mc_norm"] = 0.1*initialval["b_mc_norm"];
    initialerr["scaleFactor" ] = 0.1*initialval["scaleFactor"]; 

    //Setting up likelihood.
    int BuffLow  = 20;
    int BuffHigh = 20;

    BinnedNLLH lh; 
    lh.SetBufferAsOverflow(false);
    lh.SetBuffer(0,BuffLow,BuffHigh);
    lh.SetDataDist(fakeData); // initialise with the data set

    // Add global systematic. Order in which global systematics are applied is the same as the order in which they are applied.
    lh.AddSystematic(scale);

    // Add EDs to lh, by default any global systematic are applied.
    lh.AddDist(mcPdfs.at(0));
    lh.AddDist(mcPdfs.at(1));

    Minuit min;
    min.SetMethod("Simplex");
    min.SetMaxCalls(10000000);
    min.SetMinima(minima);
    min.SetMaxima(maxima);
    min.SetInitialValues(initialval);
    min.SetInitialErrors(initialerr);

    std::cout << "function2 : About to Fit" << std::endl;
    FitResult result = min.Optimise(&lh);
    result.SetPrintPrecision(4);
    result.Print();
    ParameterDict bestResult = result.GetBestFit();

    // Plot Result
    {
        BinnedED BiHolder = mcPdfs.at(0);
        BinnedED PoHolder = mcPdfs.at(1);
        BinnedED BiResult;
        BinnedED PoResult;

        Scale* scale = new Scale("scale");
        scale->SetScaleFactor(bestResult.at("scaleFactor"));
        scale->SetAxes(axes);
        scale->SetTransformationObs(obsSet);
        scale->SetDistributionObs(obsSet);
        scale->Construct();

        PoResult = scale->operator()( PoHolder );
        BiResult = scale->operator()( BiHolder );

        BiResult.Scale(bestResult.at("a_mc_norm"));
        PoResult.Scale(bestResult.at("b_mc_norm"));

        TH1D fakeDataHist;
        TH1D BiFit;
        TH1D PoFit;
        BiFit.Sumw2();
        PoFit.Sumw2();
        fakeDataHist = DistTools::ToTH1D(fakeData);
        BiFit= DistTools::ToTH1D(BiResult);
        PoFit= DistTools::ToTH1D(PoResult);

        TH1D FullFit("FullFit","",BiFit.GetNbinsX(),BiFit.GetXaxis()->GetXmin(),BiFit.GetXaxis()->GetXmax());
        FullFit.Sumw2();

        FullFit.Add(&BiFit,&PoFit);

        TLegend* leg =new TLegend(0.8,0.7,1,0.9); 
        leg->AddEntry(&fakeDataHist,"FakeData","lf"); 
        leg->AddEntry(&BiFit,"a fit result","lf"); 
        leg->AddEntry(&PoFit,"b fit result","lf"); 

        TCanvas* diff = new TCanvas("diff","",800,800);
        diff->cd(); 
        // -------------- Top panel 
        gStyle->SetOptStat(kFALSE);  
        TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1); 
        pad1->Draw(); 
        pad1->SetLogy(1); 
        pad1->cd(); 
        pad1->SetGrid(kTRUE); 
        pad1->SetBottomMargin(0.00); 
        gPad->RedrawAxis();  

        fakeDataHist.SetTitle("Independent Systematic Fit");
        fakeDataHist.GetYaxis()->SetTitle(Form("Counts"));
        fakeDataHist.GetYaxis()->SetTitleOffset(1.5); 
        fakeDataHist.SetFillColorAlpha(kGreen,0.5); 

        fakeDataHist.Draw(); 
        FullFit.SetFillColorAlpha(kRed,0.5); 
        BiFit.SetLineColor(kRed);
        BiFit.SetLineWidth(3);
        BiFit.Draw("same e"); 
        PoFit.SetLineColor(kBlue);
        PoFit.SetLineWidth(3);
        PoFit.Draw("same e"); 

        TH1D hist1 =  DistTools::ToTH1D(pdf1);
        TH1D hist2 =  DistTools::ToTH1D(pdf2);
        hist1.Scale(4000);
        hist2.Scale(4000);

        hist1.SetFillColorAlpha(kRed,0.5); 
        hist1.SetLineWidth(2);
        hist1.Draw("same");
        hist2.SetFillColorAlpha(kBlue,0.5); 
        hist2.SetLineWidth(2);
        hist2.Draw("same");

        leg->AddEntry(&hist1,"a pdf","lf"); 
        leg->AddEntry(&hist2,"b pdf","lf"); 
        leg->Draw(); 

        TPaveText pt(0.7,0.2,1.0,0.6,"NDC");
        pt.AddText(Form("a norm = %.2f" ,bestResult["a_mc_norm"]));
        pt.AddText(Form("b norm = %.2f",bestResult["b_mc_norm"]));
        pt.AddText(Form("scaleFactor = %.2f",bestResult["scaleFactor"]));
        pt.SetFillColor(kWhite);
        pt.SetShadowColor(kWhite);
        pt.Draw();
        diff->cd(); 

        // -------------- Bottom panel 
        TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.3); 
        pad2->SetTopMargin(0.00); 
        pad2->Draw(); 
        pad2->cd(); 
        pad2->SetBottomMargin(0.3); 
        pad2->SetGrid(kTRUE); 
        gStyle->SetOptStat(kFALSE);  

        TH1D * fracDiff= diffHist(&fakeDataHist,&FullFit); 
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

        diff->Print("systematicFitExample_noGroup.png"); 
        delete diff;
    }
}

void
function3(){

    Rand::SetSeed(0);
    AxisCollection axes;
    axes.AddAxis(BinAxis("axis1", 0, 75 ,300));
    ObsSet  obsSet(0);

    Gaussian gaus1(10, 0.5);
    Gaussian gaus2(25, 1);

    // Setting up artificial scale with a scaleFactor = 1.5.
    Scale* scale = new Scale("scale");
    scale->SetScaleFactor(1.5);
    scale->SetAxes(axes);
    scale->SetTransformationObs(obsSet);
    scale->SetDistributionObs(obsSet);
    scale->Construct();

    // Setting up artificial convolution with a mean = 0 and std = 1.
    Convolution* conv_a = new Convolution("conv_a");
    Gaussian* gaus_a = new Gaussian(0,1,"gaus_a"); 
    gaus_a->RenameParameter("means_0","gaus_a_1");
    gaus_a->RenameParameter("stddevs_0","gaus_a_2");

    conv_a->SetFunction(gaus_a);
    conv_a->SetAxes(axes);
    conv_a->SetTransformationObs(obsSet);
    conv_a->SetDistributionObs(obsSet);
    conv_a->Construct();

    // Setting up artificial convolution with a mean = 0 and std = 2.
    Convolution* conv_b = new Convolution("conv_b");
    Gaussian * gaus_b = new Gaussian(0,2,"gaus_b"); 
    gaus_b->RenameParameter("means_0","gaus_b_1");
    gaus_b->RenameParameter("stddevs_0","gaus_b_2");
    conv_b->SetFunction(gaus_b);
    conv_b->SetAxes(axes);
    conv_b->SetTransformationObs(obsSet);
    conv_b->SetDistributionObs(obsSet);
    conv_b->Construct();

    // Making PDFs. 
    BinnedED pdf1("a_mc", DistTools::ToHist(gaus1, axes));
    BinnedED pdf2("b_mc", DistTools::ToHist(gaus2, axes));

    // Making fake data. 
    // Order of the systematics is important and should be the same order as
    // added to the likelihood.
    BinnedED pdf3 = conv_a->operator()(scale->operator()(pdf1));
    BinnedED pdf4 = conv_b->operator()(scale->operator()(pdf2));

    pdf1.SetObservables(0);
    pdf2.SetObservables(0);    
    pdf3.SetObservables(0);
    pdf4.SetObservables(0);

    pdf1.Scale(100000);
    pdf2.Scale(100000);    

    std::vector<BinnedED> dataPdfs;
    dataPdfs.push_back(pdf3);
    dataPdfs.push_back(pdf4);
    padPDFs(dataPdfs);

    BinnedEDGenerator dataGen;
    dataGen.SetPdfs(dataPdfs);
    std::vector<double> rates(2,100000);
    dataGen.SetRates(rates);

    // BinnedED fakeData= dataGen.ExpectedRatesED();
    // Fluctuate data.
    BinnedED fakeData= dataGen.PoissonFluctuatedED();

    pdf1.Normalise();
    pdf2.Normalise();    
    std::vector<BinnedED> mcPdfs;
    mcPdfs.push_back(pdf1);
    mcPdfs.push_back(pdf2);

    padPDFs(mcPdfs);

    // Setting optimisation limits
    ParameterDict minima;
    minima["a_mc_norm"] = 10; 
    minima["b_mc_norm"] = 10; 
    minima["scaleFactor" ] = 0.1;
    minima["gaus_a_1" ] = -5;
    minima["gaus_a_2" ] = 0;  
    minima["gaus_b_1" ] = -5;
    minima["gaus_b_2" ] = 0;  

    ParameterDict maxima;
    maxima["a_mc_norm"] = 200000;
    maxima["b_mc_norm"] = 200000;
    maxima["scaleFactor" ] = 2.0;    
    maxima["gaus_a_1" ] = 5;    
    maxima["gaus_a_2" ] = 5;     
    maxima["gaus_b_1" ] = 5.;   
    maxima["gaus_b_2" ] = 5;     

    // ParameterDict initialval;
    // Rand rand;
    // initialval["a_mc_norm"] = rand.UniformRange(minima["a_mc_norm"],maxima["a_mc_norm"]); 
    // initialval["b_mc_norm"] = rand.UniformRange(minima["b_mc_norm"],maxima["b_mc_norm"]); 
    // initialval["gaus_a_1" ] = rand.UniformRange(minima["gaus_a_1" ],maxima["gaus_a_1" ]); 
    // initialval["gaus_a_2" ] = rand.UniformRange(minima["gaus_a_2" ],maxima["gaus_a_2" ]); 
    // initialval["gaus_b_1" ] = rand.UniformRange(minima["gaus_b_1" ],maxima["gaus_b_1" ]); 
    // initialval["gaus_b_2" ] = rand.UniformRange(minima["gaus_b_2" ],maxima["gaus_b_2" ]); 

    ParameterDict initialval;
    initialval["a_mc_norm"] = 90000; 
    initialval["b_mc_norm"] = 90000; 
    initialval["scaleFactor"]  = 1.4;   
    initialval["gaus_a_1"]  = 0.;   
    initialval["gaus_a_2"]  = 1;    
    initialval["gaus_b_1"]  = 0.;    
    initialval["gaus_b_2"]  = 2;    

    ParameterDict initialerr;
    initialerr["a_mc_norm"] = 0.1*initialval["a_mc_norm"];
    initialerr["b_mc_norm"] = 0.1*initialval["b_mc_norm"];
    initialerr["scaleFactor" ] = 0.1*initialval["scaleFactor"]; 
    initialerr["gaus_a_1" ] = 0.1*initialval["gaus_a_1"]; 
    initialerr["gaus_a_2" ] = 0.1*initialval["gaus_a_2"]; 
    initialerr["gaus_b_1" ] = 0.1*initialval["gaus_b_1"]; 
    initialerr["gaus_b_2" ] = 0.1*initialval["gaus_b_2"]; 

    //Setting up likelihood.
    int BuffLow  = 20;
    int BuffHigh = 20;

    BinnedNLLH lh; 
    lh.SetBufferAsOverflow(false);
    lh.SetBuffer(0,BuffLow,BuffHigh);
    lh.SetDataDist(fakeData); // initialise with the data set

    // Add the scale as a global systematic.
    // i.e. The scale is the first systematic applied to all distributions. 
    // In general global systematics are applied (in the order they are added)
    // before any other group systematics.
    lh.AddSystematic(scale);
    // Adding systematics to groups. The order in which the systematics are
    // added to the group is the same as the order there are applied.
    lh.AddSystematic(conv_a,"aGroup");
    lh.AddSystematic(conv_b,"bGroup");

    // Associate EDs with a vector of groups. The order of the vector dictates
    // the order in which (multiple) groups are applied. 
    lh.AddDist(mcPdfs.at(0),std::vector<std::string>(1,"aGroup"));
    lh.AddDist(mcPdfs.at(1),std::vector<std::string>(1,"bGroup"));

    Minuit min;
    min.SetMethod("Simplex");
    min.SetMaxCalls(10000000);
    min.SetMinima(minima);
    min.SetMaxima(maxima);
    min.SetInitialValues(initialval);
    min.SetInitialErrors(initialerr);

    std::cout << "function3 : About to Fit" << std::endl;
    FitResult result = min.Optimise(&lh);
    result.SetPrintPrecision(4);
    result.Print();
    ParameterDict bestResult = result.GetBestFit();

    // Plot Result
    {
        BinnedED BiHolder = mcPdfs.at(0);
        BinnedED PoHolder = mcPdfs.at(1);
        BinnedED BiResult;
        BinnedED PoResult;

        Scale* scale = new Scale("scale");
        scale->SetScaleFactor(bestResult.at("scaleFactor"));
        scale->SetAxes(axes);
        scale->SetTransformationObs(obsSet);
        scale->SetDistributionObs(obsSet);
        scale->Construct();

        Convolution BiSmearer("aConv");
        BiSmearer.SetFunction(new Gaussian(bestResult.at("gaus_a_1"),bestResult.at("gaus_a_2")));
        BiSmearer.SetAxes(axes);
        BiSmearer.SetTransformationObs(obsSet);
        BiSmearer.SetDistributionObs(obsSet);
        BiSmearer.Construct();

        Convolution PoSmearer("bConv");
        PoSmearer.SetFunction(new Gaussian(bestResult.at("gaus_b_1"),bestResult.at("gaus_b_2")));
        PoSmearer.SetAxes(axes);
        PoSmearer.SetTransformationObs(obsSet);
        PoSmearer.SetDistributionObs(obsSet);
        PoSmearer.Construct();

        BiResult = BiSmearer(scale->operator()( BiHolder ));
        PoResult = PoSmearer(scale->operator()( PoHolder ));

        BiResult.Scale(bestResult.at("a_mc_norm"));
        PoResult.Scale(bestResult.at("b_mc_norm"));

        TH1D fakeDataHist;
        TH1D BiFit;
        TH1D PoFit;
        BiFit.Sumw2();
        PoFit.Sumw2();
        fakeDataHist = DistTools::ToTH1D(fakeData);
        BiFit= DistTools::ToTH1D(BiResult);
        PoFit= DistTools::ToTH1D(PoResult);

        TH1D FullFit("FullFit","",BiFit.GetNbinsX(),BiFit.GetXaxis()->GetXmin(),BiFit.GetXaxis()->GetXmax());
        FullFit.Sumw2();

        FullFit.Add(&BiFit,&PoFit);

        TLegend* leg =new TLegend(0.8,0.7,1,0.9); 
        leg->AddEntry(&fakeDataHist,"FakeData","lf"); 
        leg->AddEntry(&BiFit,"a fit result","lf"); 
        leg->AddEntry(&PoFit,"b fit result","lf"); 

        TCanvas* diff = new TCanvas("diff","",800,800);
        diff->cd(); 
        // -------------- Top panel 
        gStyle->SetOptStat(kFALSE);  
        TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1); 
        pad1->Draw(); 
        pad1->SetLogy(1); 
        pad1->cd(); 
        pad1->SetGrid(kTRUE); 
        pad1->SetBottomMargin(0.00); 
        gPad->RedrawAxis();  

        fakeDataHist.SetTitle("Independent Systematic Fit");
        fakeDataHist.GetYaxis()->SetTitle(Form("Counts"));
        fakeDataHist.GetYaxis()->SetTitleOffset(1.5); 
        fakeDataHist.SetFillColorAlpha(kGreen,0.5); 

        fakeDataHist.Draw(); 
        FullFit.SetFillColorAlpha(kRed,0.5); 
        BiFit.SetLineColor(kRed);
        BiFit.SetLineWidth(3);
        BiFit.Draw("same e"); 
        PoFit.SetLineColor(kBlue);
        PoFit.SetLineWidth(3);
        PoFit.Draw("same e"); 

        TH1D hist1 =  DistTools::ToTH1D(pdf1);
        TH1D hist2 =  DistTools::ToTH1D(pdf2);
        hist1.Scale(4000);
        hist2.Scale(4000);

        hist1.SetFillColorAlpha(kRed,0.5); 
        hist1.SetLineWidth(2);
        hist1.Draw("same");
        hist2.SetFillColorAlpha(kBlue,0.5); 
        hist2.SetLineWidth(2);
        hist2.Draw("same");

        leg->AddEntry(&hist1,"a pdf","lf"); 
        leg->AddEntry(&hist2,"b pdf","lf"); 
        leg->Draw(); 

        TPaveText pt(0.7,0.2,1.0,0.7,"NDC");
        pt.AddText(Form("a norm = %.2f" ,bestResult["a_mc_norm"]));
        pt.AddText(Form("b norm = %.2f",bestResult["b_mc_norm"]));
        pt.AddText(Form("scaleFactor = %.2f",bestResult["scaleFactor"]));
        pt.AddText(Form("a conv mean = %.2f",bestResult["gaus_a_1"]));
        pt.AddText(Form("a conv RMS= %.2f",bestResult["gaus_a_2"]));
        pt.AddText(Form("b conv mean = %.2f",bestResult["gaus_b_1"]));
        pt.AddText(Form("b conv RMS = %.2f",bestResult["gaus_b_2"]));
        pt.SetFillColor(kWhite);
        pt.SetShadowColor(kWhite);
        pt.Draw();
        diff->cd(); 

        // -------------- Bottom panel 
        TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.3); 
        pad2->SetTopMargin(0.00); 
        pad2->Draw(); 
        pad2->cd(); 
        pad2->SetBottomMargin(0.3); 
        pad2->SetGrid(kTRUE); 
        gStyle->SetOptStat(kFALSE);  

        TH1D * fracDiff= diffHist(&fakeDataHist,&FullFit); 
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

        diff->Print("systematicFitExample_combined.png"); 
        delete diff;
    }
}

TH1D*
diffHist(TH1D * h1,TH1D * h2){
    double minBin=h1->GetXaxis()->GetXmin();
    double maxBin=h1->GetXaxis()->GetXmax();
    double numOfBins=h1->GetNbinsX();

    TH1D* rhist = new TH1D("rhist","",numOfBins,minBin,maxBin);
    for(double i=0;i<numOfBins;i++){
        double h1cont=h1->GetBinContent(i);
        double h2cont=h2->GetBinContent(i);
        double weight;
        double error;
        if (h1cont!=0 && h1cont-h2cont!=0) {
            weight= h2cont/h1cont;
            error= weight * sqrt(1/sqrt(h2cont) + 1/sqrt(h1cont));
        } else {
            // How else do you fix this?
            weight= 10000000;
            error = 0;
        }
        rhist->SetBinContent(i,weight);
        rhist->SetBinError(i,error);
    }
    return rhist;
}

void
padPDFs(std::vector<BinnedED>& binnedEDList){
    std::cout<<"Padding Now"<<std::endl;
    for(int i=0;i<binnedEDList.size();i++){
        std::vector<double> bincontent =binnedEDList.at(i).GetBinContents();
        std::vector<double> new_bincontent;
        for(int j =0; j<bincontent.size();j++){
            if(bincontent.at(j)==0)
                new_bincontent.push_back(std::numeric_limits< double >::min());
            else
                new_bincontent.push_back(bincontent[j]);
        }
        binnedEDList[i].SetBinContents(new_bincontent);
    }
}
