// #include <iostream>
// #include <NuOsc.h>
//
// int main(int argc, char *argv[])
// {
//   std::cout << "jfk" << std::endl;
//   NuOsc nuOsc("nuOsc");
//   return 0;
// }

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
#include <NuOsc.h>
#include <SurvProb.h>
#include <TH1D.h>

TH1D* diffHist(TH1D * h1,TH1D * h2);
void padPDFs(std::vector<BinnedED>& binnedEDList);

const std::string UnOscMCfile   = "/data/snoplus/blakei/antinu/mc/ntuples/unoscBruce_oxsx.root";
const std::string UnOscTreeName  = "nt";

const std::string dataFile = "/data/snoplus/blakei/antinu/mc/ntuples/Oscdata.root";
const std::string dataTreeName = "nt";

void
function(){

    Rand::SetSeed(0);
    ObsSet  obsSet(0); // Only interested in first bit of data ntuple

    Gaussian gaus1(6, 1);

    AxisCollection axes;
    axes.AddAxis(BinAxis("energy", 2, 10 ,64));


    NuOsc* data_osc = new NuOsc("data_osc"); //Oscillation systematic
    SurvProb data(7e-5,0.3,0.02,238,"data"); // Surv Prob function, with intial parameters delm21,ssqqr12, ssqr13 and NB PRECISE BASELINE for reactor pdf

    data_osc->SetFunction(&data);
    data_osc->SetAxes(axes);
    data_osc->SetTransformationObs(obsSet);
    data_osc->SetDistributionObs(obsSet);
    data_osc->Construct();
    SparseMatrix sp = data_osc->GetResponse();  //not essential
    sp.PrintMat();                              //not essential
    TH1D hist("hist","",80,2,10);
    for (int i = 0; i < sp.GetNRows(); ++i) {   //not essential
      hist.Fill(2+(i*(10-2)/80.),sp.GetComponent(i,i));
    }

    TCanvas c1;
    hist.Draw();
    c1.Print("iwann.png");

    BinnedED pdf1("a_mc", DistTools::ToHist(gaus1, axes));
    BinnedED pdf3 = data_osc->operator()(pdf1);
    pdf1.SetObservables(0);
    pdf3.SetObservables(0);

    pdf1.Scale(100000);

    std::vector<BinnedED> dataPdfs;
    dataPdfs.push_back(pdf3);
    padPDFs(dataPdfs);

    BinnedEDGenerator dataGen;
    dataGen.SetPdfs(dataPdfs);
    std::vector<double> rates(1,100000);
    dataGen.SetRates(rates);

    // BinnedED fakeData= dataGen.ExpectedRatesED();
    BinnedED fakeData= dataGen.PoissonFluctuatedED();

    pdf1.Normalise();
    std::vector<BinnedED> mcPdfs;
    mcPdfs.push_back(pdf1);

    padPDFs(mcPdfs);


    NuOsc* conv_a = new NuOsc("conv_a");

    // Gaussian* gaus_a = new Gaussian(0,1,"gaus_a");
    SurvProb* gaus_a = new SurvProb(0.1,1,0.5,240,"gaus_a");
    gaus_a->RenameParameter("delmsqr21_0","d21");
    gaus_a->RenameParameter("sinsqrtheta12_0","s12");
    gaus_a->RenameParameter("sinsqrtheta13_0","s13");

    conv_a->SetFunction(gaus_a);

    conv_a->SetAxes(axes);
    conv_a->SetTransformationObs(obsSet);
    conv_a->SetDistributionObs(obsSet);
    std::cout << "iwan" << std::endl;
    conv_a->Construct();
    std::cout << "iwan2" << std::endl;

    // Setting optimisation limits
    ParameterDict minima;
    minima["a_mc_norm"] = 10;
    minima["d21" ] = 0.;
    minima["s12" ] = 0;
    minima["s13" ] = 0;

    ParameterDict maxima;
    maxima["a_mc_norm"] = 200000;
    maxima["d21" ] = 0.0001;
    maxima["s12" ] = 1;
    maxima["s13" ] = 1;

    // ParameterDict initialval;
    // Rand rand;
    // initialval["a_mc_norm"] = rand.UniformRange(minima["a_mc_norm"],maxima["a_mc_norm"]);
    // initialval["gaus_a_1" ] = rand.UniformRange(minima["gaus_a_1" ],maxima["gaus_a_1" ]);
    // initialval["gaus_a_2" ] = rand.UniformRange(minima["gaus_a_2" ],maxima["gaus_a_2" ]);

    ParameterDict initialval;
    initialval["a_mc_norm"] = 90000;
    initialval["d21" ]  = 0.00008;
    initialval["s12" ]  = 0.9;
    initialval["s13" ]  = 0.5;

    ParameterDict initialerr;
    initialerr["a_mc_norm"] = 0.1*initialval["a_mc_norm"];
    initialerr["d21" ] = 0.1*initialval["d21"];
    initialerr["s12" ] = 0.1*initialval["s12"];
    initialerr["s13" ] = 0.1*initialval["s13"];

    //Setting up likelihood.
    int BuffLow  = 3;
    int BuffHigh = 3;

    BinnedNLLH lh;
    lh.SetBufferAsOverflow(false);
    lh.SetBuffer(0,BuffLow,BuffHigh);
    lh.SetDataDist(fakeData); // initialise with the data set
    //Add systematics to groups.
    lh.AddSystematic(conv_a);

    // Associate EDs with a vector of groups.
    lh.AddDist(mcPdfs.at(0));

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

    // Plot Result
    {
        BinnedED BiHolder = mcPdfs.at(0);
        BinnedED BiResult;

        NuOsc BiSmearer("aConv");
        BiSmearer.SetFunction(new SurvProb(bestResult.at("d21"),bestResult.at("s12"),bestResult.at("s13"),240));
        BiSmearer.SetAxes(axes);
        BiSmearer.SetTransformationObs(obsSet);
        BiSmearer.SetDistributionObs(obsSet);
        BiSmearer.Construct();


        BiResult = BiSmearer( BiHolder );

        BiResult.Scale(bestResult.at("a_mc_norm"));

        TH1D fakeDataHist;
        TH1D BiFit;
        BiFit.Sumw2();
        fakeDataHist = DistTools::ToTH1D(fakeData);
        BiFit= DistTools::ToTH1D(BiResult);

        TH1D FullFit("FullFit","",BiFit.GetNbinsX(),BiFit.GetXaxis()->GetXmin(),BiFit.GetXaxis()->GetXmax());
        FullFit.Sumw2();

        FullFit.Add(&BiFit);

        TLegend* leg =new TLegend(0.8,0.7,1,0.9);
        leg->AddEntry(&fakeDataHist,"FakeData","lf");
        leg->AddEntry(&BiFit,"a fit result","lf");

        TCanvas* diff = new TCanvas("diff","",800,800);
        diff->cd();
        // -------------- Top panel
        gStyle->SetOptStat(kFALSE);
        TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
        pad1->Draw();
        //pad1->SetLogy(1);
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

        TH1D hist1 =  DistTools::ToTH1D(pdf1);
        hist1.Scale(4000);

        hist1.SetFillColorAlpha(kRed,0.5);
        hist1.SetLineWidth(2);
        hist1.Draw("same");

        leg->AddEntry(&hist1,"a pdf","lf");
        leg->Draw();

        TPaveText pt(0.7,0.2,1.0,0.6,"NDC");
        pt.AddText(Form("a norm = %.2f" ,bestResult["a_mc_norm"]));
        pt.AddText(Form("#Delta m_{21} = %.2f",bestResult["d21"]));
        pt.AddText(Form("#theta_{12} = %.2f",bestResult["s12"]));
        pt.AddText(Form("#theta_{13} = %.2f",bestResult["s13"]));
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

        diff->Print("simpleSystematicFitExample.png");
    }

}

int main(int argc, char *argv[]){

    function();
    return 0;

}

TH1D* diffHist(TH1D * h1,TH1D * h2){
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
            // weight= (h1cont-h2cont)/h1cont;
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

void padPDFs(std::vector<BinnedED>& binnedEDList){
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
