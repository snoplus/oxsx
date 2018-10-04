// A simple fit in energy for signal and a background
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
#include <ROOTNtuple.h>

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
#include <ROOTMultiPlot.h>
#include <TRandom3.h>

// data (ntuples) to load
//const std::string bgMCfile    = "/data/snoplus/blakei/antinu/mc/ntuples/Oscbg.root";
//const std::string bgTreeName  = "nt";
const std::string UnOscfile   = "/data/snoplus/blakei/antinu/mc/ntuples/test/UnOscBruceflux1000_oxsx.root";
const std::string UnOscTreeName = "nt";

const std::string dataFile1 = "/data/snoplus/blakei/antinu/mc/ntuples/test/OscBruceflux1000ds21_7.4e-05_ss12_0.297_ss13_0.0215_oxsx.root";
const std::string dataTreeName1 = "nt";
const std::string dataFile2 = "/data/snoplus/blakei/antinu/mc/ntuples/test/OscDarlingtonflux1000ds21_7.4e-05_ss12_0.297_ss13_0.0215_oxsx.root";
const std::string dataTreeName2 = "nt";
const std::string dataFile3 = "/data/snoplus/blakei/antinu/mc/ntuples/test/OscPickeringflux1000ds21_7.4e-05_ss12_0.297_ss13_0.0215_oxsx.root";
const std::string dataTreeName3 = "nt";

double dist1 = 240.22;
double dist2 = 349.147;
double dist3 = 340.37;

int numexps = 1; 

double Emin = 2;
double Emax = 8;
int numbins = 60;

int BuffLow  = 5;
int BuffHigh = 5;
  
std::vector<double> d21bestvals;
std::vector<double> s12bestvals;
std::vector<double> s13bestvals;
std::vector<double> Brucenormbestvals;
std::vector<double> Darlingtonnormbestvals;
std::vector<double> Pickeringnormbestvals;

std::vector<double> reactorDistances;

//Want array of Reactor names/core names which can be; thrown into Rat to find distances,
//used to declare variables and objects

std::vector<std::string> Reactors;

void LHFit(){
  
  Reactors.push_back("BRUCE");
  Reactors.push_back("DARLINGTON");
  Reactors.push_back("PICKERING");

  reactorDistances.push_back(dist1);
  reactorDistances.push_back(dist2);
  reactorDistances.push_back(dist3);
  
  int numPdfs = Reactors.size();
  gStyle->SetOptStat(0);
  
  TRandom3 *r1 = new TRandom3();
  r1->SetSeed(0);
  
  ////////////////////
  // 1. Set Up PDFs //
  ////////////////////
  
  // Only interested in first bit of data ntuple
  ObsSet dataRep(0);

  // Set up binning
  AxisCollection axes;
  axes.AddAxis(BinAxis("ParKE", Emin, Emax, numbins));
  
  BinnedED  dataSetPdf("dataSetPdf",axes);
  BinnedED  dataSetPdf1("dataSetPdf1",axes);
  BinnedED  dataSetPdf2("dataSetPdf2",axes);
  BinnedED  dataSetPdf3("dataSetPdf3",axes);

  dataSetPdf.SetObservables(dataRep);
  dataSetPdf1.SetObservables(dataRep);
  dataSetPdf2.SetObservables(dataRep);
  dataSetPdf3.SetObservables(dataRep);
  ROOTNtuple dataNtp1(dataFile1, dataTreeName1);
  ROOTNtuple dataNtp2(dataFile2, dataTreeName2);
  ROOTNtuple dataNtp3(dataFile3, dataTreeName3);
  for(size_t i = 0; i < dataNtp1.GetNEntries(); i++)
    dataSetPdf1.Fill(dataNtp1.GetEntry(i));
  for(size_t i = 0; i < dataNtp2.GetNEntries(); i++)
    dataSetPdf2.Fill(dataNtp2.GetEntry(i));
  for(size_t i = 0; i < dataNtp3.GetNEntries(); i++)
    dataSetPdf3.Fill(dataNtp3.GetEntry(i));

  dataSetPdf1.Normalise();
  dataSetPdf2.Normalise();
  dataSetPdf3.Normalise();
  
  dataSetPdf1.Scale(41500);
  dataSetPdf2.Scale(9450);
  dataSetPdf3.Scale(9500);
  
  dataSetPdf.Add(dataSetPdf1,1);
  dataSetPdf.Add(dataSetPdf2,1);
  dataSetPdf.Add(dataSetPdf3,1);
  //poisson stat fluctutate input data
  for(int i = 0; i < dataSetPdf.GetNBins(); i++)
    {
      dataSetPdf.SetBinContent(i, r1->Poisson(dataSetPdf.GetBinContent(i)));
    }
  char name[100];
  //BinnedED *reactorPdf[numPdfs];
  std::vector<BinnedED*> reactorPdf(numPdfs);
  ROOTNtuple reactorNtp(UnOscfile, UnOscTreeName);
  NuOsc *reactorSystematic[numPdfs];
  ParameterDict minima;
  ParameterDict maxima;
  ParameterDict initialval;
  ParameterDict initialerr;

  BinnedNLLH lhFunction;
  lhFunction.SetBufferAsOverflow(true);
  lhFunction.SetBuffer(0,BuffLow,BuffHigh);
  lhFunction.SetDataDist(dataSetPdf); // initialise withe the data set
  
  minima["d21"] = 0.;
  minima["s12"] = 0.1;
  minima["s13"] = 0.01;
  maxima["d21" ] = 0.0001;
  maxima["s12" ] = 0.5;
  maxima["s13" ] = 0.05;
  initialval["d21" ]  = 7.4e-5;
  initialval["s12" ]  = 0.3;
  initialval["s13" ]  = 0.02;
  initialerr["d21" ] = 0.1*initialval["d21"];
  initialerr["s12" ] = 0.1*initialval["s12"];
  initialerr["s13" ] = 0.1*initialval["s13"];

  for (int i = 0; i< numPdfs; i++){
    sprintf(name,"ReactorPdf%d",i);
    reactorPdf[i] = new BinnedED(name,axes);
    reactorPdf[i]->SetObservables(0);
    for(size_t j = 0; j < reactorNtp.GetNEntries(); j++){
      reactorPdf[i]->Fill(reactorNtp.GetEntry(j));
    }
    reactorPdf[i]->Normalise();
    
    // create and fill the oscillated systematics
    sprintf(name,"ReactorSystematic%d",i);
    reactorSystematic[i] = new NuOsc(name);
    
    sprintf(name,"ReactorSurvival%d",i);
    SurvProb* survprob = new SurvProb(0.1,0.1,0.1,reactorDistances[i],name); // Surv Prob function, with intial parameters delm21,ssqqr12, ssqr13 and NB PRECISE BASELINE for reactor pdf
    survprob->RenameParameter("delmsqr21_0","d21");
    survprob->RenameParameter("sinsqrtheta12_0","s12");
    survprob->RenameParameter("sinsqrtheta13_0","s13");
    reactorSystematic[i]->SetFunction(survprob);
    reactorSystematic[i]->SetAxes(axes);
    reactorSystematic[i]->SetTransformationObs(dataRep);
    reactorSystematic[i]->SetDistributionObs(dataRep);
    
    // Setting optimisation limits
    sprintf(name,"ReactorPdf%d_norm",i);
    minima[name] = 0;
    maxima[name] = 100000;
    initialval[name] = 50000;
    initialerr[name] = 0.1*initialval[name];

    sprintf(name,"group%d",i);
    lhFunction.AddSystematic(reactorSystematic[i],name);
    lhFunction.AddDist(*reactorPdf[i],std::vector<std::string>(1,name));
  }

  std::cout << "Initialised Pdfs" << std::endl;

  //lhFunction.SetConstraint("ReactorPdf0_norm",41000,5000);
  lhFunction.SetConstraint("ReactorPdf1_norm",9450,2000);
  lhFunction.SetConstraint("ReactorPdf2_norm",9500,2000);
  
  //lhFunction.SetConstraint("d21",7.37e-5,1.6e-6);
  lhFunction.SetConstraint("s12",0.297,0.016);
  lhFunction.SetConstraint("s13",0.0215,0.009);
  
  std::cout << "Built LH function " << std::endl;

  ////////////
  // 4. Fit //
  ////////////
  Minuit min;
  min.SetMethod("Migrad");
  min.SetMaxCalls(10000000);
  min.SetMinima(minima);
  min.SetMaxima(maxima);
  min.SetInitialValues(initialval);
  min.SetInitialErrors(initialerr);

  FitResult fitResult = min.Optimise(&lhFunction);
  ParameterDict bestFit = fitResult.GetBestFit();
  fitResult.Print();
  
  /////////////////////////////////////////////
  ////////        Fit Result        ///////////
  /////////////////////////////////////////////
  
  BinnedED Result("Result",axes);
  BinnedED TotalResult("TotalResult",axes);
  ROOTMultiPlot* Plot = new ROOTMultiPlot;

  TPaveText pt(0.75,0.35,1.0,0.65,"NDC");
  for (int i = 0; i< numPdfs; i++){
    Result.Empty();
    NuOsc OscResult("OscResult");
    OscResult.SetFunction(new SurvProb(bestFit.at("d21"),bestFit.at("s12"),bestFit.at("s13"),reactorDistances[i]));
    OscResult.SetAxes(axes);
    OscResult.SetTransformationObs(dataRep);
    OscResult.SetDistributionObs(dataRep);
    OscResult.Construct();

    Result.Add(OscResult(*reactorPdf[i]),1);
    sprintf(name,"ReactorPdf%d_norm",i);
    Result.Scale(bestFit.at(name));

    Plot->AddPdf(Result, name);
    TotalResult.Add(Result,1);
    
    pt.AddText(Form("norm = %.5f" ,bestFit[name]));
  }
  Plot->SaveAs("/home/blakei/oxsx/examples/Result.root");
 
  pt.AddText(Form("#Delta m_{21} = %.6f",bestFit["d21"]));
  pt.AddText(Form("#theta_{12} = %.3f",bestFit["s12"]));
  pt.AddText(Form("#theta_{13} = %.4f",bestFit["s13"]));
  pt.SetFillColor(kWhite);
  pt.SetShadowColor(kWhite);
 
  //Brucenormbestvals.push_back(bestFit.at("BruceUnOscPdf_norm"));
  //Darlingtonnormbestvals.push_back(bestFit.at("DarlingtonUnOscPdf_norm"));
  //Pickeringnormbestvals.push_back(bestFit.at("PickeringUnOscPdf_norm"));
  //d21bestvals.push_back(bestFit.at("d21"));
  //s12bestvals.push_back(bestFit.at("s12"));
  //s13bestvals.push_back(bestFit.at("s13"));

  TH1D DataHist;
  TH1D FitHist;
  
  DataHist = DistTools::ToTH1D(dataSetPdf);
  FitHist = DistTools::ToTH1D(TotalResult);
  
  TH1D FullFit("FullFit","",FitHist.GetNbinsX(),FitHist.GetXaxis()->GetXmin(),FitHist.GetXaxis()->GetXmax());
  FullFit.Add(&FitHist);
  DataHist.Sumw2();

  TLegend* leg = new TLegend(0.75,0.8,1.0,1.0);
  leg->AddEntry(&DataHist,"Data","lf");
  leg->AddEntry(&FitHist,"Fit Result","lf");

  TCanvas* c1 = new TCanvas("c1");
  c1->cd();  
  DataHist.SetTitle("Data to Fit");  
  DataHist.GetYaxis()->SetTitle(Form("Counts"));
  DataHist.Draw();
  FitHist.SetLineColor(kRed);
  FitHist.SetLineWidth(3);
  FitHist.Draw("same e");
  leg->Draw();
  
  pt.Draw();
  c1->cd();
  
  TFile * fitout = new TFile("/home/blakei/oxsx/examples/FitOut1.root","RECREATE");
  c1->Write();
  FullFit.Write();
  DataHist.Write();
  
  fitout->Close();
  
  return;
}


int main(){
  //TRandom3 *r1 = new TRandom3();
  //r1->SetSeed(0);
  //TH1D* fits12vals = new TH1D("fits12vals","fits12vals",150 ,0.15 ,0.3);

  for (int i = 0; i < numexps; i++){
    LHFit();
  }
  //for (int i = 0; i < numexps; i++){
  //fits12vals->Fill(s12bestvals[i]);
  //}

  //TFile* LHfits = new TFile("/home/blakei/oxsx/examples/LHfits1.root","RECREATE");
  //fits12vals->Write();
  
  //LHfits->Close();
  return 0;
}
