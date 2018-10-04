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

#include <TRandom3.h>

// data (ntuples) to load
//const std::string bgMCfile    = "/data/snoplus/blakei/antinu/mc/ntuples/Oscbg.root";
//const std::string bgTreeName  = "nt";
const std::string UnOscfile   = "/data/snoplus/blakei/antinu/mc/ntuples/test/UnOscDarlingtonflux1000_oxsx.root";
const std::string UnOscTreeName = "nt";

const std::string dataFile = "/data/snoplus/blakei/antinu/mc/ntuples/test/OscDarlingtonflux1000ds21_7.4e-05_ss12_0.297_ss13_0.0215_oxsx.root";
const std::string dataTreeName = "nt";

double dist = 349;

int numexps = 1; 

double Emin = 2;
double Emax = 8;
int numbins = 60;
  
std::vector<double> d21bestvals;
std::vector<double> s12bestvals;
std::vector<double> s13bestvals;
std::vector<double> normbestvals;

void LHFit(){
  
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
  dataSetPdf.SetObservables(dataRep);
  ROOTNtuple dataNtp(dataFile, dataTreeName);
  for(size_t i = 0; i < dataNtp.GetNEntries(); i++)
    dataSetPdf.Fill(dataNtp.GetEntry(i));

  //poisson stat fluctutate input data
  /*for(int i = 0; i < dataSetPdf.GetNBins(); i++)
    {
      dataSetPdf.SetBinContent(i, r1->Poisson(dataSetPdf.GetBinContent(i)));
      }*/
  BinnedED UnOscPdf("UnOscPdf",axes);
  UnOscPdf.SetObservables(0);
  ROOTNtuple UnOscNtp(UnOscfile, UnOscTreeName);
  for(size_t i = 0; i < UnOscNtp.GetNEntries(); i++)
    UnOscPdf.Fill(UnOscNtp.GetEntry(i));
  
  UnOscPdf.Normalise();

  std::vector<BinnedED> mcPdfs;
  mcPdfs.push_back(UnOscPdf);

  std::cout << "Initialised Pdfs" << std::endl;

  ////////////////////////////////////////////

  NuOsc* osc_data = new NuOsc("osc_data"); //Oscillation systematic
  SurvProb* survprob = new SurvProb(0.1,0.1,0.1,dist,"survprob"); // Surv Prob function, with intial parameters 
  //delm21,ssqqr12, ssqr13 (init. w/ any value >0) and NB PRECISE BASELINE for reactor pdf
  survprob->RenameParameter("delmsqr21_0","d21");
  survprob->RenameParameter("sinsqrtheta12_0","s12");
  survprob->RenameParameter("sinsqrtheta13_0","s13");

  osc_data->SetFunction(survprob);
  osc_data->SetAxes(axes);
  osc_data->SetTransformationObs(dataRep);
  osc_data->SetDistributionObs(dataRep);
  //std::cout << "constructing.." << std::endl;
  //osc_data->Construct();
  //std::cout << "constructed" << std::endl;

  //SparseMatrix sp = osc_data->GetResponse();  //not essential
  //sp.PrintMat();                              //not essential
  //TH1D hist("hist","",numbins,Emin,Emax);
  //for (int i = 0; i < sp.GetNRows(); ++i) {   //not essential
  //hist.Fill(2+(i*(10-2)/80.),sp.GetComponent(i,i));
  //}

  /////////////////////////////////////////////
  // 3. Set Up LH function & fit parameters  //
  /////////////////////////////////////////////
  
  // Setting optimisation limits
  ParameterDict minima;
  minima["UnOscPdf_norm"]= 0;
  minima["d21" ] = 0.;
  minima["s12" ] = 0.1;
  minima["s13" ] = 0.01;

  ParameterDict maxima;
  maxima["UnOscPdf_norm"]= 100000;
  maxima["d21" ] = 0.0001;
  maxima["s12" ] = 0.5;
  maxima["s13" ] = 0.05;
  /*
  ParameterDict initialval;
  initialval["UnOscPdf_norm"] = rand.UniformRange(minima["UnOscPdf_norm"],maxima["UnOscPdf_norm"]);
  initialval["d21" ] = rand.UniformRange(minima["d21" ],maxima["d21" ]);
  initialval["s12" ] = rand.UniformRange(minima["s12" ],maxima["s12" ]);
  initialval["s13" ] = rand.UniformRange(minima["s13" ],maxima["s13" ]);
  */
  //TRandom3 *r1 = new TRandom3();
  //double Rand = r1->Rndm();
  //std::cout<<"intial s12 value: "<<minima["s12" ]+ (Rand*(maxima["s12" ]-minima["s12" ]))<<std::endl;

  ParameterDict initialval;
  initialval["UnOscPdf_norm"]= 9000;
  initialval["d21" ]  = 7.4e-5;
  initialval["s12" ]  = 0.3;
  initialval["s13" ]  = 0.02;
  
  ParameterDict initialerr;
  initialerr["UnOscPdf_norm"]= 0.1*initialval["UnOscPdf_norm"];;
  initialerr["d21" ] = 0.1*initialval["d21"];
  initialerr["s12" ] = 0.1*initialval["s12"];
  initialerr["s13" ] = 0.1*initialval["s13"];

  int BuffLow  = 5;
  int BuffHigh = 5;
  
  BinnedNLLH lhFunction;
  lhFunction.SetDataDist(dataSetPdf); // initialise withe the data set
  lhFunction.SetBufferAsOverflow(true);
  lhFunction.SetBuffer(0,BuffLow,BuffHigh);
  lhFunction.AddSystematic(osc_data,"groupA");
  lhFunction.AddDist(mcPdfs.at(0),std::vector<std::string>(1,"groupA"));
  //lhFunction.SetConstraint("UnOscPdf_norm",6000,1000);
  //lhFunction.SetConstraint("d21",7.37e-5,1.6e-6);
  lhFunction.SetConstraint("s12",0.297,0.016);
  lhFunction.SetConstraint("s13",0.0215,0.009);
  
  std::cout << "Built LH function " << std::endl;

  ////////////
  // 4. Fit //
  ////////////
  Minuit min;
  min.SetMethod("Migrad");
  min.SetMaxCalls(1000000);
  min.SetMinima(minima);
  min.SetMaxima(maxima);
  min.SetInitialValues(initialval);
  min.SetInitialErrors(initialerr);

  FitResult fitResult = min.Optimise(&lhFunction);
  ParameterDict bestFit = fitResult.GetBestFit();
  fitResult.Print();

  /////////////////////////////////////////////

  d21bestvals.push_back(bestFit.at("d21"));
  s12bestvals.push_back(bestFit.at("s12"));
  s13bestvals.push_back(bestFit.at("s13"));
  normbestvals.push_back(bestFit.at("UnOscPdf_norm"));

  BinnedED ResultHolder = mcPdfs.at(0);
  BinnedED Result;

  NuOsc OscResult("Oscillated");
  OscResult.SetFunction(new SurvProb(bestFit.at("d21"),bestFit.at("s12"),bestFit.at("s13"),dist));
  OscResult.SetAxes(axes);
  OscResult.SetTransformationObs(dataRep);
  OscResult.SetDistributionObs(dataRep);
  OscResult.Construct();

  Result = OscResult(ResultHolder);
  Result.Scale(bestFit.at("UnOscPdf_norm"));

  TH1D DataHist;
  TH1D FitHist;
  
  //DataHist.SetStats(kFALSE);
  //FitHist.SetStats(kFALSE);
  DataHist = DistTools::ToTH1D(dataSetPdf);
  FitHist = DistTools::ToTH1D(Result);
  
  TH1D FullFit("FullFit","",FitHist.GetNbinsX(),FitHist.GetXaxis()->GetXmin(),FitHist.GetXaxis()->GetXmax());
  
  FullFit.Add(&FitHist);
  //FullFit.Sumw2();
  
  DataHist.Sumw2();

  TLegend* leg = new TLegend(0.75,0.8,1.0,1.0);
  leg->AddEntry(&DataHist,"Data","lf");
  leg->AddEntry(&FitHist,"Fit Result","lf");

  //TCanvas* c1 = new TCanvas("c1","",800,800);
  TCanvas* c1 = new TCanvas("c1");
  c1->cd();
  
  gStyle->SetOptStat(0);
  //TPad *pad1 = new TPad("pad1","pad1", 0,0.3,1,1);
  //pad1->Draw();
  //pad1->cd();
  //pad1->SetBottomMargin(0.0);
  //gPad->RedrawAxis();

  DataHist.SetTitle("Data to Fit");  
  DataHist.GetYaxis()->SetTitle(Form("Counts"));  
  //DataHist.GetYaxis()->SetTitleOffset(1.5);
  //DataHist.SetFillColorAlpha(kGreen,0.5);
  
  DataHist.Draw();
  //FullFit.SetFillColorAlpha(kRed,0.5);
  FitHist.SetLineColor(kRed);
  FitHist.SetLineWidth(3);
  FitHist.Draw("same e");

  leg->Draw();
  
  TPaveText pt(0.75,0.35,1.0,0.65,"NDC");
  pt.AddText(Form("norm = %.5f" ,bestFit["UnOscPdf_norm"]));
  pt.AddText(Form("#Delta m_{21} = %.6f",bestFit["d21"]));
  pt.AddText(Form("#theta_{12} = %.3f",bestFit["s12"]));
  pt.AddText(Form("#theta_{13} = %.4f",bestFit["s13"]));
  pt.SetFillColor(kWhite);
  pt.SetShadowColor(kWhite);
  pt.Draw();
  c1->cd();
  
  TFile * fitout = new TFile("/home/blakei/oxsx/examples/FitOut.root","RECREATE");
  c1->Write();
  FullFit.Write();
  DataHist.Write();
  
  fitout->Close();
  
  return;
}


int main(){
  //TRandom3 *r1 = new TRandom3();
  //r1->SetSeed(0);
  TH1D* fits12vals = new TH1D("fits12vals","fits12vals",150 ,0.15 ,0.3);

  for (int i = 0; i < numexps; i++){
    LHFit();
    //double Rand = r1->Poisson(50);
  //std::cout<<Rand<<std::endl;
  }
  for (int i = 0; i < numexps; i++){
    fits12vals->Fill(s12bestvals[i]);
  }

  TFile* LHfits = new TFile("/home/blakei/oxsx/examples/LHfits.root","RECREATE");
  fits12vals->Write();
  
  LHfits->Close();
  return 0;
}
