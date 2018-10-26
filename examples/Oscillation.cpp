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
  
std::vector<double> d21bestvals;
std::vector<double> s12bestvals;
std::vector<double> s13bestvals;
std::vector<double> Brucenormbestvals;
std::vector<double> Darlingtonnormbestvals;
std::vector<double> Pickeringnormbestvals;

//Want array of Reactor names/core names which can be; thrown into Rat to find distances,
//used to declare variables and objects

std::vector<std::string> Reactors;

void LHFit(){
  
  Reactors.push_back("BRUCE");
  Reactors.push_back("DARLINGTON");
  Reactors.push_back("PICKERING");

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

  /*typedef std::map<std::string, BinnedED> ReacPdfs;
  ReacPdfs[Reactors[0]]; //= BinnedED BruceUnOscPdf("BruceUnOscPdf",axes);
  ReacPdfs[Reactors[1]]; //= BinnedED BruceUnOscPdf("BruceUnOscPdf",axes);
  ReacPdfs[Reactors[2]]; //=BinnedED BruceUnOscPdf("BruceUnOscPdf",axes);
  
  for (ReacPdfs::const_iterator it = ReacPdfs.begin(); it != ReacPdfs.end(); ++i){
    (it->second)(Reactors[0],axes);
    }*/

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

  ROOTMultiPlot* Plot1 = new ROOTMultiPlot;
  //Plot1->SetStacked(true);
  Plot1->AddPdf(dataSetPdf1, "Bruce");
  Plot1->AddPdf(dataSetPdf2, "Darl");
  Plot1->AddPdf(dataSetPdf3, "Pick");
  
  //poisson stat fluctutate input data
  for(int i = 0; i < dataSetPdf.GetNBins(); i++)
    {
      dataSetPdf.SetBinContent(i, r1->Poisson(dataSetPdf.GetBinContent(i)));
      }
  
  BinnedED BruceUnOscPdf("BruceUnOscPdf",axes);
  BinnedED DarlingtonUnOscPdf("DarlingtonUnOscPdf",axes);
  BinnedED PickeringUnOscPdf("PickeringUnOscPdf",axes);
  
  BruceUnOscPdf.SetObservables(0);
  DarlingtonUnOscPdf.SetObservables(0);
  PickeringUnOscPdf.SetObservables(0);
  
  ROOTNtuple UnOscNtp(UnOscfile, UnOscTreeName);
  for(size_t i = 0; i < UnOscNtp.GetNEntries(); i++){
    BruceUnOscPdf.Fill(UnOscNtp.GetEntry(i));
    DarlingtonUnOscPdf.Fill(UnOscNtp.GetEntry(i));
    PickeringUnOscPdf.Fill(UnOscNtp.GetEntry(i));
    }

  BruceUnOscPdf.Normalise();
  DarlingtonUnOscPdf.Normalise();
  PickeringUnOscPdf.Normalise();

  std::vector<BinnedED> mcPdfs;
  mcPdfs.push_back(BruceUnOscPdf);
  mcPdfs.push_back(DarlingtonUnOscPdf);
  mcPdfs.push_back(PickeringUnOscPdf);

  std::cout << "Initialised Pdfs" << std::endl;

  ////////////////////////////////////////////

  NuOsc* Bruceosc_data = new NuOsc("Bruceosc_data"); //Oscillation systematic
  SurvProb* Brucesurvprob = new SurvProb(0.1,0.1,0.1,dist1,"Brucesurvprob"); // Surv Prob function, with intial parameters delm21,ssqqr12, ssqr13 and NB PRECISE BASELINE for reactor pdf
  Brucesurvprob->RenameParameter("delmsqr21_0","d21");
  Brucesurvprob->RenameParameter("sinsqrtheta12_0","s12");
  Brucesurvprob->RenameParameter("sinsqrtheta13_0","s13");
  Bruceosc_data->SetFunction(Brucesurvprob);
  Bruceosc_data->SetAxes(axes);
  Bruceosc_data->SetTransformationObs(dataRep);
  Bruceosc_data->SetDistributionObs(dataRep);
  
  NuOsc* Darlingtonosc_data = new NuOsc("Darlingtonosc_data"); 
  SurvProb* Darlingtonsurvprob = new SurvProb(0.1,0.1,0.1,dist2,"Darlingtonsurvprob");
  Darlingtonsurvprob->RenameParameter("delmsqr21_0","d21");
  Darlingtonsurvprob->RenameParameter("sinsqrtheta12_0","s12");
  Darlingtonsurvprob->RenameParameter("sinsqrtheta13_0","s13");
  Darlingtonosc_data->SetFunction(Darlingtonsurvprob);
  Darlingtonosc_data->SetAxes(axes);
  Darlingtonosc_data->SetTransformationObs(dataRep);
  Darlingtonosc_data->SetDistributionObs(dataRep);
  
  NuOsc* Pickeringosc_data = new NuOsc("Pickeringosc_data"); 
  SurvProb* Pickeringsurvprob = new SurvProb(0.1,0.1,0.1,dist3,"Pickeringsurvprob");
  Pickeringsurvprob->RenameParameter("delmsqr21_0","d21");
  Pickeringsurvprob->RenameParameter("sinsqrtheta12_0","s12");
  Pickeringsurvprob->RenameParameter("sinsqrtheta13_0","s13");
  Pickeringosc_data->SetFunction(Pickeringsurvprob);
  Pickeringosc_data->SetAxes(axes);
  Pickeringosc_data->SetTransformationObs(dataRep);
  Pickeringosc_data->SetDistributionObs(dataRep);
  
  /////////////////////////////////////////////
  // 3. Set Up LH function & fit parameters  //
  /////////////////////////////////////////////
  
  // Setting optimisation limits
  ParameterDict minima;
  minima["BruceUnOscPdf_norm"]= 0;
  minima["DarlingtonUnOscPdf_norm"]= 0;
  minima["PickeringUnOscPdf_norm"]= 0;
  minima["d21" ] = 0.;
  minima["s12" ] = 0.1;
  minima["s13" ] = 0.01;

  ParameterDict maxima;
  maxima["BruceUnOscPdf_norm"]= 100000;
  maxima["DarlingtonUnOscPdf_norm"]= 50000;
  maxima["PickeringUnOscPdf_norm"]= 100000;
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
  initialval["BruceUnOscPdf_norm"]= 40000;
  initialval["DarlingtonUnOscPdf_norm"]= 9000;
  initialval["PickeringUnOscPdf_norm"]= 9000;
  initialval["d21" ]  = 7.4e-5;
  initialval["s12" ]  = 0.3;
  initialval["s13" ]  = 0.02;
  
  ParameterDict initialerr;
  initialerr["BruceUnOscPdf_norm"]= 0.1*initialval["BruceUnOscPdf_norm"];;
  initialerr["DarlingtonUnOscPdf_norm"]= 0.1*initialval["DarlingtonUnOscPdf_norm"];;
  initialerr["PickeringUnOscPdf_norm"]= 0.1*initialval["PickeringUnOscPdf_norm"];;
  initialerr["d21" ] = 0.1*initialval["d21"];
  initialerr["s12" ] = 0.1*initialval["s12"];
  initialerr["s13" ] = 0.1*initialval["s13"];

  int BuffLow  = 5;
  int BuffHigh = 5;
  
  BinnedNLLH lhFunction;
  lhFunction.SetBufferAsOverflow(true);
  lhFunction.SetBuffer(0,BuffLow,BuffHigh);
  lhFunction.SetDataDist(dataSetPdf); // initialise withe the data set

  lhFunction.AddSystematic(Bruceosc_data,"groupA");
  lhFunction.AddSystematic(Darlingtonosc_data,"groupB");
  lhFunction.AddSystematic(Pickeringosc_data,"groupC");

  lhFunction.AddDist(mcPdfs.at(0),std::vector<std::string>(1,"groupA"));
  lhFunction.AddDist(mcPdfs.at(1),std::vector<std::string>(1,"groupB"));
  lhFunction.AddDist(mcPdfs.at(2),std::vector<std::string>(1,"groupC"));

  //lhFunction.SetConstraint("BruceUnOscPdf_norm",41000,5000);
  lhFunction.SetConstraint("DarlingtonUnOscPdf_norm",9450,2000);
  lhFunction.SetConstraint("PickeringUnOscPdf_norm",9500,2000);
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
  Brucenormbestvals.push_back(bestFit.at("BruceUnOscPdf_norm"));
  Darlingtonnormbestvals.push_back(bestFit.at("DarlingtonUnOscPdf_norm"));
  Pickeringnormbestvals.push_back(bestFit.at("PickeringUnOscPdf_norm"));
  d21bestvals.push_back(bestFit.at("d21"));
  s12bestvals.push_back(bestFit.at("s12"));
  s13bestvals.push_back(bestFit.at("s13"));

  BinnedED Result1Holder = mcPdfs.at(0);
  BinnedED Result1;

  NuOsc BruceOscResult("Oscillated");
  BruceOscResult.SetFunction(new SurvProb(bestFit.at("d21"),bestFit.at("s12"),bestFit.at("s13"),dist1));
  BruceOscResult.SetAxes(axes);
  BruceOscResult.SetTransformationObs(dataRep);
  BruceOscResult.SetDistributionObs(dataRep);
  BruceOscResult.Construct();

  Result1 = BruceOscResult(Result1Holder);
  Result1.Scale(bestFit.at("BruceUnOscPdf_norm"));

  BinnedED Result2Holder = mcPdfs.at(1);
  BinnedED Result2;

  NuOsc DarlingtonOscResult("Oscillated");
  DarlingtonOscResult.SetFunction(new SurvProb(bestFit.at("d21"),bestFit.at("s12"),bestFit.at("s13"),dist2));
  DarlingtonOscResult.SetAxes(axes);
  DarlingtonOscResult.SetTransformationObs(dataRep);
  DarlingtonOscResult.SetDistributionObs(dataRep);
  DarlingtonOscResult.Construct();

  Result2 = DarlingtonOscResult(Result2Holder);
  Result2.Scale(bestFit.at("DarlingtonUnOscPdf_norm"));

  BinnedED Result3Holder = mcPdfs.at(2);
  BinnedED Result3;

  NuOsc PickeringOscResult("Oscillated");
  PickeringOscResult.SetFunction(new SurvProb(bestFit.at("d21"),bestFit.at("s12"),bestFit.at("s13"),dist3));
  PickeringOscResult.SetAxes(axes);
  PickeringOscResult.SetTransformationObs(dataRep);
  PickeringOscResult.SetDistributionObs(dataRep);
  PickeringOscResult.Construct();

  Result3 = PickeringOscResult(Result3Holder);
  Result3.Scale(bestFit.at("PickeringUnOscPdf_norm"));
  
  ROOTMultiPlot* Plot = new ROOTMultiPlot;
  //Plot->SetStacked(true);
  Plot->AddPdf(Result1, "Bruce");
  Plot->AddPdf(Result2, "Darlington");
  Plot->AddPdf(Result3, "Pickering");
  Plot->SaveAs("/home/blakei/oxsx/examples/Result.root");
  
  //Plot1->SetStacked(false);
  //Plot1->AddPdf(Result1, "BruceFit");
  //Plot1->AddPdf(Result2, "DarlingtonFit");
  //Plot1->AddPdf(Result3, "PickeringFit");
  Plot1->SaveAs("/home/blakei/oxsx/examples/Test.root");
  TH1D DataHist;
  TH1D FitHist;
  
  DataHist = DistTools::ToTH1D(dataSetPdf);
  Result1.Add(Result2,1);
  Result1.Add(Result3,1);
  FitHist = DistTools::ToTH1D(Result1);
  
  TH1D FullFit("FullFit","",FitHist.GetNbinsX(),FitHist.GetXaxis()->GetXmin(),FitHist.GetXaxis()->GetXmax());
  
  FullFit.Add(&FitHist);
  
  DataHist.Sumw2();

  TLegend* leg = new TLegend(0.75,0.8,1.0,1.0);
  leg->AddEntry(&DataHist,"Data","lf");
  leg->AddEntry(&FitHist,"Fit Result","lf");

  TCanvas* c1 = new TCanvas("c1");
  c1->cd();
  
  gStyle->SetOptStat(0);
  
  DataHist.SetTitle("Data to Fit");  
  DataHist.GetYaxis()->SetTitle(Form("Counts"));  
  
  DataHist.Draw();
  FitHist.SetLineColor(kRed);
  FitHist.SetLineWidth(3);
  FitHist.Draw("same e");

  leg->Draw();
  
  TPaveText pt(0.75,0.35,1.0,0.65,"NDC");
  pt.AddText(Form("norm = %.5f" ,bestFit["BruceUnOscPdf_norm"]));
  pt.AddText(Form("norm = %.5f" ,bestFit["DarlingtonUnOscPdf_norm"]));
  pt.AddText(Form("norm = %.5f" ,bestFit["PickeringUnOscPdf_norm"]));
  pt.AddText(Form("#Delta m_{21} = %.6f",bestFit["d21"]));
  pt.AddText(Form("#theta_{12} = %.3f",bestFit["s12"]));
  pt.AddText(Form("#theta_{13} = %.4f",bestFit["s13"]));
  pt.SetFillColor(kWhite);
  pt.SetShadowColor(kWhite);
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
  TH1D* fits12vals = new TH1D("fits12vals","fits12vals",150 ,0.15 ,0.3);

  for (int i = 0; i < numexps; i++){
    LHFit();
    //double Rand = r1->Poisson(50);
  //std::cout<<Rand<<std::endl;
  }
  for (int i = 0; i < numexps; i++){
    fits12vals->Fill(s12bestvals[i]);
  }

  TFile* LHfits = new TFile("/home/blakei/oxsx/examples/LHfits1.root","RECREATE");
  fits12vals->Write();
  
  LHfits->Close();
  return 0;
}
