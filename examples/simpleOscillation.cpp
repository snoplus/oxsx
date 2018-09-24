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


// data (ntuples) to load
//const std::string bgMCfile    = "/data/snoplus/blakei/antinu/mc/ntuples/Oscbg.root";
const std::string UnOscfile   = "/data/snoplus/blakei/antinu/mc/ntuples/unoscBruce_oxsx.root";
//const std::string bgTreeName  = "nt";
const std::string UnOscTreeName = "nt";

const std::string dataFile = "/data/snoplus/blakei/antinu/mc/ntuples/Oscdata.root";
const std::string dataTreeName = "nt";

int main(){
  ////////////////////
  // 1. Set Up PDFs //
  ////////////////////

  // Set up binning
  AxisCollection axes;
  axes.AddAxis(BinAxis("ParKE", 2, 8, 50));

  // Only interested in first bit of data ntuple
  ObsSet dataRep(0);
  
  BinnedED  dataSetPdf("dataSetPdf",axes);
  dataSetPdf.SetObservables(dataRep);
  ROOTNtuple dataNtp(dataFile, dataTreeName);
  for(size_t i = 0; i < dataNtp.GetNEntries(); i++)
    dataSetPdf.Fill(dataNtp.GetEntry(i));

  BinnedED UnOscPdf("UnOscPdf",axes);
  UnOscPdf.SetObservables(0);
  ROOTNtuple UnOscNtp(UnOscfile, UnOscTreeName);
  for(size_t i = 0; i < UnOscNtp.GetNEntries(); i++)
    UnOscPdf.Fill(UnOscNtp.GetEntry(i));

  std::cout << "Initialised Pdfs" << std::endl;

  //BinnedED pdf1("UnOsc", DistTools::ToHist(gaus1, axes));
  //BinnedED oscillatedpdf = survprob->operator()(UnOscPdf);
  //oscillatedpdf.SetObservables(0);
  
  // Set up pdf with these bins in this observable
  //BinnedED bgPdf("bgPdf",axes);
  //bgPdf.SetObservables(dataRep);
  //BinnedED  OscPdf("UnOscPdf",axes);
  //signalPdf.SetObservables(dataRep);

  ///////////////////////
  // 2. Fill with data //
  ///////////////////////
  //ROOTNtuple bgMCNtp(bgMCfile, bgTreeName);
  
  //for(size_t i = 0; i < bgMCNtp.GetNEntries(); i++)
  //bgPdf.Fill(bgMCNtp.GetEntry(i));
  //for(size_t i = 0; i < signalMCNtp.GetNEntries(); i++)
  //signalPdf.Fill(signalMCNtp.GetEntry(i));
  
  //bgPdf.Normalise();
  //signalPdf.Normalise();
  UnOscPdf.Normalise();
  
  std::vector<BinnedED> mcPdfs;
  mcPdfs.push_back(UnOscPdf);

  std::cout << "Filled pdfs " << std::endl;

  ////////////////////////////////////////////

  NuOsc* osc_data = new NuOsc("osc_data"); //Oscillation systematic
  SurvProb* survprob = new SurvProb(7e-5,0.3,0.02,238,"survprob"); // Surv Prob function, with intial parameters delm21,ssqqr12, ssqr13 and NB PRECISE BASELINE for reactor pdf

  survprob->RenameParameter("delmsqr21_0","d21");
  survprob->RenameParameter("sinsqrtheta12_0","s12");
  survprob->RenameParameter("sinsqrtheta13_0","s13");

  osc_data->SetFunction(survprob);

  osc_data->SetAxes(axes);
  osc_data->SetTransformationObs(dataRep);
  osc_data->SetDistributionObs(dataRep);
  std::cout << "constructing.." << std::endl;
  osc_data->Construct();
  std::cout << "constructed" << std::endl;

  /////////////////////////////////////////////
  // 3. Set Up LH function & fit parameters  //
  /////////////////////////////////////////////
  
  // Setting optimisation limits
  ParameterDict minima;
  minima["UnOscPdf_norm"]= 0;
  minima["d21" ] = 0.;
  minima["s12" ] = 0.05;
  minima["s13" ] = 0.001;

  ParameterDict maxima;
  maxima["UnOscPdf_norm"]= 100000;
  maxima["d21" ] = 0.0001;
  maxima["s12" ] = 0.5;
  maxima["s13" ] = 0.5;

  // ParameterDict initialval;
  // Rand rand;
  // initialval["a_mc_norm"] = rand.UniformRange(minima["a_mc_norm"],maxima["a_mc_norm"]);
  // initialval["gaus_a_1" ] = rand.UniformRange(minima["gaus_a_1" ],maxima["gaus_a_1" ]);
  // initialval["gaus_a_2" ] = rand.UniformRange(minima["gaus_a_2" ],maxima["gaus_a_2" ]);

  ParameterDict initialval;
  initialval["UnOscPdf_norm"]= 4000;
  initialval["d21" ]  = 7e-5;
  initialval["s12" ]  = 0.3;
  initialval["s13" ]  = 0.02;

  ParameterDict initialerr;
  initialerr["UnOscPdf_norm"]= 270;
  initialerr["d21" ] = 0.1*initialval["d21"];
  initialerr["s12" ] = 0.1*initialval["s12"];
  initialerr["s13" ] = 0.1*initialval["s13"];

  BinnedNLLH lhFunction;
  //BinnedOscNLLH lhFunction;
  lhFunction.SetDataDist(dataSetPdf); // initialise withe the data set
  //lhFunction.AddPdf(bgPdf);  //const std::string& name_,
  lhFunction.AddDist(mcPdfs.at(0));
  lhFunction.AddSystematic(osc_data);
  //lhFunction.AddPdf("osc", signalPdf);
  std::cout << "Built LH function " << std::endl;

  ////////////
  // 4. Fit //
  ////////////
  Minuit min;
  min.SetMethod("Migrad");
  min.SetMaxCalls(1000000);
  min.SetInitialValues(initialval);
  min.SetInitialErrors(initialerr);
  min.SetMinima(minima);
  min.SetMaxima(maxima);

  min.Optimise(&lhFunction);

  FitResult fitResult = min.GetFitResult();

  ParameterDict bestFit = fitResult.GetBestFit();
  fitResult.Print();

  /////////////////////////////////////////////

  BinnedED ResultHolder = mcPdfs.at(0);
  BinnedED Result;

  NuOsc OscResult("Oscillated");
  OscResult.SetFunction(new SurvProb(bestFit.at("d21"),bestFit.at("s12"),bestFit.at("s13"),238));
  OscResult.SetAxes(axes);
  OscResult.SetTransformationObs(dataRep);
  OscResult.SetDistributionObs(dataRep);
  OscResult.Construct();

  Result = OscResult(ResultHolder);
  Result.Scale(bestFit.at("UnOscPdf_norm"));

  TH1D DataHist;
  TH1D FitHist;

  DataHist = DistTools::ToTH1D(dataSetPdf);
  FitHist = DistTools::ToTH1D(Result);
  
  TH1D FullFit("FullFit","",FitHist.GetNbinsX(),FitHist.GetXaxis()->GetXmin(),FitHist.GetXaxis()->GetXmax());
  TH1D datahist("datahist","",FitHist.GetNbinsX(),FitHist.GetXaxis()->GetXmin(),FitHist.GetXaxis()->GetXmax());

  FullFit.Add(&FitHist);

  TCanvas* c1 = new TCanvas("c1");
  datahist.Draw("same");
  FullFit.Draw("same");
  
  TFile * fitout = new TFile("/home/blakei/FitOut.root","RECREATE");
  c1->Write();
  FullFit.Write();
  datahist.Write();
  
  fitout->Close();

  // save results
  //result.SaveAs("simpleFit_result_%d.txt");

  return 0;
}
