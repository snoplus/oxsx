// A simple fit in energy for signal and a background
#include <string>
#include <vector>
#include <math.h>
#include <Rand.h>
#include <fstream>
#include <iostream>

#include <TH1D.h>
#include <TH2D.h>
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


void readInfoFile(const std::string &runInfoFileName, std::vector<std::string> &reactorNames, std::vector<double> &distances, std::vector<std::string> &reactorTypes, std::vector<int> &nCores, std::vector<double> &powers ) {
    // Read couchDB run-info text file
    std::ifstream in;
    in.open(runInfoFileName.c_str());
    std::cout << "opening file: " << runInfoFileName.c_str() << std::endl;

    std::fill(reactorNames.begin(), reactorNames.end(), "");
    std::fill(distances.begin(), distances.end(), 0.);
    std::fill(reactorTypes.begin(), reactorTypes.end(), "");
    std::fill(nCores.begin(), nCores.end(), 0);
    std::fill(powers.begin(), powers.end(), 0.);

    std::string reactorName,distance,reactorType,nCore,power;
    int lineNo = 0;

    // read until end of file.
    while(in.good()){
        std::getline(in,reactorName,',');
        std::getline(in,distance,',');
        std::getline(in,reactorType,',');
        std::getline(in,nCore,',');
        std::getline(in,power,'\n');

        if (lineNo>0){ //skip csv header
            if (strcmp(reactorName.c_str(),"")!=0) {
                reactorNames.push_back(reactorName);
                distances.push_back(atof(distance.c_str()));
                reactorTypes.push_back(reactorType.c_str());
                nCores.push_back(atoi(nCore.c_str()));
                powers.push_back(atof(power.c_str()));

                //std::cout << "v: reactorName: " << reactorNames[lineNo-1] << ", distance: " << distances[lineNo-1] << ", reactorType: " << reactorTypes[lineNo-1] << ", nCore: " << nCores[lineNo-1] << ", power: " << powers[lineNo-1] << std::endl; //debug check ('-1' for header)
            }
        }
        lineNo++;
    }
    in.close();
}

double LHFit(TH2D *h2, const std::string UnOscfile, const std::string dataFile, int numPdfs, std::vector<double> &reactorDistances, double param_d21, double param_s12, double param_s13){

  char name[100];

  TRandom3 *r1 = new TRandom3();
  r1->SetSeed(0);

  ////////////////////
  // 1. Set Up PDFs //
  ////////////////////

  // Only interested in first bit of data ntuple
  ObsSet dataRep(0);

  // Set up binning
  AxisCollection axes;
  double Emin = 2;
  double Emax = 8;
  int numbins = 60;
  axes.AddAxis(BinAxis("ParKE", Emin, Emax, numbins));

  BinnedED dataSetPdf("dataSetPdf",axes);
  dataSetPdf.SetObservables(dataRep);
  ROOTNtuple dataNtp(dataFile, "nt");

  for(size_t i = 0; i < dataNtp.GetNEntries(); i++)
    dataSetPdf.Fill(dataNtp.GetEntry(i));

  //poisson stat fluctutate input data
  for(int i = 0; i < dataSetPdf.GetNBins(); i++)
      dataSetPdf.SetBinContent(i, r1->Poisson(dataSetPdf.GetBinContent(i)));

  std::vector<BinnedED*> reactorPdf(numPdfs);
  ROOTNtuple reactorNtp(UnOscfile, "nt");
  NuOsc *reactorSystematic[numPdfs];
  ParameterDict minima;
  ParameterDict maxima;
  ParameterDict initialval;
  ParameterDict initialerr;

  BinnedNLLH lhFunction;
  lhFunction.SetBufferAsOverflow(true);
  int Buff = 5;
  lhFunction.SetBuffer(0,Buff,Buff);
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

    for(size_t j = 0; j < reactorNtp.GetNEntries(); j++)
      reactorPdf[i]->Fill(reactorNtp.GetEntry(j));

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

  //Fit
  Minuit min;
  min.SetMethod("Migrad");
  min.SetMaxCalls(10000000);
  min.SetMinima(minima);
  min.SetMaxima(maxima);
  min.SetInitialValues(initialval);
  min.SetInitialErrors(initialerr);

  //Fit Result
  FitResult fitResult = min.Optimise(&lhFunction);
  ParameterDict bestFit = fitResult.GetBestFit();
  fitResult.Print();

  double lhmax = -1000000000;
  double lhmin = 1000000000;
  lhFunction.SetParameters(bestFit);
  double lhval =(-1)*lhFunction.Evaluate();
  if (lhval > lhmax)
    lhmax = lhval;
  if (lhval < lhmin)
    lhmin = lhval;
  std::cout<<"LH val at best fit: "<<lhval<<std::endl;
  return lhval;
}

int main(int argc, char *argv[]) {

    if (argc != 5) {
        std::cout<<"4 arguments expected."<<std::endl;
    }
    else{
        std::stringstream argParser;
        const std::string &UnOscfile = argv[1];
        const std::string &dataFile = argv[2];
        const std::string &infoFile = argv[3];
        const std::string &outFile = argv[4];

        std::vector<std::string> reactorNames;
        std::vector<double> distances;
        std::vector<std::string> reactorTypes;
        std::vector<int> nCores;
        std::vector<double> powers;
        readInfoFile(infoFile, reactorNames, distances, reactorTypes, nCores, powers); // get reactor information
        int numPdfs = reactorNames.size();

        printf("\n-------------------------------------------------------------------------------------\n");
        //double lhValue = LHFit(UnOscfile,dataFile,numPdfs,distances,d21,s12,0.0215)
        printf("-------------------------------------------------------------------------------------\n");
    }
    return 1.;//lhValue;
}
