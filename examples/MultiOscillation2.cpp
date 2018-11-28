// A fit in energy for signal and a background
#include <string>
#include <vector>
#include <math.h>
#include <Rand.h>
#include <fstream>
#include <iostream>

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

double LHFit(const std::string UnOscfile, const std::string dataFile, int numPdfs, std::vector<double> &reactorDistances, double param_d21, double param_s12, double param_s13){

  char name[100];

  // Only interested in first bit of data ntuple
  ObsSet dataRep(0);

  // Set up binning
  AxisCollection axes;
  double Emin = 2;
  double Emax = 8;
  int numbins = 60;
  axes.AddAxis(BinAxis("ParKE", Emin, Emax, numbins));

  // create and fill data ntp and pdf
  BinnedED dataSetPdf("dataSetPdf", axes);
  dataSetPdf.SetObservables(dataRep);
  ROOTNtuple dataNtp(dataFile, "nt");
  for(size_t i = 0; i < dataNtp.GetNEntries(); i++)
    dataSetPdf.Fill(dataNtp.GetEntry(i));
  dataSetPdf.Normalise();

  // create and fill simulated ntp and pdf
  BinnedED *reactorPdf0 = new BinnedED(name, axes);
  reactorPdf0->SetObservables(0);
  ROOTNtuple reactorNtp(UnOscfile, "nt");
  for(size_t i = 0; i < reactorNtp.GetNEntries(); i++)
    reactorPdf0->Fill(reactorNtp.GetEntry(i));
  reactorPdf0->Normalise();

  ParameterDict minima;
  ParameterDict maxima;
  ParameterDict initialval;
  ParameterDict initialerr;

  BinnedNLLH lhFunction;
  lhFunction.SetBufferAsOverflow(true);
  int Buff = 5;
  lhFunction.SetBuffer(0, Buff, Buff);
  lhFunction.SetDataDist(dataSetPdf); // initialise with the data set

  // loop over all reactor pdfs
  for (int i = 0; i < numPdfs; i++){
    sprintf(name, "ReactorPdf%d", i);
    BinnedED *reactorPdf = new BinnedED(name, axes);
    reactorPdf->SetObservables(0);

    NuOsc reactorSystematic("reactorSystematic");
    reactorSystematic.SetFunction(new SurvProb(param_d21, param_s12, param_s13, reactorDistances[i]));
    reactorSystematic.SetAxes(axes);
    reactorSystematic.SetTransformationObs(dataRep);
    reactorSystematic.SetDistributionObs(dataRep);
    reactorSystematic.Construct();

    reactorPdf->Add(reactorSystematic(*reactorPdf0), 1);
    reactorPdf->Normalise();

    // Setting optimisation limits
    sprintf(name,"ReactorPdf%d_norm", i);
    minima[name] = 0;
    maxima[name] = 100000;
    initialval[name] = 50000;
    initialerr[name] = 0.1*initialval[name];

    lhFunction.AddDist(*reactorPdf);
  }

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
  lhFunction.SetParameters(bestFit);
  double lhval =(-1)*lhFunction.Evaluate();
  return lhval;
}

int main(int argc, char *argv[]) {

  if (argc != 8){
      std::cout<<"Error: 7 arguments expected."<<std::endl;
      return 1; // failed
  }
  else{

    //std::stringstream argParser;
    const std::string &UnOscfile = argv[1];
    const std::string &dataFile = argv[2];
    const std::string &infoFile = argv[3];
    const std::string &outFile = argv[4];
    double d_21 = atof(argv[5]);
    double s_12 = atof(argv[6]);
    double s_13 = atof(argv[7]);

    printf("LHFitting:: del_21:%.7f, sin2_12:%.7f, sin2_13:%.7f\n",d_21,s_12,s_13);

    std::vector<std::string> reactorNames;
    std::vector<double> distances;
    std::vector<std::string> reactorTypes;
    std::vector<int> nCores;
    std::vector<double> powers;
    readInfoFile(infoFile, reactorNames, distances, reactorTypes, nCores, powers); // get reactor information
    int numPdfs = reactorNames.size();

    double lhValue = LHFit(UnOscfile,dataFile,numPdfs,distances,d_21,s_12,s_13);
    //TRandom3 *myRand = new TRandom3() ;
    //myRand->SetSeed(0);
    //double lhValue = myRand->Gaus(50,5);
    int fitValidity = 1; //make this do something useful
    printf("LHValue:%.7f\n",lhValue);

    //Write fit coefficients to txt file
    printf("writing to: %s\n",outFile.c_str());
    FILE *fOut = fopen(outFile.c_str(),"w");

    printf("fit valid: %d\n",fitValidity);
    fprintf(fOut,"fit valid: %d\n", fitValidity);

    printf("d21,s12,s13,lhValue\n");
    fprintf(fOut,"d21,s12,s13,lhValue\n", d_21, s_12, s_13, lhValue);
    printf("%.7f,%.7f,%.7f,%.7f\n", d_21, s_12, s_13, lhValue);
    fprintf(fOut,"%.7f,%.7f,%.7f,%.7f\n", d_21, s_12, s_13, lhValue);

    fclose(fOut);

    return 0; // completed successfully
  }
}
