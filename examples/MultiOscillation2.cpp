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
#include <TRandom3.h>
#include <TH1D.h>

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

void readInfoFile(const std::string &runInfoFileName, std::vector<std::string> &reactorNames, std::vector<Double_t> &distances, std::vector<std::string> &reactorTypes, std::vector<ULong64_t> &nCores, std::vector<Double_t> &powers ) {
    // Read couchDB run-info text file
    std::ifstream in;
    in.open(runInfoFileName.c_str());
    //std::cout << "opening file: " << runInfoFileName.c_str() << std::endl;

    std::fill(reactorNames.begin(), reactorNames.end(), "");
    std::fill(distances.begin(), distances.end(), 0.);
    std::fill(reactorTypes.begin(), reactorTypes.end(), "");
    std::fill(nCores.begin(), nCores.end(), 0);
    std::fill(powers.begin(), powers.end(), 0.);

    std::string reactorName,distance,reactorType,nCore,power;
    ULong64_t lineNo = 0;

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

Double_t LHFit(const std::string &inPath, std::vector<std::string> &reactorNames, std::vector<Double_t> &distances, std::vector<std::string> &reactorTypes, std::vector<ULong64_t> &nCores, std::vector<Double_t> &powers, Double_t param_d21, Double_t param_s12, Double_t param_s13){

    char name[100];

    printf("Reading data from: %s\n", inPath.c_str());

    ULong64_t flux = 100;
    sprintf(name, "/data/snoplusmc/lidgard/OXSX_kamland/flux%llu/all_flux%llu_day360_cleanround1_oxsx.root", flux, flux);
    const std::string dataFile = std::string(name);
    printf("Loading data pdf #%llu: %s\n",0,dataFile.c_str());

    // Only interested in first bit of data ntuple
    ObsSet dataRep(0);

    // Set up binning
    AxisCollection axes;
    Double_t Emin = 2;
    Double_t Emax = 8;
    Int_t numbins = 60;
    axes.AddAxis(BinAxis("mc_neutrino_energy", Emin, Emax, numbins));

    // // create and fill data ntp and pdf
    BinnedED dataSetPdf("dataSetPdf", axes);
    dataSetPdf.SetObservables(dataRep);
    ROOTNtuple dataNtp(dataFile, "nt");
    for(ULong64_t j = 0; j < dataNtp.GetNEntries(); j++)
        dataSetPdf.Fill(dataNtp.GetEntry(j));
    dataSetPdf.Normalise();

    NuOsc *reactorSystematic;
    ParameterDict minima;
    ParameterDict maxima;
    ParameterDict initialval;
    ParameterDict initialerr;

    BinnedNLLH lhFunction;
    lhFunction.SetBufferAsOverflow(true);
    int Buff = 5;
    lhFunction.SetBuffer(0, Buff, Buff);
    lhFunction.SetDataDist(dataSetPdf); // initialise with the data set
    
    const ULong64_t n_pdf = reactorNames.size();
    BinnedED *reactorPdf[n_pdf];

    //loop over all reactor pdfs
    for (ULong64_t i = 0; i < n_pdf; i++){
        
        // create and fill simulated ntp and pdf
        flux = 100;
        sprintf(name, "flux%llu/%s_flux%llu_day360_cleanround1_oxsx.root", flux, reactorNames[i].c_str(), flux);
        std::string UnOscfile = inPath + std::string(name);
        printf("Loading unoscillated pdf #%llu: %s\n", i, UnOscfile.c_str());
        
        sprintf(name, "ReactorPdf%llu", i);
        printf("ReactorPdf%llu: %s\n", i, reactorNames[i].c_str());
        reactorPdf[i] = new BinnedED(name, axes);
        reactorPdf[i]->SetObservables(0);
        ROOTNtuple reactorNtp(UnOscfile, "nt");
        for(ULong64_t j = 0; j < reactorNtp.GetNEntries(); j++)
            reactorPdf[i]->Fill(reactorNtp.GetEntry(j));
        reactorPdf[i]->Normalise();

        NuOsc reactorSystematic("reactorSystematic");
        SurvProb *surv_prob = new SurvProb(param_d21, param_s12, distances[i]);
        surv_prob->Setsinsqrtheta13s(param_s13); // manual, fixed setting of theta13
        reactorSystematic.SetFunction(surv_prob);
        reactorSystematic.SetAxes(axes);
        reactorSystematic.SetTransformationObs(dataRep);
        reactorSystematic.SetDistributionObs(dataRep);
        reactorSystematic.Construct();

        reactorPdf[i]->Add(reactorSystematic(*reactorPdf[i]), 1);

        // Setting optimisation limits
        sprintf(name,"ReactorPdf%llu_norm", i);
        minima[name] = 0;
        maxima[name] = 1000;
        initialval[name] = 10;
        initialerr[name] = 0.1*initialval[name];

        lhFunction.AddDist(*reactorPdf[i]);
    }

    //Fit
    Minuit min;
    min.SetMethod("Migrad");
    min.SetMaxCalls(10000000);
    min.SetMinima(minima);
    min.SetMaxima(maxima);
    min.SetInitialValues(initialval);
    min.SetInitialErrors(initialerr);

    // //Fit Result
    FitResult fitResult = min.Optimise(&lhFunction);
    ParameterDict bestFit = fitResult.GetBestFit();
    fitResult.Print();
    lhFunction.SetParameters(bestFit);
    Double_t lhval =(-1)*lhFunction.Evaluate();

    //Double_t lhval = 0;
    return lhval;
}

int main(int argc, char *argv[]) {

    if (argc != 7){
        std::cout<<"Error: 6 arguments expected."<<std::endl;
        return 1; // return>0 indicates error code
    }
    else{
        const std::string &inPath = argv[1];
        const std::string &infoFile = argv[2];
        const std::string &outFile = argv[3];
        Double_t d_21 = atof(argv[4]);
        Double_t s_12 = atof(argv[5]);
        Double_t s_13 = atof(argv[6]);

        printf("LHFitting:: del_21:%.7f, sin2_12:%.7f, sin2_13:%.7f\n",d_21,s_12,s_13);

        std::vector<std::string> reactorNames;
        std::vector<Double_t> distances;
        std::vector<std::string> reactorTypes;
        std::vector<ULong64_t> nCores;
        std::vector<Double_t> powers;
        readInfoFile(infoFile, reactorNames, distances, reactorTypes, nCores, powers); // get reactor information
        //printf("numPdfs:%llu\n",reactorNames.size());

        // // print out read info
        // for (size_t i=0; i<(size_t)reactorNames.size(); i++){
            // printf("i:%llu,reactorNames[i]:%s, distance: %.5f, type: %s, nCores: %llu, power: %.5f \n",i,reactorNames[i].c_str(),distances[i],reactorTypes[i].c_str(),nCores[i],powers[i]);
        // }

        Double_t lhValue = LHFit(inPath, reactorNames, distances, reactorTypes, nCores, powers, d_21, s_12, s_13);
        //input_path_unosc_i = input_path_100+"all_flux100_day360_cleanround1_oxsx.root"
        //input_path_data_i = input_path_100+osc_label+"/all_flux100_day360_cleanround1_"+osc_label+"_oxsx.root"

        //TRandom3 *myRand = new TRandom3() ;
        //myRand->SetSeed(0);
        //Double_t lhValue = myRand->Gaus(50,5);

        ULong64_t fitValidity = 1; //make this do something useful
        printf("LHValue:%.5f\n",lhValue);

        //Write fit coefficients to txt file
        //printf("writing to: %s\n",outFile.c_str());
        FILE *fOut = fopen(outFile.c_str(),"w");

        //printf("fit valid: %llu\n",fitValidity);
        fprintf(fOut,"fit valid: %llu\n", fitValidity);

        //printf("d21,s12,s13,lhValue\n");
        fprintf(fOut,"d21,s12,s13,lhValue\n", d_21, s_12, s_13, lhValue);
        //printf("%.9f,%.7f,%.7f,%.5f\n", d_21, s_12, s_13, lhValue);
        fprintf(fOut,"%.9f,%.7f,%.7f,%.5f\n", d_21, s_12, s_13, lhValue);
        fclose(fOut);
        return 0; // completed successfully
    }
}
