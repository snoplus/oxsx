// A simple fit in energy for signal and a background
<<<<<<< HEAD
//
// See/run util/make_simple_data to create some simple fake data to run this example.
=======
>>>>>>> 56b3968c4b131a23445c6494e28bea263af24967
#include <BinnedED.h>
#include <ROOTNtuple.h>
#include <BinnedNLLH.h>
#include <GridSearch.h>
#include <ParameterDict.h>

const std::string bgMCfile    = "";
const std::string sigMCfile   = "";
const std::string bgTreeName  = "";
const std::string sigTreeName = "";

const std::string dataFile = "";
const std::string dataTreeName = "";

int main(){
    ////////////////////
    // 1. Set Up PDFs //
    ////////////////////

    // Set up binning
    AxisCollection axes;
<<<<<<< HEAD
    axes.AddAxis(BinAxis("energy", 0, 10, 10, "Energy"));
=======
    axes.AddAxis(BinAxis("energy", 2, 3, 10, "Energy"));
>>>>>>> 56b3968c4b131a23445c6494e28bea263af24967

    // Only interested in first bit of data ntuple
    ObsSet dataRep(0);

    // Set up pdf with these bins in this observable
    BinnedED bgPdf("bgPDF",axes);
    bgPdf.SetObservables(dataRep);
    BinnedED  signalPdf("signalPDF",axes);
    signalPdf.SetObservables(dataRep);

    std::cout << "Initialised Pdfs" << std::endl;

    /////////////////////////////////////
    // 2. Fill with data and normalise //
    /////////////////////////////////////

    ROOTNtuple bgMC(bgMCfile, bgTreeName);
    ROOTNtuple signalMC(sigMCfile, sigTreeName);

    for(size_t i = 0; i < bgMC.GetNEntries(); i++){
        bgPdf.Fill(bgMC.GetEntry(i));
    }

    for(size_t i = 0; i < signalMC.GetNEntries(); i++){
        signalPdf.Fill(signalMC.GetEntry(i));
    }

<<<<<<< HEAD
=======

>>>>>>> 56b3968c4b131a23445c6494e28bea263af24967
    bgPdf.Normalise();
    signalPdf.Normalise();

    std::cout << "Filled pdfs " << std::endl;

    ////////////////////////////
    // 3. Set Up LH function  //
    ////////////////////////////
    ROOTNtuple dataNt(dataFile, dataTreeName);
    BinnedNLLH lhFunction;
    lhFunction.SetDataSet(&dataNt); // initialise withe the data set
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
    lhFunction.AddDist(bgPdf);        
    lhFunction.AddDist(signalPdf);        
        
    std::cout << "Built LH function " << std::endl;        
         
    // Set up the optimisation        
    GridSearch gSearch;        
             
    std::vector<double> minima;
    minima.push_back(0);
    minima.push_back(0);
    std::vector<double> maxima;
    maxima.push_back(1000);
    maxima.push_back(1000);
    std::vector<double> stepsizes(2, 1);
         
    gSearch.SetMaxima(maxima);        
    gSearch.SetMinima(minima);        
    gSearch.SetStepSizes(stepsizes);        
             
    ////////////        
    // 4. Fit //        
    ////////////        
=======
=======
>>>>>>> 56b3968c4b131a23445c6494e28bea263af24967
=======
>>>>>>> 25c79512cf4d9978296c240c2cbcc4018b42a5dc
    lhFunction.AddPdf(bgPdf);
    lhFunction.AddPdf(signalPdf);

    std::cout << "Built LH function " << std::endl;

    // Set up the optimisation
    GridSearch gSearch;

<<<<<<< HEAD
<<<<<<< HEAD
    // Set up optimisation parameters.
    // Grid search needs a minimum, maximum and step size for each fit
    // parameter.
    // For BinnedEDs the normalisations are named by appending "_norm" to the
    // BinnedED name.
    ParameterDict minima;
    minima["bgPDF_norm"]= 800;
    minima["signalPDF_norm"]=800;

    ParameterDict maxima;
    maxima["bgPDF_norm"]= 2200;
    maxima["signalPDF_norm"]=2200;

    ParameterDict steps;
    steps["bgPDF_norm"]= 200;
    steps["signalPDF_norm"]=200;
    
    // Set optimisation parameters.
    gSearch.SetMinima(minima);
    gSearch.SetMaxima(maxima);
    gSearch.SetStepSizes(steps);
=======
=======
>>>>>>> 25c79512cf4d9978296c240c2cbcc4018b42a5dc
    ParameterDict values;
    values["minima"]= 0;
    values["maxima"]= 1000;
    values["steps"]= 2;
    
    gSearch.SetMaxima(values);
    gSearch.SetMinima(values);
    gSearch.SetStepSizes(values);
<<<<<<< HEAD
>>>>>>> 56b3968c4b131a23445c6494e28bea263af24967
=======
>>>>>>> 25c79512cf4d9978296c240c2cbcc4018b42a5dc

    ////////////
    // 4. Fit //
    ////////////
<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> 110ba72eb4c6cf0d2011a489391226981ddb4e7c
=======
>>>>>>> 56b3968c4b131a23445c6494e28bea263af24967
=======
>>>>>>> 25c79512cf4d9978296c240c2cbcc4018b42a5dc
    FitResult result = gSearch.Optimise(&lhFunction);

    ParameterDict fit = result.GetBestFit();
    result.Print();
    result.SaveAs("simpleFit_result.txt");
    return 0;
}
