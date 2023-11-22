// A simple fit in energy for signal and a background
//
// See/run util/make_simple_data to create some simple fake data to run this example.
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
    axes.AddAxis(BinAxis("energy", 0, 10, 10, "Energy"));

    // Set up pdf with these bins in this observable
    BinnedED bgPdf("bgPDF",axes);
    bgPdf.SetObservables({"energy"});
    BinnedED  signalPdf("signalPDF",axes);
    signalPdf.SetObservables({"energy"});

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

    bgPdf.Normalise();
    signalPdf.Normalise();

    std::cout << "Filled pdfs " << std::endl;

    ////////////////////////////
    // 3. Set Up LH function  //
    ////////////////////////////
    ROOTNtuple dataNt(dataFile, dataTreeName);
    BinnedNLLH lhFunction;
    lhFunction.SetDataSet(&dataNt); // initialise withe the data set
    lhFunction.AddPdf(bgPdf);        
    lhFunction.AddPdf(signalPdf);        
  
    std::cout << "Built LH function " << std::endl;

    // Set up the optimisation
    GridSearch gSearch;

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

    ////////////
    // 4. Fit //
    ////////////
    FitResult result = gSearch.Optimise(&lhFunction);

    ParameterDict fit = result.GetBestFit();
    result.Print();
    result.SaveAs("simpleFit_result.txt");
    return 0;
}
