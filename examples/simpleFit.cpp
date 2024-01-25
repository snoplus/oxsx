// A simple fit in energy for signal and a background        
#include <BinnedPdf.h>        
#include <ROOTNtuple.h>        
#include <BinnedNLLH.h>        
#include <GridSearch.h>        

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
    axes.AddAxis(PdfAxis("energy", 2, 3, 10, "Energy"));
         
    // Only interested in first bit of data ntuple        
    DataRepresentation dataRep(0);        
         
    // Set up pdf with these bins in this observable        
    BinnedPdf bgPdf(axes);      bgPdf.SetDataRep(dataRep);        
    BinnedPdf signalPdf(axes);  signalPdf.SetDataRep(dataRep);        
         
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
    FitResult result = gSearch.Optimise(&lhFunction);
         
    std::vector<double> fit = result.GetBestFit();        
    result.Print();
    result.SaveAs("simpleFit_result.txt");
    return 0;        
}
