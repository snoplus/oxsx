// A very simple llh calculation for testing Barlow-Beeston
// Based on the example fit code
//
// Idea is you can hardcode bin contents of pdf and data
// and calculate LLHs with/without BB for comparison and
// validations.
// The BinnedNLLH will calculate the statistical uncertainty
// from the generated number of events for the PDF. Here we
// imagine we've simulated 1000 events, and scaled them down
// to produce the PDF
#include <BinnedED.h>
#include <ROOTNtuple.h>
#include <BinnedNLLH.h>
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
  axes.AddAxis(BinAxis("energy", 0, 4, 4, "Energy"));

  // Only interested in first bit of data ntuple
  std::vector<std::string> dataObs;
  dataObs.push_back("energy");

  // Set up pdf with these bins in this observable
  BinnedED bgPdf("bgPDF",axes);
  bgPdf.SetObservables(dataObs);

  //Hard code bin contents of PDF
  bgPdf.SetBinContent(0,20);
  bgPdf.SetBinContent(1,35);
  bgPdf.SetBinContent(2,30);
  bgPdf.SetBinContent(3,15);
  bgPdf.Normalise();

  //Hard code generated number of events
  const int gen_rate = 1000;

  std::cout << "Initialised PDF" << std::endl;

  /////////////////////////////////////
  // 2. Fill with data and normalise //
  /////////////////////////////////////

  BinnedED bgData("bgData",axes);
  bgData.SetObservables(dataObs);

  //Hard code bin contents of data
  bgData.SetBinContent(0,25);
  bgData.SetBinContent(1,30);
  bgData.SetBinContent(2,40);
  bgData.SetBinContent(3,10);

  std::cout << "Filled data " << std::endl;

  ////////////////////////////
  // 3. Set Up LH function  //
  ////////////////////////////
  BinnedNLLH lhFunction;
  lhFunction.SetDataDist(bgData); // initialise with the data set
  lhFunction.AddPdf(bgPdf, gen_rate);
  lhFunction.SetBarlowBeeston(true);

  lhFunction.RegisterFitComponents();

  std::cout << "Built LH function " << std::endl;

  ParameterDict parameterValues;

  //Set normalisation of PDF
  parameterValues["bgPDF"] = 100;

  //Calculate and print likelihood
  lhFunction.SetParameters(parameterValues);
  std::cout << "Total LogL = " << lhFunction.Evaluate() << std::endl;

  return 0;
}
