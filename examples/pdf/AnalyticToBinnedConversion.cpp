/*
  Demonstrates how to produce a binned pdf from a function, and a defined
  binning The pdf converter  stored in DistTools simply integrates the function over the
  bin boundaries, so the function passed in has to implement
  the IntegrableFunction interface.
*/
#include <Gaussian.h>
#include <DistTools.h>
#include <BinnedED.h>

int main()
{
  // Here the example is in 1D, for ND just use an ND function and
  // binning defined in ND

  // 1. Get function
  Gaussian gaus(0, 1); // mean = 0, stdDev = 1

  // 2. Define the binning
  AxisCollection axes;
  axes.AddAxis(BinAxis("fake", 0, 10, 100, "")); // min, max, nbins
  BinnedED pdf("binnedPDF", DistTools::ToHist(gaus, axes));

  // 3. Dont forget to set the data representation if you want to use it
  // to calculate event probabilities
  const std::vector<std::string> observables = {"energy"};
  pdf.SetObservables(observables);

  Event event(Event(std::vector<double>(1, 2)));
  event.SetObservableNames(&observables);

  std::cout << pdf.Probability(event) << std::endl;
}
