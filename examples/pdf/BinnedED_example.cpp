/*
  Demonstrates the use of a BinnedED. A BinnedED is just a histogram, 
  plus a data representation that tells it where to look in an event.
  
  We first define the binning, then use it to create a BinnedED. Set the 
  data representation, and then call the pdf on a plain vector input, and 
  to evaluate the probability of an event.
*/
#include <BinnedED.h>
#include <iostream>

int main(){
    // Define the binning
    AxisCollection axes; // defines the bin boundaries
    axes.AddAxis(BinAxis("axis name", 0, 10, 100,         // min, max, nbins
                         "optional nice name for latex"));

    axes.AddAxis(BinAxis("axis name2", 0, 10, 100, 
                         "optional nice name for latex2"));
    
    // Create the PDF
    BinnedED binnedED("", axes);

    // Dimensionality matches the binning
    std::cout << "Pdf is "<< binnedED.GetNDims() << " dimensions" 
              << std::endl;

    // Total bin number is calculated automatically
    std::cout << "with " << binnedED.GetNBins() << " bins"
              << std::endl;

    //  every bin is assigned a unique ID
    std::cout << "bin (5, 5) is equivalent to bin #" 
              << binnedED.FlattenIndices(std::vector<size_t>{5, 5})
              << std::endl;

    std::vector<size_t> binIndices = binnedED.UnpackIndices(505);
    std::cout << "Global bin #505 is equivilent to (" 
              << binIndices.at(0) << "," << binIndices.at(1) << ")"
              << std::endl;

    // and has an associated bin content, stored in a vector
    std::cout << "there are " << binnedED.GetBinContents().size()
              << " bin content values"
              << std::endl;

    // its now ready to call on plain old data for testing
    binnedED.Fill(std::vector<double>(2, 1));
    binnedED.Normalise();
    std::cout << binnedED.Probability(std::vector<double>(2,1)) << std::endl;

    // If we also want to use it on events, we set the data representation
    // our pdf should extract the "obs0" and "obs1" observables and
    // ignore the rest.
    // obviously there should be the same number of indices in the data rep
    // as there are axes in the pdf
    std::vector<std::string> observables = {"obs0", "obs1"};
    binnedED.SetObservables(observables);

    // Now its ready for use on events, make a fake one here with 10 obs
    // the call to Probability automatically selects the right indices,
    // this allows different pdfs to operate on different observables
    Event fakeEvent(std::vector<double>(10,1));
    fakeEvent.SetObservableNames(&observables);
    binnedED.Fill(fakeEvent);
    std::cout << binnedED.Probability(fakeEvent) << std::endl;
    

    // Marginalisation is possible, just pass the indices you would like to
    // _keep_. Here the options are 0 or 2 (think data representation)
    BinnedED projection = binnedED.Marginalise(std::vector<std::string>{"axis name2"});
    std::cout << "the projection along observable 0 is " 
              << projection.GetNDims() << " dimensional"
              << std::endl;


    return 0;
}
