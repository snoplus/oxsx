/*
  Demonstrates the use of a BinnedPdf. A BinnedPdf is just a histogram, 
  plus a data representation that tells it where to look in an event.
  
  We first define the binning, then use it to create a binnedPdf. Set the 
  data representation, and then call the pdf on a plain vector input, and 
  to evaluate the probability of an event
*/
#include <BinnedPdf.h>
#include <iostream>

int main(){
    // Define the binning
    AxisCollection axes; // defines the bin boundaries
    axes.AddAxis(PdfAxis("axis name", 0, 10, 100,         // min, max, nbins
                         "optional nice name for latex"));

    axes.AddAxis(PdfAxis("axis name2", 0, 10, 100, 
                         "optional nice name for latex2"));
    
    // Create the PDF
    BinnedPdf binnedPdf(axes);

    // Dimensionality matches the binning
    std::cout << "Pdf is "<< binnedPdf.GetNDims() << " dimensions" 
              << std::endl;

    // Total bin number is calculated automatically
    std::cout << "with " << binnedPdf.GetNBins() << " bins"
              << std::endl;

    //  every bin is assigned a unique ID
    std::cout << "bin (5, 5) is equivilent to bin #" 
              << binnedPdf.FlattenIndices(std::vector<size_t>(2, 5))
              << std::endl;

    std::vector<size_t> binIndices = binnedPdf.UnpackIndices(505);
    std::cout << "Global bin #55 is equivilent to (" 
              << binIndices.at(0) << "," << binIndices.at(1) << ")"
              << std::endl;

    // and has an associated bin content, stored in a vector
    std::cout << "there are " << binnedPdf.GetBinContents().size()
              << " bin content values"
              << std::endl;

    // its now ready to call on plain old data for testing
    binnedPdf.Fill(std::vector<double>(2, 1));
    binnedPdf.Normalise();
    std::cout << binnedPdf(std::vector<double>(2,1)) << std::endl;

    // If we also want to use it on events, we set the data representation
    // our pdf should extract the 0th and 2nd indices and ignore the rest
    // obviously there should be the same number of indices in the data rep
    // as there are axes in the pdf
    std::vector<size_t> relevantIndices;
    relevantIndices.push_back(0);
    relevantIndices.push_back(2);
    binnedPdf.SetDataRep(DataRepresentation(relevantIndices));

    // Now its ready for use on events, make a fake one here with 10 obs
    // the call to Probability automatically selects the right indices,
    // this allows different pdfs to operate on different observables
    Event fakeEvent(std::vector<double>(10,1));
    binnedPdf.Fill(fakeEvent);
    std::cout << binnedPdf.Probability(fakeEvent) << std::endl;
    

    // Marginalisation is possible, just pass the indices you would like to
    // _keep_. Here the options are 0 or 2 (think data representation)
    BinnedPdf projection = binnedPdf.Marginalise(std::vector<size_t>(1,2));
    std::cout << "the projection along index 2 is " 
              << projection.GetNDims() << " dimensional"
              << std::endl;


    return 0;
}
