/*
  The purpose of a data representation is to hand shake objects like
  pdfs, systemtics etc. that are only typically interested in a subset of 
  the observables in a dataset.

  For example, we might produce an ntuple with energy, radius, time
  and a flag to indicate fit validity. We typically will cut on the 
  fit validity and perform a fit on the other three observables.
  So, we equip each of our PDFs with a 3D data representation, that tells
  them where to get the data they need from any given event.

  Observables are referred to by their sequential position in the dataset,
  but if you have the dataset you are reading from it can produce a 
  representation from a vector of names, both are shown below

 */

#include <DataRepresentation.h>
#include <ROOTNtuple.h>
#include <BinnedPdf.h>
#include <iostream>

int main(){
    DataRepresentation dataRep1(0); 
    // if only concerned in index 0, if you were reading from a root ntuple
    // objects with this dataRep will automatically just look at the first
    // branch of an event

    // initialisation from a vector of indices
    std::vector<size_t> indices;
    indices.push_back(1);
    indices.push_back(2);
    DataRepresentation dataRep2(indices);
    // will look at observables at indices 1 and 2

    // ROOTNtuple nt("myfilename.root", "mytreename");
    std::vector<std::string> interestingObservablesForMyObject;
    interestingObservablesForMyObject.push_back("ob1");
    interestingObservablesForMyObject.push_back("ob2");

    // DataRepresentation dataRepForMyObject = nt.MakeDataRep(interestingObservablesForMyObject);
   
    // Now test it out on some other fake objects
    AxisCollection axes;
    axes.AddAxis(PdfAxis("", 0, 10, 100));
    BinnedPdf pdf(axes);
    pdf.SetDataRep(dataRep1);

    Event event(std::vector<double> (20, 1));
    std::cout << "A " << dataRep1.GetNObservables() << "D representation\n"
              << "allows a " << pdf.GetNDims() << "D pdf\n"
              << "to act on an event with " 
              << event.GetNObservables()
              << " observables\n"
              << "To give probability " << pdf.Probability(event)
              << std::endl;

    
    return 0;
}
