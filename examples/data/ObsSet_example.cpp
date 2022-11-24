/*
  The purpose of a data representation is to hand shake objects like
  pdfs, systemtics etc. that are only typically interested in a subset of 
  the observables in a dataset.

  For example, we might produce an ntuple with energy, radius, time
  and a flag to indicate fit validity. We typically will cut on the 
  fit validity and perform a fit on the other three observables.
  So, we equip each of our PDFs with a 3D data representation, that tells
  them where to get the data they need from any given event.

  The type we use here to store the data representation is the ObsSet class.

  Observables are referred to by their name.

 */

#include <ObsSet.h>
#include <ROOTNtuple.h>
#include <BinnedED.h>
#include <iostream>

int main(){
    // initialisation from a vector of strings
    std::vector<std::string> observables = {"energy", "radius"};
    ObsSet dataRep2(observables);
   
    // Now test it out on some other fake objects
    AxisCollection axes;
    axes.AddAxis(BinAxis("", 0, 10, 100));
    BinnedED pdf("", axes);
    pdf.SetObservables(observables);

    Event event(std::vector<double> (20, 1));
    std::cout << "A " << dataRep2.GetNObservables() << "D representation\n"
              << "allows a " << pdf.GetNDims() << "D pdf\n"
              << "to act on an event with " 
              << event.GetNObservables()
              << " observables\n"
              << "To give probability " << pdf.Probability(event)
              << std::endl;

    
    return 0;
}
