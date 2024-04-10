/*
  Demonstrates how to read data from a ROOTNtuple
 */
#include <ROOTNtuple.h>
#include <ObsSet.h>
#include <string>
#include <iostream>

const std::string fileName = ""; 
const std::string treeName = ""; 

// used for creating data representations, both should match an observable
// name from the root ntuple (i.e. a branch name)

int main(){
    ROOTNtuple nt(fileName, treeName);
    
    // The ntuple knows how many entries it has, and what the observables are
    std::cout << "The ntuple contains " << nt.GetNEntries() << "entries "
              << "\n each with " << nt.GetNObservables() 
              << "called:" << std::endl;

    for(size_t i = 0; i < nt.GetNObservables(); i++)
        std::cout << nt.GetObservableNames().at(i) << "\t";
    std::cout << std::endl;

    ObsSet dataRepForMyPDF = nt.GetObservableNames();
    // Call pdf.SetDataRep(dataRepForMyPDF) so that your pdf automatically
    // picks out <observableName1> and <observableName2> from events

    return 0;
}
