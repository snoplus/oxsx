/*
Fill a pdf with some data from a root ntuple. We do this two ways:
First the long way around to see what's happening, then in one line with 
DistFiller
*/

#include <BinnedED.h>
#include <ROOTNtuple.h>
#include <DistFiller.h>
#include <vector>
#include <string>

const std::string filename = "";
const std::string treename = "";
const std::string obs1 = "obs1name";
const std::string obs2 = "obs2name"; 
// this is the name inside the dataset, e.g. for a ROOTNtuple its the branch 
// name.

int main(){
    // Open up the data file
    ROOTNtuple nt(filename, treename);

    // Create a 2D binned pdf, axes named obs1, obs2
    AxisCollection axes;
    axes.AddAxis(BinAxis("obs1", 0, 10, 10));
    axes.AddAxis(BinAxis("obs2", 0, 10, 10));

    BinnedED pdf("pdf", axes);
    
    // now link up the pdf to the two observables we want, by name
    std::vector<std::string> relevantObs;
    relevantObs.push_back(obs1);
    relevantObs.push_back(obs2);

    pdf.SetObservables(nt.MakeDataRep(relevantObs));

    // Now fill em up
    for(size_t i = 0; i < nt.GetNEntries(); i++)
        pdf.Fill(nt.GetEntry(i));
    
    // Done!

    // Lets do it again, but this time using the DistFiller
    DistFiller::FillDist(pdf, nt);

    return 0;
}
