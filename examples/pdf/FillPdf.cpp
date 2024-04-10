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

    pdf.SetObservables(nt.GetObservableNames());

    // Now fill em up
    for(size_t i = 0; i < nt.GetNEntries(); i++)
        pdf.Fill(nt.GetEntry(i));
    
    // Done!

    // Lets do it again, but this time using the DistFiller
    DistFiller::FillDist(pdf, nt);

    return 0;
}
