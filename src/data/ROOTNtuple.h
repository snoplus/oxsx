/****************************************************************************************/
/* Hands out data from a simple root ntuple into an Event object for each Event */
/* reads directly to save the copy (though this happens only once anyway)               */
/****************************************************************************************/

#ifndef __OXSX_ROOT_NTUPLE__
#define __OXSX_ROOT_NTUPLE__
#include <DataSet.h>
#include <string>
#include <TFile.h>
#include <vector>

class Event;
class TNtuple;
class ROOTNtuple : public DataSet{
 public:
    ROOTNtuple(const std::string& fileName_, const std::string& treeName_);
    ~ROOTNtuple();

    Event GetEntry(size_t iEvent_) const;
    unsigned  GetNEntries() const;
    unsigned  GetNObservables() const;

    void LoadBaskets();
    void DropBaskets();

    std::vector<std::string> GetObservableNames() const; // just returns the cache

 private:
    void GatherObservableNames(); // actually works the names out
    std::vector<std::string> fObsNames;

    // no copying allowed
    ROOTNtuple(const ROOTNtuple&);
    ROOTNtuple operator=(const ROOTNtuple&);

    TFile*   fROOTFile;
    TNtuple* fNtuple;
    
    Event Assemble(size_t iEvent_) const;
};
#endif
