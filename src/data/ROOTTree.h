/****************************************************************************************/
/* Hands out data from a ROOT TTree into an Event object for each Event.                */
/* Reads directly to save the copy (though this happens only once anyway)               */
/* Requires all branches in TTree to be castable to a double.                            */
/****************************************************************************************/

#ifndef __OXSX_ROOT_TREE__
#define __OXSX_ROOT_TREE__
#include <DataSet.h>
#include <string>
#include <TFile.h>
#include <vector>

class Event;
class TTree;
class TLeaf;
class ROOTTree : public DataSet{
 public:
    ROOTTree(const std::string& fileName_, const std::string& treeName_);
    ~ROOTTree();

    Event GetEntry(size_t iEvent_) const;
    unsigned  GetNEntries() const;
    unsigned  GetNObservables() const;

    void LoadBaskets();
    void DropBaskets();

    std::vector<std::string> GetObservableNames() const; // just returns the cache

 private:
    void GatherObservableNames(); // actually works the names out
    std::vector<std::string> fObsNames;
    std::vector<TLeaf*> fLeaves;

    // no copying allowed
    ROOTTree(const ROOTTree&);
    ROOTTree operator=(const ROOTTree&);

    TFile*   fROOTFile;
    TTree* fTree;
    
    Event Assemble(size_t iEvent_) const;
};
#endif
