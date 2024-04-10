#include <ROOTTree.h>
#include <Event.h>
#include <iostream>
#include <Exceptions.h>
#include <TTree.h>
#include <TLeaf.h>
#include <Formatter.hpp>

ROOTTree::ROOTTree(const std::string& fileName_, const std::string& treeName_){
    fROOTFile = new TFile(fileName_.c_str());

    if (fROOTFile->IsZombie()){
        delete fROOTFile;
        throw IOError("ROOTTree::File Does not Exist! or is Zombie " + fileName_);
    }

    fTree = dynamic_cast<TTree*>(fROOTFile -> Get(treeName_.c_str()));

    if(!fTree){
        delete fROOTFile;
        throw IOError(Formatter()<<"ROOTTree::Tree does not exist, or isn't a TTree! tree : " << treeName_ << ", filename: "<<fileName_);
    }        
    GatherObservableNames();
}

ROOTTree::~ROOTTree(){
    if (fROOTFile)
        fROOTFile -> Close();
    delete fROOTFile;
}


std::vector<std::string>
ROOTTree::GetObservableNames() const{
    return fObsNames;
}

Event 
ROOTTree::Assemble(size_t iEvent_) const{
    // First, check we're within the range of possible events
    if (iEvent_ >= GetNEntries())
        throw NotFoundError(Formatter() << "Exceeded end of ROOT TTree"
                            << " \n\t(requested " << iEvent_ 
                            << " but only have " << GetNEntries() << ")");
    fROOTFile->cd();
    // Fill the TLeaves with the info from the ith event in the tree
    fTree -> GetEntry(iEvent_);
    // Now, get that info from those leaves!
    const size_t nObs = GetNObservables();
    std::vector<double> vals;
    vals.reserve(nObs);
    for (const TLeaf* leaf : fLeaves) {
        vals.push_back(leaf->GetValue());
    }
    // package that event's info into an OXO Event object, and return
    Event evt(vals);
    evt.SetObservableNames(&fObsNames);
    return evt;
}

Event
ROOTTree::GetEntry(size_t iEvent_) const{
    return Assemble(iEvent_);
}

void
ROOTTree::GatherObservableNames(){
    unsigned nObs = GetNObservables();
    fObsNames.clear();
    fLeaves.clear();
    fObsNames.reserve(nObs);
    fLeaves.reserve(nObs);
    for(unsigned i = 0; i < nObs; i++) {
        const std::string name = fTree->GetListOfBranches()->At(i)->GetName();
        fObsNames.push_back(name);
        fLeaves.push_back(fTree->GetLeaf(name.c_str()));
    }
}

void
ROOTTree::LoadBaskets(){
  if(!fTree)
    return;
  fTree->LoadBaskets();
}


void
ROOTTree::DropBaskets(){
  if(!fTree)
    return;
  fTree->DropBaskets();
}

unsigned
ROOTTree::GetNEntries() const{
    return fTree->GetEntries();
}

unsigned
ROOTTree::GetNObservables() const{
    return fTree->GetNbranches();
}
