#include <ROOTNtuple.h>
#include <Event.h>
#include <iostream>
#include <Exceptions.h>
#include <TNtuple.h>
#include <Formatter.hpp>

ROOTNtuple::ROOTNtuple(const std::string& fileName_, const std::string& treeName_){
    fROOTFile = new TFile(fileName_.c_str());

    if (fROOTFile->IsZombie()){
        delete fROOTFile;
        throw IOError("ROOTNtuple::File Does not Exist! or is Zombie " + fileName_);
    }

    fNtuple = dynamic_cast<TNtuple*>(fROOTFile -> Get(treeName_.c_str()));

    if(!fNtuple){
        delete fROOTFile;
        throw IOError(Formatter()<<"ROOTNtuple::Tree does not exist, or isn't an ntuple! tree : " << treeName_ << ", filename: "<<fileName_);
    }        
    GatherObservableNames();
}

ROOTNtuple::~ROOTNtuple(){
    if (fROOTFile)
        fROOTFile -> Close();
    delete fROOTFile;
}


std::vector<std::string>
ROOTNtuple::GetObservableNames() const{
    return fObsNames;
}

Event 
ROOTNtuple::Assemble(size_t iEvent_) const{
    if (iEvent_ >= GetNEntries())
        throw NotFoundError(Formatter() << "Exceeded end of ROOT NTuple"
                            << " \n\t(requested " << iEvent_ 
                            << " but only have " << GetNEntries() << ")");
    fROOTFile->cd();
    fNtuple -> GetEntry(iEvent_);
    float* vals = fNtuple -> GetArgs();
    Event ret(std::vector<double> (vals, vals + GetNObservables()));
    ret.SetObservableNames(&fObsNames);
    return ret;
}

Event
ROOTNtuple::GetEntry(size_t iEvent_) const{
    return Assemble(iEvent_);
}

void
ROOTNtuple::GatherObservableNames(){
    unsigned nObs = GetNObservables();
    fObsNames.clear();
    fObsNames.reserve(nObs);
    for(unsigned i = 0; i < nObs; i++)
        fObsNames.push_back((fNtuple -> GetListOfBranches()->At(i) -> GetName()));

}

void
ROOTNtuple::LoadBaskets(){
  if(!fNtuple)
    return;
  fNtuple->LoadBaskets();
}


void
ROOTNtuple::DropBaskets(){
  if(!fNtuple)
    return;
  fNtuple->DropBaskets();
}

unsigned
ROOTNtuple::GetNEntries() const{
    return fNtuple->GetEntries();
}

unsigned
ROOTNtuple::GetNObservables() const{
    return fNtuple->GetNvar();
}
