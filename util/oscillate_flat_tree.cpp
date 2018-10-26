#include <TFile.h>
#include <RAT/DB.hh>
#include <TH1D.h>
#include <TMath.h>
#include <string>
#include <TNtuple.h>
#include <iostream>
#include <RAT/DU/Utility.hh>
#include <TObject.h>
#include <CLHEP/Random/Randomize.h>

double NuSurvProb(double nuE, double baseline, double delmsqr21, double sinsqrtheta12, double sinsqrtheta13)
{
  double fSSqr2Theta12 = pow(sin(2.0 * TMath::ASin(sqrt(sinsqrtheta12))), 2.0);
  double fS4 = pow(sinsqrtheta13, 2.0);
  double fC4 = pow(1.0-sinsqrtheta13, 2.0);
  double scale = 1.267e3; // for nuE in [MeV] and baseline in [km]
  double sSqrDmBE = pow(sin(scale * delmsqr21 * baseline / nuE), 2.0);
  double fOscProb = (fC4 * (1.0 - fSSqr2Theta12 * sSqrDmBE) + fS4);
  return fOscProb;
}

void ntOscillate( const char* ntin, const char* ntout, double delmsqr21, double sinsqrtheta12, double sinsqrtheta13) {

  TFile *fin = new TFile(ntin);
  TTree *inT = (TTree*)fin->Get("nt");
  int nentries = inT->GetEntries();
  double survprob;
  float ParKE;
  float ReactorDistance;
  inT->SetBranchAddress("ParKE", &ParKE);
  inT->SetBranchAddress("ReactorDistance", &ReactorDistance);
  TFile *fout = new TFile(ntout,"RECREATE");
  TTree *outT = inT->CloneTree(0);
  
  for (int i = 0; i < nentries ; i++){
    inT->GetEntry(i);
    survprob = NuSurvProb(ParKE, ReactorDistance, delmsqr21, sinsqrtheta12, sinsqrtheta13);
    const double random = CLHEP::HepUniformRand();
    std::cout<<"Rand:   ";
    std::cout<<random<<std::endl;
    std::cout<<"Prob:   ";
    std::cout<<survprob<<std::endl;
    if (survprob > random){
      outT->Fill();
    }
  }
  outT->Print();
  outT->AutoSave();
  delete fin;
  delete fout;
}

int main(int argc, char* argv[])
{
  if (argc != 6){
    std::cout<<"5 arguments expected!"<<std::endl;
  }else{
    const char* rootin = argv[1];
    const char* rootout = argv[2];

    double delmsqr21 = atof(argv[3]);
    double sinsqrtheta12 = atof(argv[4]);
    double sinsqrtheta13 = atof(argv[5]);
    ntOscillate(rootin, rootout, delmsqr21, sinsqrtheta12, sinsqrtheta13);
  }
}
