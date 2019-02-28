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

Double_t NuSurvProb(Double_t nuE, Double_t baseline, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13)
{
  Double_t f_s_sqr2_theta_12 = pow(sin(2.0 * TMath::ASin(sqrt(sin_sqr_theta_12))), 2.0);
  Double_t f_s4 = pow(sin_sqr_theta_13, 2.0);
  Double_t f_c4 = pow(1.0-sin_sqr_theta_13, 2.0);
  Double_t scale = 1.267e3; // for nuE in [MeV] and baseline in [km]
  Double_t s_sqr_dm_be = pow(sin(scale * del_m_sqr_21 * baseline / nuE), 2.0);
  Double_t f_osc_prob = (f_c4 * (1.0 - f_s_sqr2_theta_12 * s_sqr_dm_be) + f_s4);
  return f_osc_prob;
}

void ntOscillate(const char* nt_in, const char* nt_out, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13) {

  TFile *f_in = new TFile(nt_in);
  TTree *in_tree = (TTree*)f_in->Get("nt");
  ULong64_t n_entries = in_tree->GetEntries();
  Double_t surv_prob, mc_energy_nu, distance;
  in_tree->SetBranchAddress("mc_neutrino_energy", &mc_energy_nu);
  in_tree->SetBranchAddress("reactor_info_distance", &distance);
  TFile *f_out = new TFile(nt_out,"RECREATE");
  TTree *out_tree = in_tree->CloneTree(0);
  
  for (ULong64_t i = 0; i < n_entries; i++){
    in_tree->GetEntry(i);
    surv_prob = NuSurvProb(mc_energy_nu, distance, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13);
    //std::cout<<mc_energy_nu<<" "<<distance<<" "<<surv_prob<<" "<<del_m_sqr_21<<" "<<sin_sqr_theta_12<<" "<<sin_sqr_theta_13<<std::endl;
    const Double_t random = CLHEP::HepUniformRand();
    //std::cout<<"Rand:   ";
    //std::cout<<random<<std::endl;
    //std::cout<<"Prob:   ";
    //std::cout<<surv_prob<<std::endl;
    if (surv_prob > random){
      out_tree->Fill();
    }
  }
  //out_tree->Print();
  out_tree->AutoSave();
  delete f_in;
  delete f_out;
}

Int_t main(int argc, char* argv[]){
  if (argc != 6){
    std::cout<<"5 arguments expected!"<<std::endl;
  }else{
    const char* root_in = argv[1];
    const char* root_out = argv[2];

    double del_m_sqr_21 = atof(argv[3]);
    double sin_sqr_theta_12 = atof(argv[4]);
    double sin_sqr_theta_13 = atof(argv[5]);
    ntOscillate(root_in, root_out, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13);
  }
}
