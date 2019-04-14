#include <TFile.h>
#include <TMath.h>
#include <string>
#include <TNtuple.h>
#include <TTree.h>
#include <iostream>
#include <TObject.h>
//#include <CLHEP/Random/Randomize.h>
#include <TRandom3.h>

Double_t NuSurvProb(Double_t nuE, Double_t baseline, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13){
    Double_t f_s_sqr2_theta_12 = pow(sin(2.0 * TMath::ASin(sqrt(sin_sqr_theta_12))), 2.0);
    Double_t f_s4 = pow(sin_sqr_theta_13, 2.0);
    Double_t f_c4 = pow(1.0-sin_sqr_theta_13, 2.0);
    Double_t scale = 1.267e3; // for nuE in [MeV] and baseline in [km]
    Double_t s_sqr_dm_be = pow(sin(scale * del_m_sqr_21 * baseline / nuE), 2.0);
    Double_t f_osc_prob = (f_c4 * (1.0 - f_s_sqr2_theta_12 * s_sqr_dm_be) + f_s4);
    return f_osc_prob;
}

void ntOscillate(TTree *in_tree, TTree *out_tree, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13) {

    ULong64_t n_entries = in_tree->GetEntries();
    Double_t surv_prob, mc_energy_nu, distance;
    in_tree->SetBranchAddress("mc_neutrino_energy", &mc_energy_nu);
    in_tree->SetBranchAddress("reactor_info_distance", &distance);
    TRandom3 *random_generator = new TRandom3();

    for (ULong64_t i = 0; i < n_entries; i++){
        in_tree->GetEntry(i);
        surv_prob = NuSurvProb(mc_energy_nu, distance, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13);

        //const Double_t random = CLHEP::HepUniformRand();
        Double_t random = random_generator->Uniform();

        if (surv_prob > random)
            out_tree->Fill();
    }
}

void write_file(const char* nt_in, const char* nt_out, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13) {

    TFile *f_in = new TFile(nt_in);
    TTree *in_tree = (TTree*)f_in->Get("nt");
    TTree *out_tree = in_tree->CloneTree(0);

    ntOscillate(in_tree, out_tree, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13);

    TFile *f_out = new TFile(nt_out, "RECREATE");
    out_tree->Write();
}
