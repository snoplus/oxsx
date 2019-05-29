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

void ntOscillate_pruned(TTree *in_tree, TNtuple *out_tree_ke, TNtuple *out_tree_prompt, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13, Double_t fixed_distance) {

    //
    // takes a TTree with multiple branches, oscillates using KE branch, fills two TNtuples each with a single branch
    //

    ULong64_t n_entries = in_tree->GetEntries();
    Double_t surv_prob, mc_energy_nu, ev_energy_p1, distance;
    in_tree->SetBranchAddress("mc_neutrino_energy", &mc_energy_nu);
    in_tree->SetBranchAddress("ev_fit_energy_p1", &ev_energy_p1);
    
    if (fixed_distance>0)
        distance = fixed_distance;
    else
        in_tree->SetBranchAddress("reactor_info_distance", &distance);

    TRandom3 *random_generator = new TRandom3();

    for (ULong64_t i = 0; i < n_entries; i++){
        in_tree->GetEntry(i);
        surv_prob = NuSurvProb(mc_energy_nu, distance, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13);

        //const Double_t random = CLHEP::HepUniformRand();
        Double_t random = random_generator->Uniform();

        if (surv_prob > random){
            //printf("ke:%0.5f ev:%0.5f diff:%0.5f\n", (Float_t)mc_energy_nu, (Float_t)ev_energy_p1, (Float_t)mc_energy_nu-(Float_t)ev_energy_p1);
            out_tree_ke->Fill((Float_t)mc_energy_nu);
            out_tree_prompt->Fill((Float_t)ev_energy_p1);
        }
    }
}

void ntOscillate_pruned(TTree *in_tree, TNtuple *out_tree_ke, TNtuple *out_tree_prompt, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13) {
    ntOscillate_pruned(in_tree, out_tree_ke, out_tree_prompt, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13, -9000);
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
    TFile *f_out = new TFile(nt_out, "RECREATE");
    TTree *out_tree = in_tree->CloneTree(0);

    ntOscillate(in_tree, out_tree, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13);

    f_out->cd();
    out_tree->Write();
}

void write_file_pruned(const char* nt_in, const char* nt_ke_out, const char* nt_prompt_out, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13, Double_t distance) {

    TFile *f_in = new TFile(nt_in);
    TTree *in_tree = (TTree*)f_in->Get("nt");
    
    // the ordering of things here to keep ROOT happy is exactly the following (TFile then TTree, TFile again TTree) 
    TFile *f_ke_out = new TFile(nt_ke_out, "RECREATE");
    TNtuple *out_tree_ke = new TNtuple("nt","Anti-neutrino processed tree", "mc_neutrino_energy");
    
    TFile *f_prompt_out = new TFile(nt_prompt_out, "RECREATE");
    TNtuple *out_tree_prompt = new TNtuple("nt","Anti-neutrino processed tree", "ev_fit_energy_p1");

    if (distance<0)
        ntOscillate_pruned(in_tree, out_tree_ke, out_tree_prompt, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13);
    else
        ntOscillate_pruned(in_tree, out_tree_ke, out_tree_prompt, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13, distance);

    f_ke_out->cd();
    out_tree_ke->Write();
    f_ke_out->Close();
    delete f_ke_out;
    
    f_prompt_out->cd();
    out_tree_prompt->Write();
    f_prompt_out->Close();
    delete f_prompt_out;
    
    f_in->cd();
    f_in->Close();
    delete f_in;
}

void write_file_pruned(const char* nt_in, const char* nt_ke_out, const char* nt_prompt_out, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13) {
    write_file_pruned(nt_in, nt_ke_out, nt_prompt_out, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13, -9000);
}

// alternate way of doing this:
// void write_file_pruned(const char* nt_in, const char* nt_ke_out, const char* nt_prompt_out, Double_t del_m_sqr_21, Double_t sin_sqr_theta_12, Double_t sin_sqr_theta_13) {

    // TFile *f_in = new TFile(nt_in);
    // TTree *in_tree = (TTree*)f_in->Get("nt");
    
    // TFile *f_ke_out = new TFile(nt_ke_out, "RECREATE");
    // TTree *out_tree = in_tree->CloneTree(0);
    // ntOscillate(in_tree, out_tree, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13);
    // out_tree->SetBranchStatus("*",0);
    // out_tree->SetBranchStatus("mc_neutrino_energy",1);
    // TNtuple *out_tree_ke = (TNtuple*)out_tree->CloneTree(0);
    // out_tree_ke->CopyEntries(out_tree);
    // out_tree_ke->Write();
    // f_ke_out->Close();
    // delete f_ke_out;
    
    // TFile *f_prompt_out = new TFile(nt_prompt_out, "RECREATE");
    // out_tree = in_tree->CloneTree(0);
    // ntOscillate(in_tree, out_tree, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13);
    // out_tree->SetBranchStatus("*",0);
    // out_tree->SetBranchStatus("ev_fit_energy_p1",1);
    // TNtuple *out_tree_prompt = (TNtuple*)out_tree->CloneTree(0);
    // out_tree_prompt->CopyEntries(out_tree);
    // out_tree_prompt->Write();
    // f_prompt_out->Close();
    // delete f_prompt_out;
    
    // f_in->Close();
    // delete f_in;
// }

