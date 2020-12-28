#include <fstream>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <string>
#include <TTree.h>
#include <TVector3.h>
#include <iostream>
#include <TObject.h>
#include <math.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TNtuple.h>
#include <TChain.h>
#include <TSystemDirectory.h>

void process_cuts(const std::string filename_input, const std::string filename_output, double FV, double z_cut1, double z_cut2, double energy1_lower, double energy1_upper, double energy2_lower, double energy2_upper, size_t delTcut_lower, size_t delTcut, double delRcut, int alphan_low_high){

  /////////////
  const bool isMC = true;
  size_t muon_nhitmin = 2000;
  ////////////

  // load input file
  char *name = new char[1000];
  sprintf(name, "%s/output",filename_input.c_str());
  TChain *tree_input = new TChain("output");
  tree_input->Add(name);
  // setup output file
  TFile *file_output = new TFile(filename_output.c_str(), "RECREATE");
  TTree *tree_output = new TTree("nt","Tagged positron + neutron events");
  
  /*TH1D h_after_cut_emc_nu("h_after_cut_emc_nu", "Parent antinu KE (MeV)", 300, 0, 9);
  TH1D h_after_cut_emc_p1("h_after_cut_emc_p1", "Particle 1 KE (MeV)", 300, 0, 9);
  TH1D h_after_cut_emc_p2("h_after_cut_emc_p2", "Particle 2 KE (MeV)", 300, 0, 1);
  TH2D h2_after_cut_emc_p2_vs_p1("h2_after_cut_emc_p2_vs_p1", "Particle 2 KE vs Particle 1 KE", 10000, 0, 10, 1000, 0, 1);*/

  TH1D *h_nhit_full = new TH1D("h_nhit_full", "h_nhit_full", 1800, 0, 9000);
  TH1D *h_nhit_prompt = new TH1D("h_nhit_prompt", "h_nhit_prompt", 1000, 0, 5000);
  TH1D *h_nhit_late = new TH1D("h_nhit_late", "h_nhit_late", 240, 0, 1200);
  TH1D *h_nhit_prompt_like = new TH1D("h_nhit_prompt_like", "h_nhit_prompt_like", 1000, 0, 5000);
  TH1D *h_nhit_late_like = new TH1D("h_nhit_late_like", "h_nhit_late_like", 240, 0, 1200);
  TH1D *h_nhit_late_like_delT = new TH1D("h_nhit_late_like_delT", "h_nhit_late_like", 240, 0, 1200);
  TH1D *h_delT = new TH1D("h_delT", "h_delT", 250, delTcut_lower, delTcut);
  TH1D *h_delR = new TH1D("h_delR", "h_delR", 200, 0, 2000);
  
  TH2D *h2_delT_vs_delR = new TH2D("h2_delR_delT", "delR vs delT ", 250, delTcut_lower, delTcut, 200, 0, 2000);
  
  TH1D *h_Z_prompt = new TH1D("h_Z_prompt", "h_Z_prompt", 1200, -6000, 6000);
  TH1D *h_Z_late = new TH1D("h_Z_late", "h_Z_late", 1200, -6000, 6000);
  TH2D *h_RvsNhit_prompt = new TH2D("h_RvsNhit_prompt", "Prompt posR vs nhits", 1000, 0, 5000, 100, 0, 6000);
  TH2D *h_RvsNhit_late = new TH2D("h_RvsNhit_late", "Late posR vs nhits", 240, 0, 1200, 100, 0, 6000);
  TH2D *h_ZvsNhit_prompt = new TH2D("h_ZvsNhit_prompt", "Prompt posZ vs nhits", 1000, 0, 5000, 100, 0, 6000);
  TH2D *h_ZvsNhit_late = new TH2D("h_ZvsNhit_late", "Late posZ vs nhits", 240, 0, 1200, 100, 0, 6000);
  TH2D *h_ZvsRho_prompt = new TH2D("h_ZvsRho_prompt", "Prompt posZ vs #rho", 60, 0, 6000, 60, 0, 6000);
  TH2D *h_ZvsRho_late = new TH2D("h_ZvsRho_late", "Late posZ vs #rho", 60, 0, 6000, 60, 0, 6000);
  TH2D *h_XY_dis = new TH2D("h_XY_dis", "Tagged Events XY Distribution", 120, -6000, 6000, 120, -6000, 6000);

  TH2D *h_EdepvsZvsNhit_late = new TH2D("h_EdepvsZvsNhit_late", "Edep vs Late posZ vs nhits", 240, 0, 1200, 100, 0, 6000);
  TH2D *h_EdepQuenchvsZvsNhit_late = new TH2D("h_EdepQuenchvsZvsNhit_late", "Quench Edep vs Late posZ vs nhits", 240, 0, 1200, 100, 0, 6000);
  TH2D *h_EdepvsZvsNhit_late_count = new TH2D("h_EdepvsZvsNhit_late_count", "Edep vs Late posZ vs nhits", 240, 0, 1200, 100, 0, 6000);
  
  // properties to load from input ntuple into output ttree

  // mc parent and neutron truth positions not in standard ntuples?
  // don't need long, lat, altitude, distance here - distances for each reactor saved in a separate txt file

  //input 
  double posX, posY, posZ, SKY , edep, edepquench, Energy;
  int days, sec, nsec, ID, owl, triggerWord;
  ULong64_t clock50, DCapplied, DCflagged;
  
  //output
  Int_t ev_index_p1, ev_index_p2, ev_index;
  Double_t ev_energy_p1, ev_energy_p2, ev_energy;
  Double_t ev_delT, ev_delR;
  ULong64_t ev_clock50;
  Bool_t ev_validity;
  
  TString *reactor_core_name = 0;
  
  Double_t mc_edep_quench, mc_energy_parent1, mc_energy_parent2, mc_energy_2, mc_energy_1;
  Int_t mc_entry, mc_pdg1, mc_pdg2;

  Double_t neutron_capture_energy = 2.2;//1.857;
  Double_t e_rem = 0.784;
  Double_t annihilation = 1.022;

  // set branches input ntuple
  tree_input->SetBranchAddress("posx", &posX);
  tree_input->SetBranchAddress("posy", &posY);
  tree_input->SetBranchAddress("posz", &posZ);
  tree_input->SetBranchAddress("uTDays", &days);
  tree_input->SetBranchAddress("uTSecs", &sec);
  tree_input->SetBranchAddress("uTNSecs", &nsec);
  tree_input->SetBranchAddress("clockCount50", &clock50);
  //////////////////////////
  tree_input->SetBranchAddress("energy", &Energy);
  tree_input->SetBranchAddress("fitValid", &ev_validity);
  //////////////////////////
  tree_input->SetBranchAddress("eventID", &ID);
  tree_input->SetBranchAddress("owlnhits", &owl);
  tree_input->SetBranchAddress("skyShine", &SKY);
  tree_input->SetBranchAddress("dcApplied", &DCapplied);
  tree_input->SetBranchAddress("dcFlagged", &DCflagged);
  tree_input->SetBranchAddress("triggerWord", &triggerWord);
  // mc variables:
  tree_input->SetBranchAddress("mcIndex", &mc_entry);
  tree_input->SetBranchAddress("evIndex", &ev_index);
  tree_input->SetBranchAddress("mcEdep", &edep);
  tree_input->SetBranchAddress("mcEdepQuenched", &edepquench);

  /*tree_input->SetBranchAddress("mcIndex", &mc_entry);
    tree_input->SetBranchAddress("mcEdepQuenched", &mc_quench_i);*/
  tree_input->SetBranchAddress("parentMeta1", &reactor_core_name);
  tree_input->SetBranchAddress("parentKE1", &mc_energy_parent1);
  tree_input->SetBranchAddress("parentKE2", &mc_energy_parent2);
  tree_input->SetBranchAddress("mcke1", &mc_energy_1);
  tree_input->SetBranchAddress("mcke2", &mc_energy_2);
  tree_input->SetBranchAddress("pdg1", &mc_pdg1);
  tree_input->SetBranchAddress("pdg2", &mc_pdg2);

  // set branches output pruned ntuple

  tree_output->Branch("mc_neutrino_energy", &mc_energy_parent1);
  tree_output->Branch("mc_positron_energy", &mc_energy_1);
  tree_output->Branch("mc_neutron_energy", &mc_energy_2);
  tree_output->Branch("reactor_core_name", &reactor_core_name);
  tree_output->Branch("entry", &mc_entry);
  // values to modify
  tree_output->Branch("ev_index_p1", &ev_index_p1);
  tree_output->Branch("ev_index_p2", &ev_index_p2);
  tree_output->Branch("ev_fit_energy_p1", &ev_energy_p1);
  tree_output->Branch("ev_fit_energy_p2", &ev_energy_p2);

  ///////////////////////////////////////////////
  /// Cleaning+Coinc. Tagging (using evindex) ///
  ///////////////////////////////////////////////

  int ntagged = 0;
  ULong64_t nsimmed = 1;
  std::vector<int> tagged214;

  int late_candidates = 0;

  std::cout<<"n_entries: "<<tree_input->GetEntries()<<std::endl;
  for (int i = 0; i < tree_input->GetEntries(); i++) {
    if (i == 0) continue;

    tree_input->GetEntry(i);
    reactor_core_name = 0;
    double x1, y1, z1, r1, rho1, sky1, edep1, edepquench1, energy1;
    long long time1;
    int date1, sec1, nsec1, id1, owl1,triggerWord1, ev_index1;
    ULong64_t clock1;

    x1 = posX;
    y1 = posY;
    rho1 = pow(x1*x1 + y1*y1, .5);
    z1 = posZ - 108;
    r1 = pow((x1*x1 + y1*y1 + z1*z1), .5);
    energy1 = Energy;
    time1 = days * 24 * 3600 * pow(10, 9) + sec * pow(10, 9) + nsec;
    clock1 = clock50;
    id1 = ID;
    owl1 = owl;
    sky1 = SKY;
    triggerWord1 = triggerWord;
      
    date1 = days;
    sec1 = sec;
    nsec1 = nsec;

    edep1 = edep;
    edepquench1 = edepquench;
    ev_index1 = ev_index;

    // muon tag
    bool muon_found = false;

    //if (((DCapplied & 0x210000000242) & DCflagged ) != (DCapplied & 0x210000000242)) continue;
    //if (((DCapplied & 0x210000003FF6) & DCflagged ) != (DCapplied & 0x210000003FF6)) continue;

    /*if (!isMC)
    if (((DCapplied & 0x21000000FFF6) & DCflagged ) != (DCapplied & 0x21000000FFF6)) continue;   //Matt Mask
    */

    //  for MC events, neutron will have evindex >= 1
    if (isMC && (ev_index1 <= 0))//{
      nsimmed += 1;

    if (isMC && (ev_index1 == 0)) continue;

    if (!ev_validity) continue;

    //if (sky1 < 1) continue;

    for (int d = 1; d < 9; d++) {	
      if (nsec1 % 10 == 0) nsec1 = nsec1 / 10;
    }
    
    if (r1 > FV || z1 > z_cut1 || z1 < z_cut2) continue;
      
    if (energy1 > energy1_lower)
      h_nhit_full->Fill(energy1);
      
    if (energy1 > energy2_upper || energy1 < energy2_lower) continue;
    
    // counting prompt event candidates
    late_candidates += 1;
    //std::cout<<" -----------> late candidate"<<std::endl;

    bool pair = false;

    for (int ii = 1; ii < (i + 1); ii++) {
	
      if (pair) break;

      bool tagged = false;
	
      for (int i214 = 0; i214 < tagged214.size(); i214++) {
        if ((i - ii) == tagged214[i214]) {
          tagged = true;
          break;
        } 
      }
	
      if (tagged == true) continue;

      tree_input->GetEntry(i - ii);

      if (!ev_validity) continue;

      //if (((DCapplied & 0x210000003FF6) & DCflagged ) != (DCapplied & 0x210000003FF6)) continue;
      if (!isMC)
        if (((DCapplied & 0x21000000FFF6) & DCflagged ) != (DCapplied & 0x21000000FFF6)) continue;   //Matt Mask

      double x2, y2, z2, r2, rho2, sky2, energy2;
      long long time2;
      int id2, owl2, triggerWord2, mc_entry2, ev_index2;
      ULong64_t clock2;
      //TString reactor_core_name2;

      x2 = posX;
      y2 = posY;
      rho2 = pow(x2*x2 + y2*y2, .5);
      z2 = posZ - 108;
      r2 = pow((x2*x2 + y2*y2 + z2*z2), .5);
      energy2 = Energy;
      time2 = days * 24 * 3600 * pow(10, 9) + sec * pow(10, 9) + nsec;
      clock2 = clock50;
      id2 = ID;
      owl2 = owl;
      sky2 = SKY;
      triggerWord2 = triggerWord;

      ev_index2 = ev_index;
      mc_entry2 = mc_entry;

      //reactor_core_name2 = reactor_core_name;
      
      ULong64_t delT = (clock1 - clock2) * 20;

      double delR = pow( ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)), 0.5);

      if (delT > delTcut) break;

      if (energy2 < energy1_lower || energy2 > energy1_upper) continue;
      
      //if (sky2 < 1) continue;

      if (r2 > FV || z2 > z_cut1 || z2 < z_cut2) continue;

      if (delR > delRcut) continue;
	
      if (delT < delTcut_lower) continue;

      ////////////////////////////
      //if (delT > 1000000) break;
      ///////////////////////////
      
      //////////
      // remove gamma events from mc:
      //if (energy2 < 3.4)//5.2)
      //std::cout<<"energy1: "<<energy2<<" ke1: "<<mc_energy_1<<" ke2: "<<mc_energy_2<<" pdg1: "<<mc_pdg1<<" pdg2: "<<mc_pdg2<<std::endl;
      bool high_low_continue = false;
      if (alphan_low_high == 1){
        if ((6.1 < mc_energy_1 < 6.2 and mc_pdg1 == 22) or\
            ((mc_energy_1+mc_energy_2 > 5) and mc_pdg1 == -11 and mc_pdg2 == 11) or
            ((mc_energy_1+mc_energy_2 > 2.513) and mc_pdg1 == -11 and mc_pdg2 == 1000020040)){
          high_low_continue = true;
          //std::cout<<"----->skipping 6.1MeV gammas (pair)  ev_energy: "<<energy2<<std::endl;
        }
      }else if (alphan_low_high == 2){
        high_low_continue = true;
        if ((6.1 < mc_energy_1 < 6.2 and mc_pdg1 == 22) or\
            ((mc_energy_1+mc_energy_2 > 5) and mc_pdg1 == -11 and mc_pdg2 == 11) or
            ((mc_energy_1+mc_energy_2 > 2.513) and mc_pdg1 == -11 and mc_pdg2 == 1000020040)){
          high_low_continue = false;
          //std::cout<<"keeping 6.1MeV gammas (pair)"<<std::endl;
        }
      }
      
      if (high_low_continue) continue;
      ////////
      
      ntagged += 1;
      //std::cout<<"-------------------------------------------------------------\n"<< \
          //"Tagged pair "<<" prompt_i = "<<(i-ii)<<" nhit: "<<energy2<<"  late_i = "<<i<<"nhit: "<<energy1<<"\n-------------------------------------------------------------"<<std::endl;

      pair = true;
      tagged214.push_back(i);
      tagged214.push_back(i-ii);

      h_nhit_prompt->Fill(energy2);
      h_nhit_late->Fill(energy1);
      h_delT->Fill(delT);
      h_delR->Fill(delR);
      h2_delT_vs_delR->Fill(delT,delR);
      h_RvsNhit_prompt->Fill(energy2, r2);
      h_ZvsNhit_prompt->Fill(energy2, z2);
      h_ZvsRho_prompt->Fill(rho2, z2);
      h_RvsNhit_late->Fill(energy1, r1);
      h_ZvsNhit_late->Fill(energy1, z1);
      h_ZvsRho_late->Fill(rho1, z1);
      h_XY_dis->Fill(x1, y1);
      h_XY_dis->Fill(x2, y2);
      h_Z_prompt->Fill(z2);
      h_Z_late->Fill(z1);
      
      /*h_after_cut_emc_nu->Fill(mc_energy_nu);
      h_after_cut_emc_p1->Fill(mc_energy_pe);
      h_after_cut_emc_p2->Fill(mc_energy_n);*/

      // for output
      mc_entry = mc_entry2;
      ev_index_p1 = ev_index2;
      ev_index_p2 = ev_index1;
      ev_energy_p1 = energy2;
      ev_energy_p2 = energy1;

      tree_output->Fill();

      if (isMC){
        h_EdepvsZvsNhit_late->Fill(energy1,z1,edep1);
        h_EdepQuenchvsZvsNhit_late->Fill(energy1,z1,edep1);
        h_EdepvsZvsNhit_late_count->Fill(energy1,z1);
      }
    }
  }

    // save file
    //file_input->Close();
    file_output->cd();

    // write objects...
    h_nhit_full->Write();
    h_nhit_prompt->Write();
    h_nhit_late->Write();
    h_delT->Write();
    h_delR->Write();
    h2_delT_vs_delR->Write();
    h_Z_prompt->Write();
    h_Z_late->Write();
    h_RvsNhit_prompt->Write();
    h_ZvsNhit_prompt->Write();
    h_ZvsRho_prompt->Write();
    h_RvsNhit_late->Write();
    h_ZvsNhit_late->Write();
    h_ZvsRho_late->Write();
    h_XY_dis->Write();    
    
    //write csv output file with event numbers
    sprintf(name, "%s.csv",filename_output.c_str());

    std::ofstream fOut(name);
    fOut<<"filename_input,filename_output,Rmax,Zmin,Zmax,energy_ep_min,energy_ep_max,energy_n_min,energy_n_max,deltaTmin,deltaTmax,deltaRmax,initial_entries,final_entries,finished\n";
    fOut<<filename_input<<", "<<filename_output<<", "<<FV<<", "<<z_cut1<<", "<<z_cut2<<", "<<energy1_lower<<", "<<energy1_upper<<", "<<energy2_lower<<", "<<energy2_upper<<", "<<delTcut_lower<<", "<<delTcut<<", "<<delRcut<<", "<<nsimmed<<", "<<ntagged<<", 1\n";
    fOut.close();

    std::cout<<"fin: "<<filename_input<<" \n fout: "<<filename_output<<"\n Rmax "<<FV<<"\n Zmin "<<z_cut2<<"\n Zmax "<<z_cut1<<"\n e1min "<<energy1_lower<<"\n e1max "<<energy1_upper<<"\n e2min "<<energy2_lower<<"\n e2max "<<energy2_upper<<"\n delT "<<delTcut_lower<<"\n delTmax "<<delTcut<<"\n delR "<<delRcut<<"\n nsimmed: "<<nsimmed<<"\n ntagged: "<<ntagged<<"\n nlatecandidates: "<<late_candidates<<std::endl;

    tree_output->AutoSave();
    file_output->Close();
    delete tree_input;
    delete file_output;
}
