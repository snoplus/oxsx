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

//needed for KL P1
double correction(double x){
  //double correc = 8.83411e-02 + (9.02287e-03)*x;
  double correc = 0.;
  return correc;
}


//AVERAGE
double Average(std::vector<double> v){      double sum=0;
  for(int i=0;i<abs(v.size());i++)
    sum+=v[i];
  return sum/(double)v.size();
}

//DEVIATION
double Uncert(std::vector<double> v, double ave){
  double E=0;
  for(int i=0;i<abs(v.size());i++)
    E+=(v[i] - ave)*(v[i] - ave);
  int sizemin1size = (v.size()-1)*v.size();
  return (sqrt(E/(double)sizemin1size));
}

double myline (double *x, double *par){
  double pdf = par[0] + par[1]*x[0];
  return pdf;
}

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

  TH1D *h_energy_full = new TH1D("h_energy_full", "h_energy_full", 1800, 0, 9000);
  TH1D *h_energy_prompt = new TH1D("h_energy_prompt", "h_energy_prompt", 1000, 0, 5000);
  TH1D *h_energy_late = new TH1D("h_energy_late", "h_energy_late", 240, 0, 1200);
  TH1D *h_energy_prompt_like = new TH1D("h_energy_prompt_like", "h_energy_prompt_like", 1000, 0, 5000);
  TH1D *h_energy_late_like = new TH1D("h_energy_late_like", "h_energy_late_like", 240, 0, 1200);
  TH1D *h_energy_late_like_delT = new TH1D("h_energy_late_like_delT", "h_energy_late_like", 240, 0, 1200);
  TH1D *h_delT = new TH1D("h_delT", "h_delT", 250, delTcut_lower, delTcut);
  TH1D *h_delR = new TH1D("h_delR", "h_delR", 200, 0, 2000);
  
  TH2D *h2_delT_vs_delR = new TH2D("h2_delR_delT", "delR vs delT ", 250, delTcut_lower, delTcut, 200, 0, 2000);
  
  TH1D *h_Z_prompt = new TH1D("h_Z_prompt", "h_Z_prompt", 1200, -6000, 6000);
  TH1D *h_Z_late = new TH1D("h_Z_late", "h_Z_late", 1200, -6000, 6000);
  TH2D *h_RvsEnergy_prompt = new TH2D("h_RvsEnergy_prompt", "Prompt posR vs energys", 1000, 0, 5000, 100, 0, 6000);
  TH2D *h_RvsEnergy_late = new TH2D("h_RvsEnergy_late", "Late posR vs energys", 240, 0, 1200, 100, 0, 6000);
  TH2D *h_ZvsEnergy_prompt = new TH2D("h_ZvsEnergy_prompt", "Prompt posZ vs energys", 1000, 0, 5000, 100, 0, 6000);
  TH2D *h_ZvsEnergy_late = new TH2D("h_ZvsEnergy_late", "Late posZ vs energys", 240, 0, 1200, 100, 0, 6000);
  TH2D *h_ZvsRho_prompt = new TH2D("h_ZvsRho_prompt", "Prompt posZ vs #rho", 60, 0, 6000, 60, 0, 6000);
  TH2D *h_ZvsRho_late = new TH2D("h_ZvsRho_late", "Late posZ vs #rho", 60, 0, 6000, 60, 0, 6000);
  TH2D *h_XY_dis = new TH2D("h_XY_dis", "Tagged Events XY Distribution", 120, -6000, 6000, 120, -6000, 6000);

  TH2D *h_EdepvsZvsEnergy_late = new TH2D("h_EdepvsZvsEnergy_late", "Edep vs Late posZ vs energys", 240, 0, 1200, 100, 0, 6000);
  TH2D *h_EdepQuenchvsZvsEnergy_late = new TH2D("h_EdepQuenchvsZvsEnergy_late", "Quench Edep vs Late posZ vs energys", 240, 0, 1200, 100, 0, 6000);
  TH2D *h_EdepvsZvsEnergy_late_count = new TH2D("h_EdepvsZvsEnergy_late_count", "Edep vs Late posZ vs energys", 240, 0, 1200, 100, 0, 6000);
  
  // properties to load from input ntuple into output ttree

  //input 
  double posX, posY, posZ, SKY , energy, edep, edepquench;
  int nhits, days, sec, nsec, ID, owl, triggerWord, evinde;
  ULong64_t clock50, DCapplied, DCflagged;
  
  //output
  Int_t ev_index, ev_index_p1, ev_index_p2;
  Double_t ev_energy_p1, ev_energy_p2;
  Double_t ev_delT, ev_delR;
  ULong64_t ev_clock50;
  Bool_t ev_validity;
  //Double_t ev_pos_x, ev_pos_y, ev_pos_z;
  
  TString *reactor_core_name = 0;
  
  Double_t mc_edep_quench, mc_energy_parent1, mc_energy_parent2, mc_energy_2, mc_energy_1;
  Int_t mc_entry, mc_pdg1, mc_pdg2;

  Double_t neutron_capture_energy = 2.2;//1.857;
  Double_t e_rem = 0.784;
  Double_t annihilation = 1.022;

  //Iwan
  /*ULong64_t numsimmed = 1;
  ULong64_t numtagged = 0;
  ULong64_t ev01 = 0;
  ULong64_t ev02 = 0;
  ULong64_t ev03 = 0;
  ULong64_t badev1 = 0;
  Double_t totalbadevdeltaT = 0.;
  
  Double_t mc_time_nanoseconds;
  */
  
  // set branches input ntuple
  tree_input->SetBranchAddress("posx", &posX);
  tree_input->SetBranchAddress("posy", &posY);
  tree_input->SetBranchAddress("posz", &posZ);
  tree_input->SetBranchAddress("uTDays", &days);
  tree_input->SetBranchAddress("uTSecs", &sec);
  tree_input->SetBranchAddress("uTNSecs", &nsec);
  tree_input->SetBranchAddress("clockCount50", &clock50);
  //////////////////////////
  tree_input->SetBranchAddress("nhitsCleaned", &nhits);
  tree_input->SetBranchAddress("energy", &energy);
  //////////////////////////
  tree_input->SetBranchAddress("fitValid", &ev_validity);
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
    energy1 = energy;
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
    /*if (!isMC){
      if (run_number_int >= 258979 && run_number_int <= 259946){
        if (nhits1 > muon_nhitmin){
          muon_found = true;
          TVector3 position = TVector3(x1,y1,z1);
          muon_pos.push_back(position);
          muon_clock50.push_back(clock1);
          muon_gtid.push_back(id1);
          muon_nhit.push_back(nhits1);
	  
          std::cout<<"\n MUON FOUND (Crates off) "<<std::endl;

          ULong64_t muon_tag_time = clock1 * 20;
          int date_muon = muon_tag_time / 86400000000000;
          int sec_muon = (muon_tag_time - date_muon * 86400000000000) / 1000000000;
          int nsec_muon = muon_tag_time - date_muon * 86400000000000 - sec_muon * 1000000000;

          std::stringstream muon_tag_time_strm;
          muon_tag_time_strm<<date_muon<<" "<<sec_muon<<"."<<nsec_muon;
          std::cout<<" pos: ("<<position.X()<<", "<<position.Y()<<", "<<position.Z()<<")  clock50_sec: "<<muon_tag_time_strm.str()<<"  gtid: "<<id1<<"  nhit: "<<nhits1<<"\n"<<std::endl;
        }
      }else{
        if (!(((DCapplied & 0x80) & DCflagged ) == (DCapplied & 0x80))){
          muon_found = true;
          TVector3 position = TVector3(x1,y1,z1);
          muon_pos.push_back(position);
          muon_clock50.push_back(clock1);
          muon_gtid.push_back(id1);
          muon_nhit.push_back(nhits1);

          std::cout<<"\n MUON FOUND"<<std::endl;

          ULong64_t muon_tag_time = clock1 * 20;
          int date_muon = muon_tag_time / 86400000000000;
          int sec_muon = (muon_tag_time - date_muon * 86400000000000) / 1000000000;
          int nsec_muon = muon_tag_time - date_muon * 86400000000000 - sec_muon * 1000000000;

          std::stringstream muon_tag_time_strm;
          muon_tag_time_strm<<date_muon<<" "<<sec_muon<<"."<<nsec_muon;
          std::cout<<" pos: ("<<position.X()<<", "<<position.Y()<<", "<<position.Z()<<")  clock50_sec: "<<muon_tag_time_strm.str()<<"  gtid: "<<id1<<"  nhit: "<<nhits1<<"\n"<<std::endl;
        }
      }
      }*/

    //if (((DCapplied & 0x210000000242) & DCflagged ) != (DCapplied & 0x210000000242)) continue;
    //if (((DCapplied & 0x210000003FF6) & DCflagged ) != (DCapplied & 0x210000003FF6)) continue;

    /*if (!isMC)
    if (((DCapplied & 0x21000000FFF6) & DCflagged ) != (DCapplied & 0x21000000FFF6)) continue;   //Matt Mask
    */

    //  for MC events, neutron will have evindex >= 1
    if (isMC && (ev_index1 <= 0))//{
      nsimmed += 1;
    /*std::cout<<"\ni: "<<i<<" evindex: "<<ev_index1<<" energy: "<<energy1<<" ("<<energy1_lower<<","<<energy1_upper<<") "<<" r1: "<<r1<<" z1: "<<z1<<" ("<<z_cut2<<","<<z_cut1<<") "<<" sky: "<<sky1<<" "<<std::endl;
    }else{
      std::cout<<"i: "<<i<<" evindex: "<<ev_index1<<" energy: "<<energy1<<" ("<<energy2_lower<<","<<energy2_upper<<") "<<" r1: "<<r1<<" z1: "<<z1<<" ("<<z_cut2<<","<<z_cut1<<") "<<" sky: "<<sky1<<" "<<std::endl;
    }*/

    if (isMC && (ev_index1 == 0)) continue;

    if (!ev_validity) continue;

    if (sky1 < 1) continue;

    for (int d = 1; d < 9; d++) {	
      if (nsec1 % 10 == 0) nsec1 = nsec1 / 10;
    }

    //std::cout<<"i: "<<i<<" evindex: "<<ev_index<<", "<<ev_index1<<" energy: "<<energy<<" r1: "<<r1<<" z1: "<<z1<<"  "<<FV<<"  "<<z_cut1<<"  "<<z_cut2<<std::endl;

    if (r1 > FV || z1 > z_cut1 || z1 < z_cut2) continue;
      
    if (energy1 > energy1_lower)
      h_energy_full->Fill(energy1);
      
    /*// counting prompt event candidates
    if (z1 > Nhit_Zcut) {
      if (nhits1 < nhit1_upper && nhits1 > nhit1_lower){
        prompt_candidates += 1;
        prompt_candidate_entries.push_back(i);
        h_nhit_prompt_like->Fill(nhits1);
      }
    } else {
      if (nhits1 < nhit1_upper_2nd && nhits1 > nhit1_lower_2nd){
        prompt_candidates += 1;
        prompt_candidate_entries.push_back(i);
        h_nhit_prompt_like->Fill(nhits1);
      }
    }*/

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

      //if (((DCapplied & 0x210000003FF6) & DCflagged ) != (DCapplied & 0x210000003FF6)) continue;
      if (!isMC)
        if (((DCapplied & 0x21000000FFF6) & DCflagged ) != (DCapplied & 0x21000000FFF6)) continue;   //Matt Mask
      
      if (!ev_validity) continue;


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
      energy2 = energy;
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
      
      if (sky2 < 1) continue;

      if (r2 > FV || z2 > z_cut1 || z2 < z_cut2) continue;

      if (delR > delRcut) continue;
	
      if (delT < delTcut_lower) continue;

      ////////////////////////////
      //if (delT > 1000000) break;
      ///////////////////////////
      
      //////////
      // remove gamma events from mc:
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

      /*tag_prompt_pos.push_back(TVector3(x2,y2,z2));
      tag_late_pos.push_back(TVector3(x1,y1,z1));

      tag_prompt_time.push_back(clock2*20.);
      tag_late_time.push_back(clock1*20.);


      ULong64_t tag_time = clock1 * 20;
      int date_1 = tag_time / 86400000000000;
      int sec_1 = (tag_time - date_1 * 86400000000000) / 1000000000;
      int nsec_1 = tag_time - date_1 * 86400000000000 - sec_1 * 1000000000;

      std::ostringstream BiPo_date_str, BiPo_sec_str, BiPo_nsec_str;
      BiPo_date_str << date_1;
      BiPo_sec_str << sec_1;
      BiPo_nsec_str << nsec_1;
      std::string event_time = BiPo_date_str.str() + " " + BiPo_sec_str.str() + "." + BiPo_nsec_str.str() + "\n";
      outfile << event_time.c_str();

      ULong64_t tag_time_2 = clock2 * 20;
      int date_2 = tag_time_2 / 86400000000000;
      int sec_2 = (tag_time_2 - date_2 * 86400000000000) / 1000000000;
      int nsec_2 = tag_time_2 - date_2 * 86400000000000 - sec_2 * 1000000000;

      std::stringstream tag_time_2_strm;
      tag_time_2_strm<<date_2<<" "<<sec_2<<"."<<nsec_2;

      std::stringstream tag_time_1_strm;
      tag_time_1_strm<<date_1<<" "<<sec_1<<"."<<nsec_1;

      std::cout<<"\n-------------------------------------------------------------\n"<< \
      "Tagged pair, event time =  "<<event_time<<" prompt_i = "<<(i-ii)<<" late_i = "<<i<<"  gtid1 = "<<id2<<" gtid2 = "<<id1<<" (clock50)day_sec_1 = "<<tag_time_2_strm.str()<<" (clock50)day_sec_2 = "<<tag_time_1_strm.str()<<" trigword1 = "<<triggerWord2<<" trigword2 = "<<triggerWord1 \
	       <<"\n-------------------------------------------------------------"<<std::endl;
      

      std::ostringstream Bi_id_str, Po_id_str;
      Bi_id_str << id2;
      Po_id_str << id1;
      std::string event_UTID = Bi_id_str.str() + " " + Po_id_str.str() + " ";
      outfile3 << event_UTID.c_str();
      */
      pair = true;
      tagged214.push_back(i);
      tagged214.push_back(i-ii);

      h_energy_prompt->Fill(energy2);
      h_energy_late->Fill(energy1);
      h_delT->Fill(delT);
      h_delR->Fill(delR);
      h2_delT_vs_delR->Fill(delT,delR);
      h_RvsEnergy_prompt->Fill(energy2, r2);
      h_ZvsEnergy_prompt->Fill(energy2, z2);
      h_ZvsRho_prompt->Fill(rho2, z2);
      h_RvsEnergy_late->Fill(energy1, r1);
      h_ZvsEnergy_late->Fill(energy1, z1);
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
        h_EdepvsZvsEnergy_late->Fill(energy1,z1,edep1);
        h_EdepQuenchvsZvsEnergy_late->Fill(energy1,z1,edep1);
        h_EdepvsZvsEnergy_late_count->Fill(energy1,z1);
      }
    }
  }
  /*
    // look at how reconstruction compares w/ MC truth
    //  dele = positronKE + 1.022MeV - reconstructedEPrompt
    std::vector<Double_t> dele_e1xvec;
    std::vector<Double_t> dele_e1yvec;
    std::vector<Double_t> dele_e1yerrvec;
    std::vector<Double_t> dele_poskexvec;
    std::vector<Double_t> dele_poskeyvec;
    std::vector<Double_t> dele_poskeyerrvec;
    std::vector<Double_t> dele_poskexvec01;
    std::vector<Double_t> dele_poskeyvec01;
    std::vector<Double_t> dele_poskeyerrvec01;
    std::vector<Double_t> dele_poskexvec02;
    std::vector<Double_t> dele_poskeyvec02;
    std::vector<Double_t> dele_poskeyerrvec02;
    std::vector<Double_t> dele_e1xerrvec;
    std::vector<Double_t> dele_poskexerrvec;
    std::vector<Double_t> dele_poskexerrvec01;
    std::vector<Double_t> dele_poskexerrvec02;

    for (ULong64_t i = 0; i < nxbins; i++){
      std::vector<Double_t> slicecontente1;
      std::vector<Double_t> slicecontent;
      std::vector<Double_t> slicecontent01;
      std::vector<Double_t> slicecontent02;
      for (ULong64_t j = 0; j < nybins; j++){
        ULong64_t pointcounte1 = delE_efit_prompt.GetBinContent(i+1,j+1);
        ULong64_t pointcount = delE_emc_p1.GetBinContent(i+1,j+1);
        ULong64_t pointcount01 = delE_emc_p1_0_1.GetBinContent(i+1,j+1);
        ULong64_t pointcount02 = delE_emc_p1_0_2.GetBinContent(i+1,j+1);
        if (pointcounte1 >= 1){
          for (ULong64_t k = 0; k < pointcounte1; k++)
            slicecontente1.push_back(((TAxis*)delE_efit_prompt.GetYaxis())->GetBinCenter(j+1));
        }
        if (pointcount >= 1){
          for (ULong64_t k = 0; k < pointcount; k++)
            slicecontent.push_back(((TAxis*)delE_emc_p1.GetYaxis())->GetBinCenter(j+1));
        }
        if (pointcount01 >= 1){
          for (ULong64_t k = 0; k < pointcount01; k++)
            slicecontent01.push_back(((TAxis*)delE_emc_p1_0_1.GetYaxis())->GetBinCenter(j+1));
        }
        if (pointcount02 >= 1){
          for (ULong64_t k = 0; k < pointcount02; k++)
            slicecontent02.push_back(((TAxis*)delE_emc_p1_0_2.GetYaxis())->GetBinCenter(j+1));
        }
      }
      ULong64_t minnum = 10; //minimum sample size used to calculate a mean and uncertainty
      if (abs(slicecontente1.size()) >= minnum){
        double mean = Average(slicecontente1);
        double uncert = Uncert(slicecontente1,mean);
        //std::cout<<"\n E1: "<<((TAxis*)delE_E1.GetXaxis())->GetBinCenter(i+1)<<" deleavg: "<<mean<<" uncert: "<<uncert<<" num in slice: "<<slicecontente1.size()<<std::endl;
        dele_e1xvec.push_back(((TAxis*)delE_efit_prompt.GetXaxis())->GetBinCenter(i+1));
        dele_e1yvec.push_back(mean);
        dele_e1yerrvec.push_back(uncert);
        dele_e1xerrvec.push_back(0.);
      }
      if (abs(slicecontent.size()) >= minnum){
        double mean = Average(slicecontent);
        double uncert = Uncert(slicecontent,mean);
        //std::cout<<" poske: "<<((TAxis*)delE_posKE.GetXaxis())->GetBinCenter(i+1)<<" deleavg: "<<mean<<" uncert: "<<uncert<<" num in slice: "<<slicecontent.size()<<std::endl;
        dele_poskexvec.push_back(((TAxis*)delE_emc_p1.GetXaxis())->GetBinCenter(i+1));
        dele_poskeyvec.push_back(mean);
        dele_poskeyerrvec.push_back(uncert);
        dele_poskexerrvec.push_back(0.);
      }
      if (abs(slicecontent01.size()) >= minnum){
        double mean01 = Average(slicecontent01);
        double uncert01 = Uncert(slicecontent01,mean01);
        //std::cout<<"0_1: poske: "<<((TAxis*)delE_posKE01.GetXaxis())->GetBinCenter(i+1)<<" deleavg: "<<mean01<<" uncert: "<<uncert01<<" num in slice: "<<slicecontent01.size()<<std::endl;
        dele_poskexvec01.push_back(((TAxis*)delE_emc_p1_0_1.GetXaxis())->GetBinCenter(i+1));
        dele_poskeyvec01.push_back(mean01);
        dele_poskeyerrvec01.push_back(uncert01);
        dele_poskexerrvec01.push_back(0.);
      }
      if (abs(slicecontent02.size()) >= minnum){
        double mean02 = Average(slicecontent02);
        double uncert02 = Uncert(slicecontent02,mean02);
        //std::cout<<"0_2: poske: "<<((TAxis*)delE_posKE02.GetXaxis())->GetBinCenter(i+1)<<" deleavg: "<<mean02<<" uncert: "<<uncert02<<" num in slice: "<<slicecontent02.size()<<std::endl;
        dele_poskexvec02.push_back(((TAxis*)delE_emc_p1_0_2.GetXaxis())->GetBinCenter(i+1));
        dele_poskeyvec02.push_back(mean02);
        dele_poskeyerrvec02.push_back(uncert02);
        dele_poskexerrvec02.push_back(0.);
      }
    }
    // Plotting (and fitting) DelE vs reconstructed EPrompt
    TCanvas* c_delekeplus_e1 = new TCanvas("c_delekeplus_e1", "positronKE + 1.022MeV - EPrompt vs EPrompt");
    TGraphErrors *gr_dele_e1 = new TGraphErrors(dele_e1xvec.size(),&dele_e1xvec[0],&dele_e1yvec[0],&dele_e1xerrvec[0],&dele_e1yerrvec[0]);

    TF1 *fitfunc = new TF1("fitfunc", myline, e1min, e1max, 2);
    fitfunc->SetParName(0,"C");
    fitfunc->SetParName(1,"m");

    double Cmax = 0.6;
    double mmax = 10;

    TRandom3 *r1 = new TRandom3();
    r1->SetSeed(0);

    fitfunc->SetParameters(r1->Rndm()*(Cmax),(r1->Rndm())*mmax);
    fitfunc->SetParLimits(0, 0., Cmax);
    fitfunc->SetParLimits(1, 0., mmax);
    gr_dele_e1->Fit(fitfunc);
    Double_t par[2];
    fitfunc->GetParameters(par);
    gStyle->SetOptFit(1);
    fitfunc->Draw("same");

    gr_dele_e1->SetTitle("(positronKE + 1.022MeV - EPrompt) vs EPrompt");
    gr_dele_e1->SetMarkerColor(kRed);
    gr_dele_e1->SetLineColor(kRed);
    gr_dele_e1->SetMarkerStyle(21);
    gr_dele_e1->Draw("AP");

    // Plotting DelE vs positron KE
    TCanvas* c_delekeplus_poske = new TCanvas("c_delekeplus_poske", "positronKE + 1.022MeV - EPrompt vs positronKE");
    TGraphErrors *gr_dele_poske = new TGraphErrors(dele_poskexvec.size(),&dele_poskexvec[0],&dele_poskeyvec[0],&dele_poskexerrvec[0],&dele_poskeyerrvec[0]);

    gr_dele_poske->SetTitle("(positronKE + 1.022MeV - EPrompt) vs positronKE");
    gr_dele_poske->SetMarkerColor(kRed);
    gr_dele_poske->SetLineColor(kRed);
    gr_dele_poske->SetMarkerStyle(21);
    gr_dele_poske->Draw("AP");

    // Plotting DelE vs positron KE (ev index 0_1 and 0_2)
    TCanvas* c_delekeplus_poske_1v2 = new TCanvas("c_delekeplus_e11v2", "positronKE + 1.022MeV - EPrompt vs EPrompt (0_1 vs 0_2)");
    TGraphErrors *gr_dele_poske01 = new TGraphErrors(dele_poskexvec01.size(),&dele_poskexvec01[0],&dele_poskeyvec01[0],&dele_poskexerrvec01[0],&dele_poskeyerrvec01[0]);
    gr_dele_poske01->SetTitle("(positronKE + 1.022MeV - EPrompt) vs  positronKE   evindex= 0_1(red) vs 0_2(blue)" );
    gr_dele_poske01->SetMarkerColor(kRed);
    gr_dele_poske01->SetLineColor(kRed);
    gr_dele_poske01->SetMarkerStyle(21);
    gr_dele_poske01->Draw("AP");
    TGraphErrors *gr_dele_poske02 = new TGraphErrors(dele_poskexvec02.size(),&dele_poskexvec02[0],&dele_poskeyvec02[0],&dele_poskexerrvec02[0],&dele_poskeyerrvec02[0]);
    gr_dele_poske02->SetMarkerColor(kBlue);
    gr_dele_poske02->SetLineColor(kBlue);
    gr_dele_poske02->SetMarkerStyle(21);
    gr_dele_poske02->Draw("P");

    //summary stats
    std::cout<<"\n \n Summary: "<<std::endl;
    std::cout<<"\n number of antinus simmed: "<<numsimmed<<std::endl;
    std::cout<<"MC partner events within deltaT: "<<numtagged<<std::endl;

    std::cout<<"\n  EVindex 0_1: "<<ev01<<std::endl;
    std::cout<<" EVindex 0_2: "<<ev02<<std::endl;
    std::cout<<" EVindex 0_3: "<<ev03<<std::endl;

    std::cout<<"\n noise events between e+ an n: "<<badev1<<std::endl;
    std::cout<<"bad ev following positron mean time (ns): "<<totalbadevdeltaT/badev1<<std::endl;
  */

    // save file
    //file_input->Close();
    file_output->cd();

    // write objects...
    h_energy_full->Write();
    h_energy_prompt->Write();
    h_energy_late->Write();
    h_delT->Write();
    h_delR->Write();
    h2_delT_vs_delR->Write();
    h_Z_prompt->Write();
    h_Z_late->Write();
    h_RvsEnergy_prompt->Write();
    h_ZvsEnergy_prompt->Write();
    h_ZvsRho_prompt->Write();
    h_RvsEnergy_late->Write();
    h_ZvsEnergy_late->Write();
    h_ZvsRho_late->Write();
    h_XY_dis->Write();    
    
    /*h_after_cut_emc_nu.Write();
    h_after_cut_emc_p1.Write();
    h_after_cut_emc_p2.Write();
    h2_after_cut_emc_p2_vs_p1.Write();

    h_after_cut_deltaT.Write();
    h_after_cut_deltaR.Write();
    h2_after_cut_deltaT_vs_deltaR.Write();
    h_after_cut_deltaT_0_1.Write();
    h_after_cut_deltaT_0_2.Write();
    h_after_cut_deltaR_0_1.Write();
    h_after_cut_deltaR_0_2.Write();

    h_after_cut_efit_prompt.Write();
    h_after_cut_efit_delayed.Write();
    h2_after_cut_efit_delayed_vs_prompt.Write();

    delE_efit_prompt.Write();
    delE_emc_p1.Write();
    delE_emc_p1_0_1.Write();
    delE_emc_p1_0_2.Write();
    c_delekeplus_e1->Write();
    c_delekeplus_poske->Write();
    c_delekeplus_poske_1v2->Write();

    deltaTimeBadEVindex1.Write();
    deltaRBadEVindex1.Write();
    */
    //write csv output file with event numbers
    sprintf(name, "%s.csv",filename_output.c_str());

    //FILE *fOut = fopen(name,"w");
    //fprintf(fOut,"filename_input,filename_output,Rmax,Zmin,Zmax,energy_ep_min,energy_ep_max,energy_n_min,energy_n_max,deltaTmin,deltaTmax,deltaRmax,initial_entries,final_entries,finished\n");
    //fprintf(fOut,"%s,%s,%f,%f,%f,%f,%f,%f,%f,%i,%i,%f,%i,%i,%i\n", filename_input.c_str().str(), filename_output.c_str().str(), FV, z_cut1, z_cut2, energy1_lower, energy1_upper, energy2_lower, energy2_upper, delTcut_lower, delTcut, delRcut, nsimmed, ntagged,1);
    //fclose(fOut);
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
