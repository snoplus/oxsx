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

/*double correction(double x){
  double correc = 8.83411e-02 + (9.02287e-03)*x;
  return correc;
}
*/

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

void process_cuts(const std::string filename_input, const std::string filename_output, Double_t energy_min, Double_t energy_max, Double_t Rmax){

  // load input file
  char *name = new char[1000];
  sprintf(name, "%s/output",filename_input.c_str());
  //TFile *file_input = TFile::Open(filename_input.c_str()); // for loading ntuples directly (old way)
  //TNtuple *tree_input = (TNtuple*)file_input->Get("output"); // for loading ntuples directly (old way)
  TChain *tree_input = new TChain("output");
  tree_input->Add(name);
  // setup output file
  TFile *file_output = new TFile(filename_output.c_str(), "RECREATE");
  TTree *tree_output = new TTree("nt","Tagged positron + neutron events");

  /*TH1D h_after_cut_emc_nu("h_after_cut_emc_nu", "Parent antinu KE (MeV)", 300, 0, 9);
  TH1D h_after_cut_emc_p1("h_after_cut_emc_p1", "Particle 1 KE (MeV)", 300, 0, 9);
  TH1D h_after_cut_emc_p2("h_after_cut_emc_p2", "Particle 2 KE (MeV)", 300, 0, 1);
  TH2D h2_after_cut_emc_p2_vs_p1("h2_after_cut_emc_p2_vs_p1", "Particle 2 KE vs Particle 1 KE", 10000, 0, 10, 1000, 0, 1);
  */
  TH1D h_after_cut_R("h_after_cut_R", "Particle Radius (mm)", 300, 0, 2);
  //TH1D h_after_cut_deltaR_0_1("h_after_cut_deltaR_0_1", "Inter-particle distance (mm)  evindex=0,1", 300, 0, 10000);
  //TH1D h_after_cut_deltaR_0_2("h_after_cut_deltaR_0_2", "Inter-particle distance (mm)  evindex=0,2", 300, 0, 10000);

  TH1D h_after_cut_efit("h_after_cut_efit", "Reconstructed Energy (MeV)", 300, 0, 9);

  double classifier_del_lhmax = 200;
  double classifier_del_lhmin = -200.;
  size_t num_bins_classifier = 400;
  TH1D h_berkeleyAlphaBeta("h_berkeleyAlphaBeta", "BerkeleyAlphaBeta", num_bins_classifier, classifier_del_lhmin, classifier_del_lhmax);
  TH1D h_alphaBeta212("h_alphaBeta212", "AlphaBeta212", num_bins_classifier, classifier_del_lhmin, classifier_del_lhmax);
  TH1D h_alphaBeta214("h_alphaBeta214", "AlphaBeta214", num_bins_classifier, classifier_del_lhmin, classifier_del_lhmax);
  
  // Comparing Reconstructed to Truth:
  /*int nxbins = 16;
  double e1min = 0.;
  double e1max = 8.;
  int nybins = 500;
  double delemin = 0.;
  double delemax = 1.;

  TH2D delE_efit_prompt("delE_efit_prompt","delE_efit_prompt",nxbins,e1min,e1max,nybins,delemin,delemax);
  TH2D delE_emc_p1("delE_emc_p1","delE_emc_p1",nxbins,e1min,e1max,nybins,delemin,delemax);
  TH2D delE_emc_p1_0_1("delE_emc_p1_0_1","delE_emc_p1_0_1",nxbins,e1min,e1max,nybins,delemin,delemax);
  TH2D delE_emc_p1_0_2("delE_emc_p1_0_2","delE_emc_p1_0_2",nxbins,e1min,e1max,nybins,delemin,delemax);
  
  // Investing in-between event (evindex = 1 when positron is 0th evindec and neutron is 2nd
  TH1D deltaTimeBadEVindex1("deltaTimeBadEVindex1","deltaTimeBadEVindex1 (ns)",20,0,1000);
  TH1D deltaRBadEVindex1("deltaRBadEVindex1","deltaRBadEVindex1 (mm)",300,0,10000);
  */

  Double_t neutron_capture_energy = 2.2;//1.857;
  Double_t e_rem = 0.784;
  Double_t annihilation = 1.022;

  //Iwan
  //ULong64_t numsimmed = 1;
  ULong64_t numtagged = 0;
  //ULong64_t ev01 = 0;
  //ULong64_t ev02 = 0;
  //ULong64_t ev03 = 0;
  //ULong64_t badev1 = 0;
  

  // properties to load from input ntuple into output ttree

  // mc parent and neutron truth positions not in standard ntuples?
  // don't need long, lat, altitude, distance here - distances for each reactor saved in a separate txt file

  //Double_t mc_pos_r, mc_pos_x, mc_pos_y, mc_pos_z;
  //input and/or output
  Int_t mc_entry, mc_entryb4;
  Double_t ev_pos_x, ev_pos_y, ev_pos_z;
  Double_t ev_energy, ev_next_energy, ev_energy_p1, ev_energy_p2;
  Int_t ev_time_seconds, ev_time_days, ev_time_nanoseconds;
  Double_t ev_time_ns, ev_next_time_ns;
  Int_t ev_nhit, ev_next_nhit;
  Bool_t ev_validity,ev_next_validity;
  Int_t ev_index, ev_next_index, ev_index_p1, ev_index_p2;
  TString *reactor_core_name = 0;

  //Double_t mc_pos_r_nu, mc_pos_x_nu, mc_pos_y_nu, mc_pos_z_nu;
  //Double_t mc_pos_n_r, mc_pos_x_n, mc_pos_y_n, mc_pos_z_n;
  //Double_t mc_pos_r_ep, mc_pos_x_ep, mc_pos_y_ep, mc_pos_z_ep;
  Double_t mc_quench_i, mc_energy_nu, mc_energy_n, mc_energy_ep;
  UInt_t mc_time_days, mc_time_seconds;
  ULong64_t ev_clock50, ev_next_clock50;
  Double_t mc_time_nanoseconds;
  Double_t ev_berkeleyAlphaBeta, ev_alphaBeta212, ev_alphaBeta214,\
           ev_next_berkeleyAlphaBeta, ev_next_alphaBeta212, ev_next_alphaBeta214;
  //Double_t latitude, longitude, altitude, distance;

  UInt_t pdg1, pdg2;
  UInt_t parentpdg1, parentpdg2;
  
  // set branches input ntuple
  tree_input->SetBranchAddress("mcIndex", &mc_entry);
  tree_input->SetBranchAddress("mcEdepQuenched", &mc_quench_i);
  tree_input->SetBranchAddress("parentKE1", &mc_energy_nu);
  tree_input->SetBranchAddress("mcke1", &mc_energy_ep);
  tree_input->SetBranchAddress("mcke2", &mc_energy_n);
  //tree_input->SetBranchAddress("mcPosr", &mc_pos_r);
  //tree_input->SetBranchAddress("mcPosx", &mc_pos_x);
  //tree_input->SetBranchAddress("mcPosy", &mc_pos_y);
  //tree_input->SetBranchAddress("mcPosz", &mc_pos_z);
  tree_input->SetBranchAddress("evIndex", &ev_index);
  tree_input->SetBranchAddress("energy", &ev_energy);
  tree_input->SetBranchAddress("fitValid", &ev_validity);
  tree_input->SetBranchAddress("posx", &ev_pos_x);
  tree_input->SetBranchAddress("posy", &ev_pos_y);
  tree_input->SetBranchAddress("posz", &ev_pos_z);
  tree_input->SetBranchAddress("nhits", &ev_nhit);
  tree_input->SetBranchAddress("uTDays", &ev_time_days);
  tree_input->SetBranchAddress("uTSecs", &ev_time_seconds);
  tree_input->SetBranchAddress("uTNSecs", &ev_time_nanoseconds);
  tree_input->SetBranchAddress("parentMeta1", &reactor_core_name);
  tree_input->SetBranchAddress("clockCount50", &ev_clock50);
  tree_input->SetBranchAddress("berkeleyAlphaBeta", &ev_berkeleyAlphaBeta);
  tree_input->SetBranchAddress("alphaBeta212", &ev_alphaBeta212);
  tree_input->SetBranchAddress("alphaBeta214", &ev_alphaBeta214);

  tree_input->SetBranchAddress("pdg1", &pdg1);
  tree_input->SetBranchAddress("pdg2", &pdg2);
  tree_input->SetBranchAddress("parentpdg1", &parentpdg1);
  tree_input->SetBranchAddress("parentpdg2", &parentpdg2);
  
  // set branches output pruned ntuple
  tree_output->Branch("entry", &mc_entry);
  //tree_output->Branch("mc_time_days", &mc_time_days);
  //tree_output->Branch("mc_time_seconds", &mc_time_seconds);
  //tree_output->Branch("mc_time_nanoseconds", &mc_time_nanoseconds);
  tree_output->Branch("mc_quench", &mc_quench_i);
  tree_output->Branch("mc_neutrino_energy", &mc_energy_nu);
  tree_output->Branch("mc_positron_energy", &mc_energy_ep);
  tree_output->Branch("mc_neutron_energy", &mc_energy_n);
  //tree_output->Branch("mc_neutrino_position_r", &mc_pos_r_nu);
  //tree_output->Branch("mc_neutrino_position_x", &mc_pos_x_nu);
  //tree_output->Branch("mc_neutrino_position_y", &mc_pos_y_nu);
  //tree_output->Branch("mc_neutrino_position_z", &mc_pos_z_nu);
  //tree_output->Branch("mc_positron_position_r", &mc_pos_r_ep);
  //tree_output->Branch("mc_positron_position_x", &mc_pos_x_ep);
  //tree_output->Branch("mc_positron_position_y", &mc_pos_y_ep);
  //tree_output->Branch("mc_positron_position_z", &mc_pos_z_ep);
  //tree_output->Branch("mc_neutron_position_r", &mc_pos_n_r);
  //tree_output->Branch("mc_neutron_position_x", &mc_pos_x_n);
  //tree_output->Branch("mc_neutron_position_y", &mc_pos_y_n);
  //tree_output->Branch("mc_neutron_position_z", &mc_pos_z_n);
  tree_output->Branch("ev_index", &ev_index);
  tree_output->Branch("ev_fit_energy", &ev_energy);
  tree_output->Branch("ev_fit_validity", &ev_validity);
  tree_output->Branch("ev_fit_position_x", &ev_pos_x);
  tree_output->Branch("ev_fit_position_y", &ev_pos_y);
  tree_output->Branch("ev_fit_position_z", &ev_pos_z);
  tree_output->Branch("ev_nhit", &ev_nhit);
  tree_output->Branch("ev_time_days", &ev_time_days);
  tree_output->Branch("ev_time_seconds", &ev_time_seconds);
  tree_output->Branch("ev_time_nanoseconds", &ev_time_nanoseconds);
  //tree_output->Branch("reactor_info_latitude", &latitude);
  //tree_output->Branch("reactor_info_longitude", &longitude);
  //tree_output->Branch("reactor_info_altitude", &altitude);
  //tree_output->Branch("reactor_info_distance", &distance);
  tree_output->Branch("reactor_core_name", &reactor_core_name);

  // values to modify
  tree_output->Branch("ev_fit_energy_p1", &ev_energy_p1);
  //tree_output->Branch("ev_fit_energy_p2", &ev_energy_p2);
  tree_output->Branch("ev_index_p1", &ev_index_p1);
  //tree_output->Branch("ev_index_p2", &ev_index_p2);


  ///////////////////////////////////////////////
  /// Cleaning+Coinc. Tagging (using evindex) ///
  ///////////////////////////////////////////////

  ULong64_t n_entries = tree_input->GetEntries();

  if (n_entries != 0){
    //go through entries and collect triggered evs which correspond to a single generated antinu MC event.
    n_entries = 500;
    for(ULong64_t i = 0; i < n_entries; i++){
      if (i % 100 == 0) std::cout <<  "    Processed entries: " << i << " (" << (double)i/n_entries*100. << "%) " << std::endl;
      tree_input->GetEntry(i);
      // reset pointers
      //std::cout<<"parentpdg1: "<<parentpdg1<<" parentpdg2: "<<parentpdg2<<"\n"<<std::endl;
      //std::cout<<"pdg1: "<<pdg1<<" pdg2: "<<pdg2<<"\n"<<std::endl;
      if (ev_index >= 0){
	if (ev_validity){
	  TVector3 ev_vecR= TVector3(ev_pos_x,ev_pos_y,ev_pos_z);
	  //std::cout<<"\n MCentry: "<<mc_entry<<" EVindex: "<<ev_index<<" E: "<<ev_energy<<" R: "<<ev_vecR.Mag()<<std::endl;
	  //std::cout<<"MCentry: "<<mc_entry<<" EVindex2: "<<ev_next_index<<" E: "<<ev_next_energy<<" R: "<<ev_next_vecR.Mag()<<std::endl;
	  //std::cout<<"delatR: "<<deltaR<<" deltaT: "<<deltaT<<"\n"<<std::endl;
	  if (ev_vecR.Mag() < Rmax){
	    //if want to correc energies
	    //double corrected_ev_energy = ev_energy + correction(ev_energy);
	    //if (energy_ep_min <= corrected_ev_energy && corrected_ev_energy <= energy_ep_max){
	    //double corrected_ev_next_energy = nextEnergy + correction(nextEnergy);
	    //if (energy_n_min <= corrected_ev_next_energy && corrected_ev_next_energy <= energy_n_max){
	    if (energy_min <= ev_energy && ev_energy <= energy_max){
	      if (ev_index == 0)
		std::cout<<"\nev index: "<<ev_index<<"  energy: "<<ev_energy<<std::endl;
	      else
		std::cout<<"ev index: "<<ev_index<<"  energy: "<<ev_energy<<std::endl;
	      numtagged += 1;
	      //std::cout<<"tagged \n"<<std::endl;
	      //std::cout<<"\n MCentry: "<<MCentry<<" EVindex: "<<evindex<<" nhits: "<<nhit<<" day: "<<days<<std::endl;
	      //std::cout<<"MCentry2: "<<nextMCentry<<" EVindex2: "<<nextevindex<<" nhits2: "<<nextnhit<<" day2: "<<nextdays<<std::endl;
	      //std::cout<<"timeD: "<<timeD<<std::endl;
	      //std::cout<<"------------------------------------------------------- "<<timeD<<std::endl;
	      //std::cout<<" EVindex: "<<evindex<<std::endl;
	      //std::cout<<" EVindex2: "<<nextevindex<<std::endl;
	      //std::cout<<" BerkAlphBet: "<<ev_berkeleyAlphaBeta<<std::endl;

	      /*h_after_cut_emc_nu.Fill(mc_energy_nu);
	      h_after_cut_emc_p1.Fill(mc_energy_ep);
	      h_after_cut_emc_p2.Fill(mc_energy_n);
	      h2_after_cut_emc_p2_vs_p1.Fill(mc_energy_ep, mc_energy_n);
	      */          

	      h_after_cut_efit.Fill(ev_energy);
	      h_after_cut_R.Fill(pow(ev_vecR.Mag()/6000.,3));
	      
	      h_berkeleyAlphaBeta.Fill(ev_berkeleyAlphaBeta);
	      h_alphaBeta212.Fill(ev_alphaBeta212);
	      h_alphaBeta214.Fill(ev_alphaBeta214);
			    
	      //delE_efit_prompt.Fill(ev_energy,mc_energy_ep + annihilation - ev_energy);
	      //delE_emc_p1.Fill(mc_energy_ep,mc_energy_ep + annihilation - ev_energy);

	      ev_energy_p1 = ev_energy;
	      ev_index_p1 = ev_index;
	      //ev_energy_p2 = ev_next_energy;
	      //ev_index_p2 = ev_next_index;
	      tree_output->Fill();
	    }
	  }
	}
      }
    }

    // finished tagging
    
    //summary stats
    std::cout<<"\n \n Summary: "<<std::endl;
    std::cout<<"MC events tagged: "<<numtagged<<std::endl;

    // save file
    //file_input->Close();
    file_output->cd();

    // write objects...
    /*h_after_cut_emc_nu.Write();
    h_after_cut_emc_p1.Write();
    h_after_cut_emc_p2.Write();
    h2_after_cut_emc_p2_vs_p1.Write();
    */
    h_after_cut_efit.Scale(1./(double)h_after_cut_efit.Integral());
    h_after_cut_R.Scale(1./(double)h_after_cut_R.Integral());    
    
    h_berkeleyAlphaBeta.Scale(1./(double)h_berkeleyAlphaBeta.Integral());
    h_alphaBeta212.Scale(1./(double)h_alphaBeta212.Integral());
    h_alphaBeta214.Scale(1./(double)h_alphaBeta214.Integral());

    h_after_cut_R.Write();

    h_after_cut_efit.Write();

    h_berkeleyAlphaBeta.Write();
    h_alphaBeta212.Write();
    h_alphaBeta214.Write();

    //write csv output file with event numbers
    /*sprintf(name, "%s.csv",filename_output.c_str());
    FILE *fOut = fopen(name,"w");
    fprintf(fOut,"filename_input,filename_output,energy_ep_min,energy_ep_max,energy_n_min,energy_n_max,deltaTmin,deltaTmax,deltaRmax,initial_entries,final_entries,finished\n");
    fprintf(fOut,"%s,%s,%f,%f,%f,%f,%f,%f,%f,%i,%llu,%llu,%i\n", filename_input.c_str(), filename_output.c_str(), energy_ep_min, energy_ep_max, energy_n_min, energy_n_max, deltaTmin, deltaTmax, promptRmax, lateRmax, deltaRmax, numsimmed, numtagged,1);
    fclose(fOut);
    */
    std::cout<<"fin: "<<filename_input<<" \n fout: "<<filename_output<<"\n emin "<<energy_min<<"\n emax "<<energy_max<<"\n Rmax "<<Rmax<<std::endl;
  }
  tree_output->AutoSave();
  file_output->Close();
  delete tree_input;
  delete file_output;
}
