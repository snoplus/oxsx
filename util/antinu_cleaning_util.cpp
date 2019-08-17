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

void process_cuts(const std::string filename_input, const std::string filename_output, Double_t energy_ep_min, Double_t energy_ep_max, Double_t energy_n_min, Double_t energy_n_max,Double_t deltaTmin, Double_t deltaTmax, Double_t promptRmax, Double_t lateRmax, Double_t deltaRmax){

    // load input file
    TFile *file_input = TFile::Open(filename_input.c_str());
    TNtuple *tree_input = (TNtuple*)file_input->Get("output");
    
    // setup output file
    TFile *file_output = new TFile(filename_output.c_str(), "RECREATE");
    TTree *tree_output = new TTree("nt","Tagged positron + neutron events");

    char *name = new char[1000];

    TH1D h_after_cut_emc_nu("h_after_cut_emc_nu", "Parent antinu KE (MeV)", 300, 0, 9);
    TH1D h_after_cut_emc_p1("h_after_cut_emc_p1", "Particle 1 KE (MeV)", 300, 0, 9);
    TH1D h_after_cut_emc_p2("h_after_cut_emc_p2", "Particle 2 KE (MeV)", 300, 0, 1);
    TH2D h2_after_cut_emc_p2_vs_p1("h2_after_cut_emc_p2_vs_p1", "Particle 2 KE vs Particle 1 KE", 10000, 0, 10, 1000, 0, 1);
    
    TH1D h_after_cut_deltaT("h_after_cut_deltaT", "Time diff (ns)", 300, 0, 1500000);
    TH1D h_after_cut_deltaR("h_after_cut_deltaR", "Inter-particle distance (mm)", 300, 0, 5000);
    TH2D h2_after_cut_deltaT_vs_deltaR("h2_after_cut_deltaR_deltaT", "deltaR vs deltaT ", 500, 0, 1000000, 300, 0, 10000);
    TH1D h_after_cut_deltaT_0_1("h_after_cut_deltaT_0_1", "Time diff (ns)   evindex=0,1", 100, 0, 1000000);
    TH1D h_after_cut_deltaT_0_2("h_after_cut_deltaT_0_2", "Time diff (ns)   evindex=0,2", 100, 0, 1000000);
    TH1D h_after_cut_deltaR_0_1("h_after_cut_deltaR_0_1", "Inter-particle distance (mm)  evindex=0,1", 300, 0, 10000);
    TH1D h_after_cut_deltaR_0_2("h_after_cut_deltaR_0_2", "Inter-particle distance (mm)  evindex=0,2", 300, 0, 10000);
    
    TH1D h_after_cut_efit_prompt("h_after_cut_efit_prompt", "Prompt Reconstructed Energy (MeV)", 300, 0, 9);
    TH1D h_after_cut_efit_delayed("h_after_cut_efit_delayed", "Delayed Reconstructed Energy (MeV)", 300, 0, 9);
    TH2D h2_after_cut_efit_delayed_vs_prompt("h2_after_cut_efit_delayed_vs_prompt", "Delayed vs Prompt Reconstructed Energy (MeV)", 300, 0, 9, 300, 0, 9);
    
    // Comparing Reconstructed to Truth:
    int nxbins = 16;
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
    

    Double_t neutron_capture_energy = 2.2;//1.857;
    Double_t e_rem = 0.784;
    Double_t annihilation = 1.022;
    
    //Iwan
    ULong64_t numsimmed = 1;
    ULong64_t numtagged = 0;
    ULong64_t ev01 = 0;
    ULong64_t ev02 = 0;
    ULong64_t ev03 = 0;
    ULong64_t badev1 = 0;
    Double_t totalbadevdeltaT = 0.;
    
    
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

    //Double_t mc_pos_r_nu, mc_pos_x_nu, mc_pos_y_nu, mc_pos_z_nu;
    //Double_t mc_pos_n_r, mc_pos_x_n, mc_pos_y_n, mc_pos_z_n;
    //Double_t mc_pos_r_ep, mc_pos_x_ep, mc_pos_y_ep, mc_pos_z_ep;
    Double_t mc_quench_i, mc_energy_nu, mc_energy_n, mc_energy_ep;
    UInt_t mc_time_days, mc_time_seconds;
    Double_t mc_time_nanoseconds;
    //Double_t latitude, longitude, altitude, distance;

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

    // values to modify
    tree_output->Branch("ev_fit_energy_p1", &ev_energy_p1);
    tree_output->Branch("ev_fit_energy_p2", &ev_energy_p2);
    tree_output->Branch("ev_index_p1", &ev_index_p1);
    tree_output->Branch("ev_index_p2", &ev_index_p2);


    ///////////////////////////////////////////////
    /// Cleaning+Coinc. Tagging (using evindex) ///
    ///////////////////////////////////////////////

    ULong64_t n_entries = tree_input->GetEntries();

    if (n_entries != 0){
      //get first event and its MCentry number
      tree_input->GetEntry(0);
      mc_entryb4 = mc_entry;  

      //go through entries and collect triggered evs which correspond to a single generated antinu MC event.
      for(ULong64_t i = 1; i < n_entries; i++){
	tree_input->GetEntry(i);
	if(abs(mc_entry - mc_entryb4)>= 1 && ev_index == 0){ //mark when passed group of triggered evs with same mcindex:
	  numsimmed += 1;
	  tree_input->GetEntry(i-1);//Get triggered event count for the antinu MC entry
	  ULong64_t EVcount = ev_index + 1; 

	  ULong64_t j = 0;  // j index for potential positron event
	  bool pairfound = false;
	  //move j index within MC index group, move onto next group when j reaches second to last event
	  while (j < EVcount - 1){
	    ULong64_t k = 1; // k index for potential neutron event
	    Double_t time_diff_condition = 0;
	    //move j index within MC index group, move on to next group if nothing within deltaT or no more events to look at
	    while((time_diff_condition < deltaTmax) && (k < EVcount - j)){
	      tree_input->GetEntry(i-EVcount+j+k); // store potential neutron event parameters
	      ev_next_validity = ev_validity;
	      ev_next_time_ns = ev_time_nanoseconds + (ev_time_seconds * pow(10, 9)) + (ev_time_days * 24 * 3600 * pow(10, 9));
	      ev_next_nhit = ev_nhit;
	      ev_next_index = ev_index;
	      ev_next_energy = ev_energy;
	      TVector3 ev_next_vecR = TVector3(ev_pos_x,ev_pos_y,ev_pos_z);
	  
	      tree_input->GetEntry(i-EVcount+j);// potential positron event
	      bool goodpair = false;
	      Double_t deltaT = 0.;
	      if (ev_validity && ev_next_validity){
		TVector3 ev_vecR= TVector3(ev_pos_x,ev_pos_y,ev_pos_z);
		ev_time_ns = ev_time_nanoseconds + (ev_time_seconds * pow(10, 9)) + (ev_time_days * 24 * 3600 * pow(10, 9));
		deltaT = std::fabs(ev_next_time_ns - ev_time_ns);
		Double_t deltaR = (ev_vecR - ev_next_vecR).Mag();
		//std::cout<<"\n MCentry: "<<mc_entry<<" EVindex: "<<ev_index<<" E: "<<ev_energy<<" R: "<<ev_vecR.Mag()<<std::endl;
		//std::cout<<"MCentry: "<<mc_entry<<" EVindex2: "<<ev_next_index<<" E: "<<ev_next_energy<<" R: "<<ev_next_vecR.Mag()<<std::endl;
		//std::cout<<"delatR: "<<deltaR<<" deltaT: "<<deltaT<<"\n"<<std::endl;
		if(deltaT > deltaTmin && deltaT < deltaTmax){
		  if (ev_vecR.Mag() < promptRmax){
		    if (ev_next_vecR.Mag() < lateRmax){
		      if (deltaR < deltaRmax){
			// if want to apply nhit cuts
			//if (nhit1Min <= nhit  && nhit <= nhit1Max){
			//if (nhit2Min <= nextnhit && nextnhit <= nhit2Max){
		
			//if want to correc energies
			//double corrected_ev_energy = ev_energy + correction(ev_energy);
			//if (energy_ep_min <= corrected_ev_energy && corrected_ev_energy <= energy_ep_max){
			//double corrected_ev_next_energy = nextEnergy + correction(nextEnergy);
			//if (energy_n_min <= corrected_ev_next_energy && corrected_ev_next_energy <= energy_n_max){
			if (energy_ep_min <= ev_energy && ev_energy <= energy_ep_max){
			  if (energy_n_min <= ev_next_energy && ev_next_energy <= energy_n_max){ //|| (4. <= nextEnergycorrec && nextEnergycorrec <= 5.8)){ //if want to allow for carbon capture?

			    numtagged += 1;		
			    //std::cout<<"tagged \n"<<std::endl;
			    //std::cout<<"\n MCentry: "<<MCentry<<" EVindex: "<<evindex<<" nhits: "<<nhit<<" day: "<<days<<std::endl;
			    //std::cout<<"MCentry2: "<<nextMCentry<<" EVindex2: "<<nextevindex<<" nhits2: "<<nextnhit<<" day2: "<<nextdays<<std::endl;
			    //std::cout<<"timeD: "<<timeD<<std::endl;
			    //std::cout<<"------------------------------------------------------- "<<timeD<<std::endl;
			    //std::cout<<" EVindex: "<<evindex<<std::endl;
			    //std::cout<<" EVindex2: "<<nextevindex<<std::endl;
			    //std::cout<<"Energy: "<<Energy<<" EDep: "<<EDep<<" QuenchEDep: "<<QuenchEDep<<std::endl;
			    h_after_cut_emc_nu.Fill(mc_energy_nu);
			    h_after_cut_emc_p1.Fill(mc_energy_ep);
			    h_after_cut_emc_p2.Fill(mc_energy_n);
			    h2_after_cut_emc_p2_vs_p1.Fill(mc_energy_ep, mc_energy_n);

			    h_after_cut_deltaT.Fill(deltaT);
			    h_after_cut_deltaR.Fill(deltaR);
			    h2_after_cut_deltaT_vs_deltaR.Fill(deltaR, deltaT);
			    h_after_cut_efit_prompt.Fill(ev_energy);
			    h_after_cut_efit_delayed.Fill(ev_next_energy);
			    h2_after_cut_efit_delayed_vs_prompt.Fill(ev_energy, ev_next_energy);

			    delE_efit_prompt.Fill(ev_energy,mc_energy_ep + annihilation - ev_energy);
			    delE_emc_p1.Fill(mc_energy_ep,mc_energy_ep + annihilation - ev_energy);

			    ev_energy_p1 = ev_energy;
			    ev_index_p1 = ev_index;
			    ev_energy_p2 = ev_next_energy;
			    ev_index_p2 = ev_next_index;
			    tree_output->Fill();
			  
			    if (ev_index == 0 && ev_next_index == 1){
			      ev01 += 1;
			      h_after_cut_deltaT_0_1.Fill(deltaT);
			      h_after_cut_deltaR_0_1.Fill(deltaR);

			      delE_emc_p1_0_1.Fill(mc_energy_ep,mc_energy_ep + annihilation - ev_energy);
			    }
			    
			    if (ev_index == 0 && ev_next_index == 2){
			      ev02 += 1;
			      h_after_cut_deltaT_0_2.Fill(deltaT);
			      h_after_cut_deltaR_0_2.Fill(deltaR);

			      delE_emc_p1_0_2.Fill(mc_energy_ep,mc_energy_ep + annihilation - ev_energy);

			      Double_t ev_time0_ns = ev_time_ns;
			      TVector3 ev_vecR0 = ev_vecR;
			      tree_input->GetEntry(i-EVcount+j+k-1);
			      if (ev_validity){
				badev1 += 1;
				Double_t ev_time1_ns = ev_time_nanoseconds + (ev_time_seconds * pow(10, 9)) + (ev_time_days * 24 * 3600 * pow(10, 9));
				TVector3 ev_vecR1 = TVector3(ev_pos_x,ev_pos_y,ev_pos_z);
				Double_t deltaR01 = (ev_vecR0 - ev_vecR1).Mag();
				Double_t deltaT01 = std::fabs(ev_time1_ns - ev_time0_ns);
			      
				deltaTimeBadEVindex1.Fill(deltaT01);
				deltaRBadEVindex1.Fill(deltaR01);
				totalbadevdeltaT += deltaT01;
			      }
			  
			    }
			    if (ev_index == 0 && ev_next_index == 3)
			      ev03 += 1;

			    goodpair = true;  //pair passed quality cuts, if not continue search in group
			    pairfound = true; // a pair was found
			    k += 100; // cancel sub search
			  }
			}
		      }
		    }
		  }
		}
	      } // if pair failed quality cuts, continue sub group search
	      if (!goodpair){
		k += 1;
		if (!ev_validity){ //if positron/j_index not a valid fit, move onto next sub group
		  k += 100;
		}else{ //else if good, check that it satisfies time diff, if not move onto next sub group
		  time_diff_condition = deltaT;
		}
	      }
	    }
	    if (pairfound){
	      j += 100;
	    }else{
	      j += 1;
	    }
	  }
	}
	mc_entryb4 = mc_entry;
      }
    
      // have to repeat process for very last mc entry
      ULong64_t lastentry = (n_entries - 1);
      tree_input->GetEntry(lastentry);

      ULong64_t EVcount = ev_index + 1; 

      ULong64_t j = 0;  // j index for potential positron event
      bool pairfound = false;
      //move j index within MC index group, move onto next group when j reaches second to last event
      while (j < EVcount - 1){
	ULong64_t k = 1; // k index for potential neutron event
	Double_t time_diff_condition = 0;
	//move j index within MC index group, move on to next group if nothing within deltaT or no more events to look at
	while((time_diff_condition < deltaTmax) && (k < EVcount - j)){
	  tree_input->GetEntry(lastentry-EVcount+j+k+1); // store potential neutron event parameters
	  ev_next_validity = ev_validity;
	  ev_next_time_ns = ev_time_nanoseconds + (ev_time_seconds * pow(10, 9)) + (ev_time_days * 24 * 3600 * pow(10, 9));
	  ev_next_nhit = ev_nhit;
	  ev_next_index = ev_index;
	  ev_next_energy = ev_energy;
	  TVector3 ev_next_vecR = TVector3(ev_pos_x,ev_pos_y,ev_pos_z);
	  
	  tree_input->GetEntry(lastentry+1-EVcount+j);// potential positron event
	  bool goodpair = false;
	  Double_t deltaT = 0.;
	  if (ev_validity && ev_next_validity){
	    TVector3 ev_vecR= TVector3(ev_pos_x,ev_pos_y,ev_pos_z);
	    ev_time_ns = ev_time_nanoseconds + (ev_time_seconds * pow(10, 9)) + (ev_time_days * 24 * 3600 * pow(10, 9));
	    deltaT = std::fabs(ev_next_time_ns - ev_time_ns);
	    Double_t deltaR = (ev_vecR - ev_next_vecR).Mag();
	    if(deltaT > deltaTmin && deltaT < deltaTmax){
	      if (ev_vecR.Mag() < promptRmax){
		if (ev_next_vecR.Mag() < lateRmax){
		  if (deltaR < deltaRmax){
		    // if want to apply nhit cuts
		    //if (nhit1Min <= nhit  && nhit <= nhit1Max){
		    //if (nhit2Min <= nextnhit && nextnhit <= nhit2Max){
		
		    //if want to correc energies
		    //double corrected_ev_energy = ev_energy + correction(ev_energy);
		    //if (energy_ep_min <= corrected_ev_energy && corrected_ev_energy <= energy_ep_max){
		    //double corrected_ev_next_energy = nextEnergy + correction(nextEnergy);
		    //if (energy_n_min <= corrected_ev_next_energy && corrected_ev_next_energy <= energy_n_max){
		    if (energy_ep_min <= ev_energy && ev_energy <= energy_ep_max){
		      if (energy_n_min <= ev_next_energy && ev_next_energy <= energy_n_max){ //|| (4. <= nextEnergycorrec && nextEnergycorrec <= 5.8)){ //if want to allow for carbon capture?
			
			numtagged += 1;		
		      
			h_after_cut_emc_nu.Fill(mc_energy_nu);
			h_after_cut_emc_p1.Fill(mc_energy_ep);
			h_after_cut_emc_p2.Fill(mc_energy_n);
			h2_after_cut_emc_p2_vs_p1.Fill(mc_energy_ep, mc_energy_n);

			h_after_cut_deltaT.Fill(deltaT);
			h_after_cut_deltaR.Fill(deltaR);
			h2_after_cut_deltaT_vs_deltaR.Fill(deltaR, deltaT);
			h_after_cut_efit_prompt.Fill(ev_energy);
			h_after_cut_efit_delayed.Fill(ev_next_energy);
			h2_after_cut_efit_delayed_vs_prompt.Fill(ev_energy, ev_next_energy);

			delE_efit_prompt.Fill(ev_energy,mc_energy_ep + annihilation - ev_energy);
			delE_emc_p1.Fill(mc_energy_ep,mc_energy_ep + annihilation - ev_energy);

			ev_energy_p1 = ev_energy;
			ev_index_p1 = ev_index;
			ev_energy_p2 = ev_next_energy;
			ev_index_p2 = ev_next_index;
			tree_output->Fill();
			  
			if (ev_index == 0 && ev_next_index == 1){
			  ev01 += 1;
			  h_after_cut_deltaT_0_1.Fill(deltaT);
			  h_after_cut_deltaR_0_1.Fill(deltaR);

			  delE_emc_p1_0_1.Fill(mc_energy_ep,mc_energy_ep + annihilation - ev_energy);
			  
			}
			
			if (ev_index == 0 && ev_next_index == 2){
			  ev02 += 1;
			  h_after_cut_deltaT_0_2.Fill(deltaT);
			  h_after_cut_deltaR_0_2.Fill(deltaR);

			  delE_emc_p1_0_2.Fill(mc_energy_ep,mc_energy_ep + annihilation - ev_energy);

			  Double_t ev_time0_ns = ev_time_ns;
			  TVector3 ev_vecR0 = ev_vecR;
			  tree_input->GetEntry(lastentry+1-EVcount+j+k-1);
			  if (ev_validity){
			    badev1 += 1;
			    Double_t ev_time1_ns = ev_time_nanoseconds + (ev_time_seconds * pow(10, 9)) + (ev_time_days * 24 * 3600 * pow(10, 9));
			    TVector3 ev_vecR1 = TVector3(ev_pos_x,ev_pos_y,ev_pos_z);
			    Double_t deltaR01 = (ev_vecR0 - ev_vecR1).Mag();
			    Double_t deltaT01 = std::fabs(ev_time1_ns - ev_time0_ns);
			      
			    deltaTimeBadEVindex1.Fill(deltaT01);
			    deltaRBadEVindex1.Fill(deltaR01);
			    totalbadevdeltaT += deltaT01;
			  }
			  
			}
			
			if (ev_index == 0 && ev_next_index == 3)
			  ev03 += 1;

			
			goodpair = true;  //pair passed quality cuts, if not continue search in group
			pairfound = true; // a pair was found
			k += 100; // cancel sub search
		      }
		    }
		  }
		}
	      }
	    }
	  } // if pair failed quality cuts, continue sub group search
	  if (!goodpair){
	    k += 1;
	    if (!ev_validity){ //if positron/j_index not a valid fit, move onto next sub group
	      k += 100;
	    }else{ //else if good, check that it satisfies time diff, if not move onto next sub group
	      time_diff_condition = deltaT;
	    }
	  }
	}
	if (pairfound){
	  j += 100;
	}else{
	  j += 1;
	}
      }
      // finished tagging


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

    
      // save file
      file_input->Close();
      file_output->cd();

      // write objects...
      h_after_cut_emc_nu.Write();
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
    
      //write csv output file with event numbers
      sprintf(name, "%s.csv",filename_output.c_str());
      FILE *fOut = fopen(name,"w");
      fprintf(fOut,"filename_input,filename_output,energy_ep_min,energy_ep_max,energy_n_min,energy_n_max,deltaTmin,deltaTmax,deltaRmax,initial_entries,final_entries,finished\n");
      fprintf(fOut,"%s,%s,%f,%f,%f,%f,%f,%f,%f,%i,%llu,%llu,%i\n", filename_input.c_str(), filename_output.c_str(), energy_ep_min, energy_ep_max, energy_n_min, energy_n_max, deltaTmin, deltaTmax, promptRmax, lateRmax, deltaRmax, numsimmed, numtagged,1);
      fclose(fOut);

      std::cout<<"fin: "<<filename_input<<" \n fout: "<<filename_output<<"\n e1min "<<energy_ep_min<<"\n e1max "<<energy_ep_max<<"\n e2min "<<energy_n_min<<"\n e2max "<<energy_n_max<<"\n delT "<<deltaTmin<<"\n delTmax "<<deltaTmax<<"\n R1 "<<promptRmax<<"\n R2 "<<lateRmax<<"\n delR "<<deltaRmax<<std::endl;
    }
    tree_output->AutoSave();
    file_output->Close();
    
}
