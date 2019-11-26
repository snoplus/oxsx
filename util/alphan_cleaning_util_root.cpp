#include <TFile.h>
#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/FitResult.hh>
#include <RAT/DB.hh>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TMath.h>
#include <TCanvas.h>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <TNtuple.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <RAT/DU/Utility.hh>
#include <RAT/TrackNav.hh>
#include <RAT/TrackCursor.hh>
#include <RAT/TrackNode.hh>
#include "TView.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include <RAT/DB.hh>
#include <TObject.h>
#include <algorithm>
#include <TRandom3.h>

void process_cuts(const std::string filename_input, const std::string filename_output, Double_t energy_ep_min, Double_t energy_ep_max, Double_t energy_n_min, Double_t energy_n_max,Double_t deltaTmin, Double_t deltaTmax, Double_t promptRmax, Double_t lateRmax, Double_t deltaRmax) {
  
  TH1D h_edep("h_edep", "Deposited Energy (MeV)", 500, 0, 15);
  TH1D h_equench("h_equench", "Quenched Energy (MeV)", 500, 0, 15);
  TH1D h_part1ke("h_part1ke", "Particle 1 KE (neutron) (MeV)", 500, 0, 15);
  TH1D h_part2ke("h_part2ke", "Particle 2 (alpha) KE (MeV)", 1000, 0, 10);
  TH1D h_edep_01("h_edep", "Deposited Energy (MeV)", 500, 0, 15);
  TH1D h_equench_01("h_equench", "Quenched Energy (MeV)", 500, 0, 15);
  TH1D h_part1ke_01("h_part1ke_01", "Particle 1 KE (neutron) (MeV)", 500, 0, 15);
  TH1D h_part2ke_01("h_part2ke_01", "Particle 2 (alpha) KE (MeV)", 500, 0, 15);
  TH1D h_edep_02("h_edep", "Deposited Energy (MeV)", 500, 0, 15);
  TH1D h_equench_02("h_equench", "Quenched Energy (MeV)", 500, 0, 15);
  TH1D h_part1ke_02("h_part1ke_02", "Particle 1 KE (neutron) (MeV)", 500, 0, 15);
  TH1D h_part2ke_02("h_part2ke_02", "Particle 2 (alpha) KE (MeV)", 500, 0, 15);
  TH2D h_ke1vske2("h_ke1vske2", "Particle 1 KE (neutron) vs alpha ke(MeV)", 500, 0, 15, 500, 0, 15);
  TH2D h_ke1vske2_01("h_ke1vske2_01", "Particle 1 KE (neutron) vs alpha ke(MeV)", 500, 0, 15, 500, 0, 15);
  TH2D h_ke1vske2_02("h_ke1vske2_02", "Particle 1 KE (neutron) vs alpha ke(MeV)", 500, 0, 15, 500, 0, 15);
  
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

  size_t time_res_bins = 850;
  double h_time_res_min = -250;
  double h_time_res_max = 600;

  TH1D Prompt_index0_1TimeRes("Prompt_index0_1TimeRes","Prompt_index0_1TimeRes",time_res_bins,h_time_res_min,h_time_res_max);
  TH1D Late_index0_1TimeRes("Late_index0_1TimeRes","Late_index0_1TimeRes",time_res_bins,h_time_res_min,h_time_res_max);
  TH1D Prompt_index0_2TimeRes("Prompt_index0_2TimeRes","Prompt_index0_2TimeRes",time_res_bins,h_time_res_min,h_time_res_max);
  TH1D Late_index0_2TimeRes("Late_index0_2TimeRes","Late_index0_2TimeRes",time_res_bins,h_time_res_min,h_time_res_max);
  
  TH1D BadEvPos1in1_index0_2TimeRes("BadEvPos1in1_index0_2TimeRes","BadEvPos1in1_index0_2TimeRes",time_res_bins,h_time_res_min,h_time_res_max);
  TH1D BadEvPos2in1_index0_2TimeRes("BadEvPos2in1_index0_2TimeRes","BadEvPos2in1_index0_2TimeRes",time_res_bins,h_time_res_min,h_time_res_max);

  double classifier_del_lhmax = 200;
  double classifier_del_lhmin = -200.;
  size_t num_bins_classifier = 600;

  TH1D h_prompt_berkeleyAlphaBeta("h_prompt_erkeleyAlphaBeta", "BerkeleyAlphaBeta", num_bins_classifier, classifier_del_lhmin, classifier_del_lhmax);
  TH1D h_prompt_berkeleyAlphaBeta_01("h_prompt_berkeleyAlphaBeta_01", "BerkeleyAlphaBeta", num_bins_classifier, classifier_del_lhmin, classifier_del_lhmax);
  TH1D h_prompt_berkeleyAlphaBeta_02("h_prompt_berkeleyAlphaBeta_02", "BerkeleyAlphaBeta", num_bins_classifier, classifier_del_lhmin, classifier_del_lhmax);

  TH1D h_prompt_ITR("h_prompt_ITR", "ITR", 300, -1.5, 1.5  );
  TH1D h_prompt_ITR_01("h_prompt_ITR_01", "ITR", 300, -1.5, 1.5  );
  TH1D h_prompt_ITR_02("h_prompt_ITR_02", "ITR", 300, -1.5, 1.5  );

  TH1D h_prompt_Beta14("h_prompt_Beta14", "Beta14", 400, -1, 3);
  TH1D h_prompt_Beta14_01("h_prompt_Beta14_01", "Beta14", 400, -1, 3);
  TH1D h_prompt_Beta14_02("h_prompt_Beta14_02", "Beta14", 400, -1, 3);

  TH1D h_late_berkeleyAlphaBeta("h_late_erkeleyAlphaBeta", "BerkeleyAlphaBeta", num_bins_classifier, classifier_del_lhmin, classifier_del_lhmax);
  TH1D h_late_berkeleyAlphaBeta_01("h_late_berkeleyAlphaBeta_01", "BerkeleyAlphaBeta", num_bins_classifier, classifier_del_lhmin, classifier_del_lhmax);
  TH1D h_late_berkeleyAlphaBeta_02("h_late_berkeleyAlphaBeta_02", "BerkeleyAlphaBeta", num_bins_classifier, classifier_del_lhmin, classifier_del_lhmax);

  TH1D h_late_ITR("h_late_ITR", "ITR", 300, -1.5, 1.5  );
  TH1D h_late_ITR_01("h_late_ITR_01", "ITR", 300, -1.5, 1.5  );
  TH1D h_late_ITR_02("h_late_ITR_02", "ITR", 300, -1.5, 1.5  );

  TH1D h_late_Beta14("h_late_Beta14", "Beta14", 400, -1, 3);
  TH1D h_late_Beta14_01("h_late_Beta14_01", "Beta14", 400, -1, 3);
  TH1D h_late_Beta14_02("h_late_Beta14_02", "Beta14", 400, -1, 3);

  
  int Fitted, nextFitted;
  Double_t timeD;
  
  Int_t mc_entry, mc_entryb4;
  Double_t ev_pos_x, ev_pos_y, ev_pos_z;
  Double_t ev_energy, ev_next_energy, ev_energy_p1, ev_energy_p2;
  Int_t ev_time_seconds, ev_time_days, ev_time_nanoseconds;
  Double_t ev_time_ns, ev_next_time_ns;
  Int_t ev_nhit, ev_next_nhit;
  Bool_t ev_validity,ev_next_validity;
  Int_t ev_index, ev_next_index, ev_index_p1, ev_index_p2;
  TString *reactor_core_name = 0;

  Double_t mc_quench_i, mc_energy_nu, mc_energy_n, mc_energy_ep;
  UInt_t mc_time_days, mc_time_seconds;
  ULong64_t ev_clock50;

  TVector3 ev_vecR, ev_next_vecR;

  Double_t ev_berkeleyAB_lh, ev_next_berkeleyAB_lh;
  Double_t ev_beta14 = 0.;
  Double_t ev_next_beta14 = 0.;
  Double_t ev_itr, ev_next_itr;

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

  Double_t mc_time_nanoseconds;

  RAT::DU::DSReader ds(filename_input.c_str());
  
  size_t n_entries = ds.GetEntryCount();
  std::cout<<"num of events: "<<n_entries<<std::endl;

  RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path
  const RAT::DU::GroupVelocity& groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity(); // To get the group velocity
  const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo(); // The PMT positions etc..

  for(size_t ie = 0; abs(ie) < abs(n_entries); ie++){  
    
    if (ie % 100 == 0) std::cout <<  "    Processed events: " << ie << " (" << (double)ie/n_entries*100. << "%) " << std::endl;
    
    const RAT::DS::Entry& rds = ds.GetEntry(ie);
    const RAT::DS::MC& rMC = rds.GetMC();
    //const RAT::DS::MCParticle& rmcparent = rMC.GetMCParent(0);

    bool is_neutron = false;
    bool is_alpha = false;

    double equench = rMC.GetScintQuenchedEnergyDeposit();
    double edep = rMC.GetScintEnergyDeposit();
    double part1ke;
    double part2ke;

    size_t particle_count = rMC.GetMCParticleCount();
    const RAT::DS::MCParticle& rmcparticle1 = rMC.GetMCParticle(0);
    size_t particle1pdg = rmcparticle1.GetPDGCode();
    
    if (particle_count > 1){
      const RAT::DS::MCParticle& rmcparticle2 = rMC.GetMCParticle(1);
      size_t particle2pdg = rmcparticle2.GetPDGCode();
      if (particle2pdg == 1000020040){
	part2ke = rmcparticle2.GetKineticEnergy();
	is_alpha = true;
      }
    }
    //std::cout<<"part1 pdg: "<<particle1pdg<<" particle2 pdg: "<<particle2pdg<<std::endl;

    if (particle1pdg == 2112){
      part1ke = rmcparticle1.GetKineticEnergy();
      is_neutron = true;
    }

    //double ParKE = rmcparent.GetKineticEnergy();
    int EVcount = rds.GetEVCount();
    //evcounttot += EVcount;
    size_t j = 0;  // j index for potential positron event
    bool pairfound = false;
    //move j index within MC index group, move onto next group when j reaches second to last event
    while (abs(j) < abs(EVcount - 1)){
      int k = 1; // k index for potential neutron event
      float time_diff_condition = 0;
      //move j index within MC index group, move on to next group if nothing within deltaT or no more events to look at
      while((time_diff_condition < deltaTmax) && (k < abs(EVcount - j))){
	const RAT::DS::EV& rEV2 = rds.GetEV(j+k);// neutron event
	bool BadFit2 = true;
	if( rEV2.FitResultExists("scintFitter")) {
	  if( rEV2.GetFitResult("scintFitter").GetVertex(0).ContainsPosition() ) {
	    if( rEV2.GetFitResult("scintFitter").GetVertex(0).ValidPosition() ) {
	      if( rEV2.GetFitResult("scintFitter").GetVertex(0).ContainsEnergy() ) {
		if( rEV2.GetFitResult("scintFitter").GetVertex(0).ValidEnergy() ) {
		  
		  BadFit2 = false;
		  nextFitted = 1;
		  RAT::DS::FitResult fResult = rEV2.GetFitResult("scintFitter");
		  RAT::DS::FitVertex fVertex = fResult.GetVertex(0);
		  ev_time_days = rEV2.GetUniversalTime().GetDays();
		  ev_time_seconds = rEV2.GetUniversalTime().GetSeconds();
		  ev_time_nanoseconds = rEV2.GetUniversalTime().GetNanoSeconds();
		  ev_next_time_ns = ev_time_nanoseconds + (ev_time_seconds * pow(10, 9)) + (ev_time_days * 24 * 3600 * pow(10, 9));
		  ev_next_energy = fVertex.GetEnergy();
		  ev_next_vecR = fVertex.GetPosition();

		  const RAT::DS::ClassifierResult &result = rEV2.GetClassifierResult("BerkeleyAlphaBeta:scintFitter");
		  ev_next_berkeleyAB_lh = result.GetClassification("likelihood");
		  const RAT::DS::ClassifierResult &result_itr = rEV2.GetClassifierResult("ITR:scintFitter");
		  ev_next_itr = result_itr.GetClassification("ITR");
		  //const RAT::DS::ClassifierResult &result_beta14 = rEV2.GetClassifierResult("isotropy:scintFitter");
		  //ev_next_beta14 = result_beta14.GetClassification("snobeta14");
		  
		}
	      }
	    }
	  }
	}
	if (BadFit2)
	  nextFitted = -1;
	
	const RAT::DS::EV& rEV1 = rds.GetEV(j);// positron event
	bool BadFit1 = true;
	  
	if( rEV1.FitResultExists("scintFitter")) {
	  if( rEV1.GetFitResult("scintFitter").GetVertex(0).ContainsPosition() ) {
	    if( rEV1.GetFitResult("scintFitter").GetVertex(0).ValidPosition() ) {
	      if( rEV1.GetFitResult("scintFitter").GetVertex(0).ContainsEnergy() ) {
		if( rEV1.GetFitResult("scintFitter").GetVertex(0).ValidEnergy() ) {
		  BadFit1 = false;
		  Fitted = 1;
		  RAT::DS::FitResult fResult = rEV1.GetFitResult("scintFitter");
		  RAT::DS::FitVertex fVertex = fResult.GetVertex(0);
		  ev_time_days = rEV1.GetUniversalTime().GetDays();
		  ev_time_seconds = rEV1.GetUniversalTime().GetSeconds();
		  ev_time_nanoseconds = rEV1.GetUniversalTime().GetNanoSeconds();
		  ev_time_ns = ev_time_nanoseconds + (ev_time_seconds * pow(10, 9)) + (ev_time_days * 24 * 3600 * pow(10, 9));
		  ev_energy = fVertex.GetEnergy();
		  ev_vecR = fVertex.GetPosition();

		  const RAT::DS::ClassifierResult &result = rEV1.GetClassifierResult("BerkeleyAlphaBeta:scintFitter");
		  ev_berkeleyAB_lh = result.GetClassification("likelihood");
		  const RAT::DS::ClassifierResult &result_itr = rEV1.GetClassifierResult("ITR:scintFitter");
		  ev_itr = result_itr.GetClassification("ITR");
		  //const RAT::DS::ClassifierResult &result_beta14 = rEV1.GetClassifierResult("Beta14:scintFitter");
		  //ev_beta14 = result_beta14.GetClassification("isotropy");
		}
	      }
	    }
	  }
	}
	if (BadFit1)
	  Fitted = -1;
	  
	bool goodpair = false;
	timeD = 0;
	//std::cout<<"MCentry2: "<<ie<<" EVindex2: "<<j+k<<" nhits2: "<<nextnhit<<" day2: "<<nextdays<<std::endl;
	//std::cout<<"\n MCentry: "<<ie<<" EVindex: "<<j<<" nhits: "<<nhit<<" day: "<<days<<std::end;

	if (Fitted > 0 && nextFitted > 0){
	  timeD = std::fabs(ev_next_time_ns - ev_time_ns);
	  if(timeD > 400 && timeD < deltaTmax){
	    if (ev_vecR.Mag() < promptRmax){
	      if (ev_next_vecR.Mag() < lateRmax){
		Double_t deltaR = (ev_vecR - ev_next_vecR).Mag();
		if (deltaR < deltaRmax){
		  if (energy_ep_min <= ev_energy && ev_energy <= energy_ep_max){
		    if (energy_n_min <= ev_next_energy && ev_next_energy <= energy_n_max){
		      numtagged += 1;
		      const RAT::DS::CalPMTs& calibratedPMTs2 = rEV2.GetCalPMTs();
		      for( size_t iPMT = 0; iPMT < calibratedPMTs2.GetCount(); iPMT++ ){
			const RAT::DS::PMTCal& pmtCal = calibratedPMTs2.GetPMT( iPMT );
			lightPath.CalcByPosition( ev_next_vecR, pmtInfo.GetPosition( pmtCal.GetID() ) );
			double distInInnerAV = lightPath.GetDistInInnerAV();
			double distInAV = lightPath.GetDistInAV();
			double distInWater = lightPath.GetDistInWater();
			  
			const double transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater ); // Assumes a 400nm photon
			// Time residuals estimate the photon emission time relative to the event start so subtract off the transit time
			// hit times are relative to the trigger time, which will depend on event time and detector position so correct for that to line up events
			// The 390ns corrects for the electronics delays and places the pulse in the middle of the window
			if (j+k == 1 && j == 0){
			  Late_index0_1TimeRes.Fill( pmtCal.GetTime() - transitTime - 390 + rds.GetMCEV(j+k).GetGTTime());
			}
			if (j+k == 2 && j == 0){
			  Late_index0_2TimeRes.Fill( pmtCal.GetTime() - transitTime - 390 + rds.GetMCEV(j+k).GetGTTime());
			}
		      }
			
		      const RAT::DS::CalPMTs& calibratedPMTs1 = rEV1.GetCalPMTs();
		      for( size_t iPMT = 0; iPMT < calibratedPMTs1.GetCount(); iPMT++ ){
			const RAT::DS::PMTCal& pmtCal = calibratedPMTs1.GetPMT( iPMT );
			lightPath.CalcByPosition( ev_vecR, pmtInfo.GetPosition( pmtCal.GetID() ) );
			double distInInnerAV = lightPath.GetDistInInnerAV();
			double distInAV = lightPath.GetDistInAV();
			double distInWater = lightPath.GetDistInWater();
			      
			const double transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater ); // Assumes a 400nm photon
			// Time residuals estimate the photon emission time relative to the event start so subtract off the transit time
			// hit times are relative to the trigger time, which will depend on event time and detector position so correct for that to line up events
			// The 390ns corrects for the electronics delays and places the pulse in the middle of the window
			if (j+k == 1 && j == 0)
			  Prompt_index0_1TimeRes.Fill( pmtCal.GetTime() - transitTime - 390 + rds.GetMCEV(j).GetGTTime());
			if (j+k == 2 && j == 0)
			  Prompt_index0_2TimeRes.Fill( pmtCal.GetTime() - transitTime - 390 + rds.GetMCEV(j).GetGTTime());
		      }
			    
		      if (j+k == 2 && j == 0){
			
			const RAT::DS::CalPMTs& calibratedPMTs = rds.GetEV(1).GetCalPMTs();
			for( size_t iPMT = 0; iPMT < calibratedPMTs.GetCount(); iPMT++ ){
			  const RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetPMT( iPMT );
			  lightPath.CalcByPosition( ev_vecR, pmtInfo.GetPosition( pmtCal.GetID() ) );
			  double distInInnerAV = lightPath.GetDistInInnerAV();
			  double distInAV = lightPath.GetDistInAV();
			  double distInWater = lightPath.GetDistInWater();
			  double transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater ); // Assumes a 400nm photon
			  BadEvPos1in1_index0_2TimeRes.Fill( pmtCal.GetTime() - transitTime - 390 + rds.GetMCEV(j).GetGTTime());
			    
			  lightPath.CalcByPosition( ev_next_vecR, pmtInfo.GetPosition( pmtCal.GetID() ) );
			  distInInnerAV = lightPath.GetDistInInnerAV();
			  distInAV = lightPath.GetDistInAV();
			  distInWater = lightPath.GetDistInWater();
			  transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater ); // Assumes a 400nm photon
			  BadEvPos2in1_index0_2TimeRes.Fill( pmtCal.GetTime() - transitTime - 390 + rds.GetMCEV(j).GetGTTime());
			}
		      }
		      //h_after_cut_emc_nu.Fill(rmcparent.GetKineticEnergy());
		      //h_after_cut_emc_p1.Fill(rmcparticle1.GetKineticEnergy());
		      //h_after_cut_emc_p2.Fill(rmcparticle2.GetKineticEnergy());
		      //h2_after_cut_emc_p2_vs_p1.Fill(rmcparticle1.GetKineticEnergy(), rmcparticle2.GetKineticEnergy());

		      //std::cout<<"\nEdep: "<<edep<<std::endl;
		      //std::cout<<"EQuench: "<<equench<<std::endl;
		      //std::cout<<"neuton ke: "<<part1ke<<std::endl;
		      //std::cout<<"alpha ke: "<<part2ke<<std::endl;
		      h_edep.Fill(edep);
		      h_equench.Fill(equench);
		      if (is_neutron){
			h_part1ke.Fill(part1ke);
		      }
		      if (is_alpha){
			h_part2ke.Fill(part2ke);
		      }
		      if (is_alpha && is_neutron)
			h_ke1vske2.Fill(part2ke, part1ke);
		      
		      h_after_cut_deltaT.Fill(timeD);
		      h_after_cut_deltaR.Fill(deltaR);
		      h2_after_cut_deltaT_vs_deltaR.Fill(deltaR, timeD);
		      h_after_cut_efit_prompt.Fill(ev_energy);
		      h_after_cut_efit_delayed.Fill(ev_next_energy);
		      h2_after_cut_efit_delayed_vs_prompt.Fill(ev_energy, ev_next_energy);


		      h_prompt_berkeleyAlphaBeta.Fill(ev_berkeleyAB_lh);
		      h_prompt_ITR.Fill(ev_itr);
		      h_prompt_Beta14.Fill(ev_beta14);

		      h_late_berkeleyAlphaBeta.Fill(ev_next_berkeleyAB_lh);
		      h_late_ITR.Fill(ev_next_itr);
		      h_late_Beta14.Fill(ev_next_beta14);
//delE_efit_prompt.Fill(ev_energy,rmcparticle1.GetKineticEnergy() + annihilation - ev_energy);
		      //delE_emc_p1.Fill(rmcparticle1.GetKineticEnergy(),rmcparticle1.GetKineticEnergy() + annihilation - ev_energy);

		      if (j+k == 2 && j == 0){
			ev02 += 1;
			h_prompt_berkeleyAlphaBeta_02.Fill(ev_berkeleyAB_lh);
			h_prompt_ITR_02.Fill(ev_itr);
			h_prompt_Beta14_02.Fill(ev_beta14);
			h_late_berkeleyAlphaBeta_02.Fill(ev_next_berkeleyAB_lh);
			h_late_ITR_02.Fill(ev_next_itr);
			h_late_Beta14_02.Fill(ev_next_beta14);

			h_edep_02.Fill(edep);
			h_equench_02.Fill(equench);
			if (is_neutron){
			  h_part1ke_02.Fill(part1ke);
			}
			if (is_alpha){
			  h_part2ke_02.Fill(part2ke);
			}
			if (is_alpha && is_neutron)
			  h_ke1vske2_02.Fill(part2ke, part1ke);

		      }
		      if (j+k == 1 && j == 0){
			ev01 += 1;
			h_prompt_berkeleyAlphaBeta_01.Fill(ev_berkeleyAB_lh);
			h_prompt_ITR_01.Fill(ev_itr);
			h_prompt_Beta14_01.Fill(ev_beta14);
			h_late_berkeleyAlphaBeta_01.Fill(ev_next_berkeleyAB_lh);
			h_late_ITR_01.Fill(ev_next_itr);
			h_late_Beta14_01.Fill(ev_next_beta14);

			h_edep_01.Fill(edep);
			h_equench_01.Fill(equench);
			if (is_neutron){
			  h_part1ke_01.Fill(part1ke);
			}
			if (is_alpha){
			  h_part2ke_01.Fill(part2ke);
			}
			if (is_alpha && is_neutron)
			  h_ke1vske2_01.Fill(part2ke, part1ke);
		      }
	
		      goodpair = true;  //pair passed quality cuts, if not continue search in group
		      pairfound = true; // a pair was found
		      k += 100; // cancel sub search
		    }
		  }
		}
	      }
	    }
	  }
	}// if pair failed quality cuts, continue sub group search
	if (!goodpair){
	  k += 1;
	  if (Fitted <= 0){ //if positron/j_index not a valid fit, move onto next sub group
	    k += 100;
	  }else{ //else if good, check that it satisfies time diff, if not move onto next sub group
	    time_diff_condition = timeD;
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

  TFile fout(filename_output.c_str(), "RECREATE");

  std::cout<<"evindex0_1: "<<ev01<<std::endl;
  std::cout<<"evindex0_2: "<<ev02<<std::endl;
  
  std::cout<<"EV partner events within deltaT: "<<numtagged<<std::endl;
  //std::cout <<  " tot events: " << evcounttot << std::endl;

  h_edep.Write();
  h_edep_01.Write();
  h_edep_02.Write();
  h_equench.Write();
  h_equench_01.Write();
  h_equench_02.Write();
  h_part1ke.Write();
  h_part1ke_01.Write();
  h_part1ke_02.Write();
  h_part2ke.Write();
  h_part2ke_01.Write();
  h_part2ke_02.Write();
  h_ke1vske2.Write();
  h_ke1vske2_01.Write();
  h_ke1vske2_02.Write();

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
  /*delE_emc_p1_0_1.Write();
  delE_emc_p1_0_2.Write();
  c_delekeplus_e1->Write();
  c_delekeplus_poske->Write();
  c_delekeplus_poske_1v2->Write();
  */
  deltaTimeBadEVindex1.Write();
  deltaRBadEVindex1.Write();


  h_prompt_berkeleyAlphaBeta.Write();
  h_prompt_berkeleyAlphaBeta_01.Scale(1./(double)h_prompt_berkeleyAlphaBeta_01.Integral());
  h_prompt_berkeleyAlphaBeta_01.Write();
  h_prompt_berkeleyAlphaBeta_02.Scale(1./(double)h_prompt_berkeleyAlphaBeta_02.Integral());
  h_prompt_berkeleyAlphaBeta_02.Write();
  h_late_berkeleyAlphaBeta_01.Scale(1./(double)h_late_berkeleyAlphaBeta_01.Integral());
  h_late_berkeleyAlphaBeta_01.Write();
  h_late_berkeleyAlphaBeta_02.Scale(1./(double)h_late_berkeleyAlphaBeta_02.Integral());
  h_late_berkeleyAlphaBeta_02.Write();
  h_prompt_ITR.Write();
  h_prompt_ITR_01.Scale(1./(double)h_prompt_ITR_01.Integral());
  h_prompt_ITR_01.Write();
  h_prompt_ITR_02.Scale(1./(double)h_prompt_ITR_02.Integral());
  h_prompt_ITR_02.Write();
  h_late_ITR_01.Scale(1./(double)h_late_ITR_01.Integral());
  h_late_ITR_01.Write();
  h_late_ITR_02.Scale(1./(double)h_late_ITR_02.Integral());
  h_late_ITR_02.Write();
  h_prompt_Beta14.Write();
  h_prompt_Beta14.Write();
  h_prompt_Beta14_01.Scale(1./(double)h_prompt_Beta14_01.Integral());
  h_prompt_Beta14_01.Write();
  h_prompt_Beta14_02.Scale(1./(double)h_prompt_Beta14_02.Integral());
  h_prompt_Beta14_02.Write();
  h_late_Beta14_01.Scale(1./(double)h_late_Beta14_01.Integral());
  h_late_Beta14_01.Write();
  h_late_Beta14_02.Scale(1./(double)h_late_Beta14_02.Integral());
  h_late_Beta14_02.Write();
  
  //oscRoot
  
  TCanvas *c_prompt01tres = new TCanvas("c_prompt01tres","time resolution [ns]");
  c_prompt01tres->SetLogy();
  Prompt_index0_1TimeRes.Scale(1./(double)Prompt_index0_1TimeRes.Integral());
  Prompt_index0_1TimeRes.Draw();
  TCanvas *c_late01tres = new TCanvas("c_late01tres","time resolution [ns]");
  c_late01tres->SetLogy();
  Late_index0_1TimeRes.Scale(1./(double)Late_index0_1TimeRes.Integral());
  Late_index0_1TimeRes.Draw();
  TCanvas *c_prompt02tres = new TCanvas("c_prompt02tres","time resolution [ns]");
  c_prompt02tres->SetLogy();
  Prompt_index0_2TimeRes.Scale(1./(double)Prompt_index0_2TimeRes.Integral());
  Prompt_index0_2TimeRes.Draw();
  TCanvas *c_late02tres = new TCanvas("c_late02tres","time resolution [ns]");
  c_late02tres->SetLogy();
  Late_index0_2TimeRes.Scale(1./(double)Late_index0_2TimeRes.Integral());
  Late_index0_2TimeRes.Draw();

  Prompt_index0_1TimeRes.Write();
  Late_index0_1TimeRes.Write();
  Prompt_index0_2TimeRes.Write();
  Late_index0_2TimeRes.Write();

  BadEvPos1in1_index0_2TimeRes.Write();
  BadEvPos2in1_index0_2TimeRes.Write();
  
  c_prompt01tres->Write();
  c_prompt02tres->Write();
  c_late01tres->Write();
  c_late02tres->Write();


  fout.Close();
    
}

 
