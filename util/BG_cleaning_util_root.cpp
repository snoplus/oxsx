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

#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/FitResult.hh>
#include <RAT/DB.hh>
#include <RAT/DU/Utility.hh>

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
  
  TH1D h_edep("h_edep", "Deposited Energy (MeV)", 300, 0, 9);
  TH1D h_equench("h_equench", "Quenched Energy (MeV)", 300, 0, 4);
  TH1D h_part1ke("h_part1ke", "Particle 1 KE (MeV)", 300, 0, 9);
  //TH1D h_part1ke("h_part1ke", "Particle 1 KE (MeV)", 300, 0, 9);
  
  TH1D h_after_cut_R("h_after_cut_R", "Particle Radius (mm)", 300, 0, 2);
  
  TH1D h_after_cut_efit("h_after_cut_efit", "Reconstructed Energy (MeV)", 300, 0, 9);

  TH1D h_after_cut_R_0("h_after_cut_R_0", "Particle Radius (mm)", 300, 0, 2);  
  TH1D h_after_cut_efit_0("h_after_cut_efit_0", "Reconstructed Energy (MeV)", 300, 0, 9);
  TH1D h_after_cut_R_1("h_after_cut_R_1", "Particle Radius (mm)", 300, 0, 2);  
  TH1D h_after_cut_efit_1("h_after_cut_efit_1", "Reconstructed Energy (MeV)", 300, 0, 9);
  TH1D h_after_cut_R_2("h_after_cut_R_2", "Particle Radius (mm)", 300, 0, 2);  
  TH1D h_after_cut_efit_2("h_after_cut_efit_2", "Reconstructed Energy (MeV)", 300, 0, 9);

  double classifier_del_lhmax = 200;
  double classifier_del_lhmin = -200.;
  size_t num_bins_classifier = 400;

  size_t time_res_bins = 850;
  double h_time_res_min = -250;
  double h_time_res_max = 600;

  TH1D *h_TimeRes = new TH1D("TimeRes","TimeRes",time_res_bins,h_time_res_min, \
			     h_time_res_max);
  TH1D h_berkeleyAlphaBeta("h_berkeleyAlphaBeta", "BerkeleyAlphaBeta", num_bins_classifier, classifier_del_lhmin, classifier_del_lhmax);

  TH1D *h_TimeRes_0 = new TH1D("TimeRes_0","TimeRes",time_res_bins,h_time_res_min,h_time_res_max);
  TH1D h_berkeleyAlphaBeta_0("h_berkeleyAlphaBeta_0", "BerkeleyAlphaBeta", num_bins_classifier, classifier_del_lhmin, classifier_del_lhmax);
  TH1D *h_TimeRes_1 = new TH1D("TimeRes_1","TimeRes",time_res_bins,h_time_res_min, \
			       h_time_res_max);
  TH1D h_berkeleyAlphaBeta_1("h_berkeleyAlphaBeta_1", "BerkeleyAlphaBeta", num_bins_classifier, classifier_del_lhmin, classifier_del_lhmax);
  TH1D *h_TimeRes_2 = new TH1D("TimeRes_2","TimeRes",time_res_bins,h_time_res_min, \
			       h_time_res_max);
  TH1D h_berkeleyAlphaBeta_2("h_berkeleyAlphaBeta_2", "BerkeleyAlphaBeta", num_bins_classifier, classifier_del_lhmin, classifier_del_lhmax);


  TH1D h_alphaBeta212("h_alphaBeta212", "AlphaBeta212", num_bins_classifier, classifier_del_lhmin, classifier_del_lhmax);
  TH1D h_alphaBeta214("h_alphaBeta214", "AlphaBeta214", num_bins_classifier, classifier_del_lhmin, classifier_del_lhmax);

  TH1D h_ITR("h_ITR", "ITR", 300, -1.5, 1.5  );
  TH1D h_Beta14("h_Beta14", "Beta14", 400, -1, 3);
  
  // load input file
  RAT::DU::DSReader dsreader(filename_input.c_str());
    

  Double_t neutron_capture_energy = 2.2;//1.857;
  Double_t e_rem = 0.784;
  Double_t annihilation = 1.022;

  ULong64_t numtagged = 0;
    

  RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path
  const RAT::DU::GroupVelocity& groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity(); // To get the group velocity
  const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo(); // The PMT positions etc..

  
  ///////////////////////////////////////////////
  /// Cleaning+Coinc. Tagging (using evindex) ///
  ///////////////////////////////////////////////

  size_t n_entries = dsreader.GetEntryCount();
  //n_entries = 500;
  for(size_t i = 0; i < n_entries; i++){
    if (i % 100 == 0) std::cout <<  "    Processed entries: " << i << " (" << (double)i/n_entries*100. << "%) " << std::endl;

    const RAT::DS::Entry& rds = dsreader.GetEntry(i);

    const RAT::DS::MC& rMC = rds.GetMC();
    //size_t ParentCount = rMC.GetMCParentCount();
    //size_t ParticleCount = rMC.GetMCParticleCount();
    //const RAT::DS::MCParticle& rmcparent = rMC.GetMCParent(0);
    const RAT::DS::MCParticle& rmcparticle1 = rMC.GetMCParticle(0);
    double part1ke = rmcparticle1.GetKineticEnergy();
    double equench = rMC.GetScintQuenchedEnergyDeposit();
    double edep = rMC.GetScintEnergyDeposit();
    //int parentpdg = rmcparent.GetPDGCode();
    //int particlepdg = rmcparticle.GetPDGCode();
    //std::cout<<"parent pdg: "<<parentpdg<<" particle pdg: "<<particlepdg<<std::endl;
    
    for( int iEV = 0; iEV < rds.GetEVCount(); iEV++ ){
      const RAT::DS::EV& rEV = rds.GetEV( iEV );

      //std::cout<<"\nEntry: "<<j<<std::endl;
      if( rEV.FitResultExists("scintFitter") == 0 ) 
	continue;
	  
      if( !rEV.GetFitResult("scintFitter").GetVertex(0).ContainsPosition() ) continue;
      if( !rEV.GetFitResult("scintFitter").GetVertex(0).ValidPosition() ) continue;
	  
      RAT::DS::FitResult fResult = rEV.GetFitResult("scintFitter");
      RAT::DS::FitVertex fVertex = fResult.GetVertex(0); 
      TVector3 fPosition = fVertex.GetPosition();
      double fEnergy = fVertex.GetEnergy();
      double Z = fPosition.Z();
      Z = Z - 108;
      fPosition = TVector3(fPosition.X(),fPosition.Y(),Z);

      const RAT::DS::ClassifierResult &result = rEV.GetClassifierResult("BerkeleyAlphaBeta:scintFitter");
      double likelihood = result.GetClassification("likelihood");

      const RAT::DS::ClassifierResult &result_itr = rEV.GetClassifierResult("ITR:scintFitter");
      double ev_itr = result_itr.GetClassification("ITR");
      //const RAT::DS::ClassifierResult &result_beta14 = rEV.GetClassifierResult("Beta14:scintFitter");
      double ev_beta14 = 0.;//result.GetClassification("isotropy");
      if (fPosition.Mag() < Rmax){
	//if want to correc energies
	//double corrected_ev_energy = ev_energy + correction(ev_energy);
	//if (energy_ep_min <= corrected_ev_energy && corrected_ev_energy <= energy_ep_max){
	//double corrected_ev_next_energy = nextEnergy + correction(nextEnergy);
	//if (energy_n_min <= corrected_ev_next_energy && corrected_ev_next_energy <= energy_n_max){
	if (energy_min <= fEnergy && fEnergy <= energy_max){
	  /*if (iEV == 0)
	    std::cout<<"\nev index: "<<iEV<<"  energy: "<<fEnergy<<std::endl;
	  else
	    std::cout<<"ev index: "<<iEV<<"  energy: "<<fEnergy<<std::endl;
	  */
	  numtagged += 1;

	  h_edep.Fill(edep);
	  h_equench.Fill(equench);
	  h_part1ke.Fill(part1ke);
	  
	  h_after_cut_efit.Fill(fEnergy);
	  h_after_cut_R.Fill(pow(fPosition.Mag()/6000.,3));
	      
	  h_berkeleyAlphaBeta.Fill(likelihood);
	  h_ITR.Fill(ev_itr);
	  h_Beta14.Fill(ev_beta14);

	  //h_alphaBeta212.Fill(ev_alphaBeta212);
	  //h_alphaBeta214.Fill(ev_alphaBeta214);
	  if (iEV == 0){
	    h_berkeleyAlphaBeta_0.Fill(likelihood);
	    h_after_cut_efit_0.Fill(fEnergy);
	    h_after_cut_R_0.Fill(pow(fPosition.Mag()/6000.,3));
	  }if (iEV == 1){
	    h_after_cut_efit_1.Fill(fEnergy);
	    h_after_cut_R_1.Fill(pow(fPosition.Mag()/6000.,3));
	    h_berkeleyAlphaBeta_1.Fill(likelihood);
	  }if (iEV == 2){
	    h_after_cut_efit_2.Fill(fEnergy);
	    h_after_cut_R_2.Fill(pow(fPosition.Mag()/6000.,3));
	    h_berkeleyAlphaBeta_2.Fill(likelihood);
	  }
	  const RAT::DS::CalPMTs& calibratedPMTs = rEV.GetPartialCalPMTs(1); //what 1 mean??
	  for( size_t iPMT = 0; iPMT < calibratedPMTs.GetCount(); iPMT++ ){
	    const RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetPMT( iPMT );
	    TVector3 pmt_pos = pmtInfo.GetPosition( pmtCal.GetID() );
	    lightPath.CalcByPosition( fPosition, pmt_pos);
	    double distInInnerAV = lightPath.GetDistInInnerAV();
	    double distInAV = lightPath.GetDistInAV();
	    double distInWater = lightPath.GetDistInWater();
	    double transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater ); // Assumes a 400nm photon
		
	    double timeresidual = pmtCal.GetTime() - transitTime - fVertex.GetTime();
	    h_TimeRes->Fill(timeresidual);
	    if (iEV == 0)
	      h_TimeRes_0->Fill(timeresidual);
	    if (iEV == 1)
	      h_TimeRes_1->Fill(timeresidual);
	    if (iEV == 2)
	      h_TimeRes_2->Fill(timeresidual);
	  }
	}
      }
    }
  }

  // finished tagging
    
  //summary stats
  std::cout<<"\n \n Summary: "<<std::endl;
  std::cout<<"MC events tagged: "<<numtagged<<std::endl;

  //std::stringstream bozo;
  //bozo<<fileout<<"_N"<<Nmin<<"_"<<Nmax<<"_Z"<<Zmin<<"_"<<Zmax<<"_R"<<Rmin<<"_"<<Rmax<<".root";
  //std::cout<<bozo.str().c_str()<<std::endl;
  TFile *file_output = new TFile( filename_output.c_str(),"RECREATE");

  h_edep.Write();
  h_equench.Write();
  h_part1ke.Write();
  
  h_after_cut_efit.Scale(1./(double)h_after_cut_efit.Integral());
  h_after_cut_R.Scale(1./(double)h_after_cut_R.Integral());    
    
  h_berkeleyAlphaBeta.Scale(1./(double)h_berkeleyAlphaBeta.Integral());
  h_alphaBeta212.Scale(1./(double)h_alphaBeta212.Integral());
  h_alphaBeta214.Scale(1./(double)h_alphaBeta214.Integral());
  
  h_ITR.Scale(1./(double)h_ITR.Integral());
  h_Beta14.Scale(1./(double)h_Beta14.Integral());
  h_ITR.Write();
  h_Beta14.Write();
  
  h_after_cut_R.Write();
  h_after_cut_R_0.Write();
  h_after_cut_R_1.Write();
  h_after_cut_R_2.Write();
  h_after_cut_efit.Write();
  h_after_cut_efit_0.Write();
  h_after_cut_efit_1.Write();
  h_after_cut_efit_2.Write();

  h_berkeleyAlphaBeta.Write();
  h_berkeleyAlphaBeta_0.Write();
  h_berkeleyAlphaBeta_1.Write();
  h_berkeleyAlphaBeta_2.Write();
  h_alphaBeta212.Write();
  h_alphaBeta214.Write();

  TCanvas *c_tres = new TCanvas("c_tres","time resolution [ns]");
  c_tres->SetLogy();
  h_TimeRes->Scale(1./(double)h_TimeRes->Integral());
  h_TimeRes->Draw();
  h_TimeRes->Write();
  TCanvas *c_tres_0 = new TCanvas("c_tres_0","time resolution [ns]");
  c_tres_0->SetLogy();
  h_TimeRes_0->Scale(1./(double)h_TimeRes_0->Integral());
  h_TimeRes_0->Draw();
  h_TimeRes_0->Write();
  TCanvas *c_tres_1 = new TCanvas("c_tres_1","time resolution [ns]");
  c_tres_1->SetLogy();
  h_TimeRes_1->Scale(1./(double)h_TimeRes_1->Integral());
  h_TimeRes_1->Draw();
  h_TimeRes_1->Write();
  TCanvas *c_tres_2 = new TCanvas("c_tres_2","time resolution [ns]");
  c_tres_2->SetLogy();
  h_TimeRes_2->Scale(1./(double)h_TimeRes_2->Integral());
  h_TimeRes_2->Draw();
  h_TimeRes_2->Write();
 
  c_tres->Write();
  c_tres_0->Write();
  c_tres_1->Write();
  c_tres_2->Write();

  std::cout<<"fin: "<<filename_input<<" \n fout: "<<filename_output<<"\n emin "<<energy_min<<"\n emax "<<energy_max<<"\n Rmax "<<Rmax<<std::endl;
  
  file_output->Close();
}
