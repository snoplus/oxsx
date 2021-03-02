// A fit in energy for signal and a background
#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include <Rand.h>
#include <fstream>
#include <iostream>

#include <TCanvas.h>
#include <ROOTNtuple.h>
#include <TRandom3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <THStack.h>

#include <TLatex.h>

#include <BinnedED.h>
#include <BinnedEDGenerator.h>
#include <SystematicManager.h>
#include <BinnedNLLH.h>
#include <FitResult.h>
#include <Minuit.h>
#include <DistTools.h>
#include <Minuit.h>
#include <Convolution.h>
#include <Scale.h>
#include <BoolCut.h>
#include <BoxCut.h>
#include <Gaussian.h>
#include <ParameterDict.h>
#include <ContainerTools.hpp>
#include <NuOsc.h>
#include <SurvProb.h>
#include "AntinuUtils.cpp"
#include "../util/oscillate_util.cpp"

#include <Event.h>
#include <TVector3.h>

Double_t CalculateDistance(TVector3 point1, TVector3 point2) {
    return (point2 - point1).Mag();
}

TVector3 LLAtoECEF(Double_t latitude, Double_t longitude, Double_t altitude) {
    // reference http://www.mathworks.co.uk/help/aeroblks/llatoecefposition.html
    const Double_t toRad = TMath::Pi()/180.;
    const Double_t Earthradius = 6378137.0; //Radius of the Earth (in meters)
    const Double_t f_f = 1./298.257223563; //Flattening factor WGS84 Model
    Double_t x, y, z;
    latitude = latitude*toRad;
    longitude = longitude*toRad;
    const Double_t f_l = TMath::ATan(TMath::Power((1. - f_f),2)*TMath::Tan(latitude))/toRad;
    const Double_t f_rs = TMath::Sqrt( TMath::Power(Earthradius,2)/(1. + (1./TMath::Power((1. - f_f),2) - 1.)*TMath::Power(TMath::Sin(f_l*toRad),2)));
    x = (f_rs*TMath::Cos(f_l*toRad)*TMath::Cos(longitude) + altitude*TMath::Cos(latitude)*TMath::Cos(longitude))/1000.; // in km 
    y = (f_rs*TMath::Cos(f_l*toRad)*TMath::Sin(longitude) + altitude*TMath::Cos(latitude)*TMath::Sin(longitude))/1000.; // in km 
    z = (f_rs*TMath::Sin(f_l*toRad) + altitude*TMath::Sin(latitude))/1000.; // in km
    TVector3 ECEF(x,y,z);
    return ECEF;
}

Double_t GetReactorDistanceLLA(Double_t latitude, Double_t longitude, Double_t altitude, TVector3 ECEF_coord) {
    return CalculateDistance(ECEF_coord, LLAtoECEF(latitude, longitude, altitude));
} 

Double_t reactor_distance_to_evs(double distance, const double P_over_D_coeff, const double avg_total_therm_power){
  double estimate_evs_yr = P_over_D_coeff*avg_total_therm_power/(distance*distance);
  return estimate_evs_yr;
}

/*Double_t nd_evs_integral_to_power(double int_evs, size_t n_detectors, std::vector<double> distance, std::vector<double> osc_loss, const double P_over_D_coeff){
  //calc in steps
  double avg_total_therm_power = int_evs/(double)pow(P_over_D_coeff, n_detectors);
  for (size_t i = 0; i < n_detectors; i++){
    avg_total_therm_power *= (distance[i]*distance[i]/osc_loss[i]);
  }
  avg_total_therm_power = pow(avg_total_therm_power, 1/(double)n_detectors);
  return avg_total_therm_power;
  }*/
Double_t nd_evs_integral_to_power(double int_evs, size_t n_detectors, std::vector<double> distance, std::vector<double> osc_loss, const double P_over_D_coeff){
  //calc in steps
  double osc_over_dist_sum = 0.;
  for (size_t i = 0; i < n_detectors; i++){
    osc_over_dist_sum += (osc_loss[i]/(distance[i]*distance[i]));
  }
  double denominator = P_over_D_coeff*osc_over_dist_sum;
  double avg_total_therm_power = int_evs/(double)denominator;
  
  return avg_total_therm_power;
}

std::vector<Double_t>
LHFit_fit(const std::string &spectrum_phwr_unosc_filepath,
	  const std::string &spectrum_pwr_unosc_filepath,
	  const std::string &spectrum_geo_uraniumthorium_unosc_filepath,
	  const std::string &spectrum_bkg_alphan_unosc_filepath,
	  const std::vector<std::string> &detector_names,
	  const std::vector<std::string> &reactor_names, std::vector<std::string> &reactor_types,
	  std::vector<std::vector<Double_t> > &distances,
      std::vector<std::vector<std::vector<Double_t> > > &distances_fit,
      std::vector<Double_t> &known_reactor_constraints, /*std::vector<Double_t> &constraint_sigmas,*/
	  TFile *file_out,
	  Double_t param_d21, Double_t param_s12, Double_t param_s13,
	  bool &fit_validity,
	  const double e_min, const double e_max, const size_t n_bins,
	  const double flux_data, const double mc_scale_factor,
	  double &fit_geo_uth_norm,
	  const double param_d21_plot_min, const double param_d21_plot_max, const double param_s12_plot_min, 
	  const double param_s12_plot_max, const double P_over_D_coeff, const std::vector<double> thermal_powers,
	  const std::string &data_path, THStack *hs_e1, THStack *hs_e2, size_t n_fit,
	  TH2D &h2_data, TH2D * h2_best_fit, const bool poisson_fluc,
	  TH1D &h_data_set_detector_1, TH1D &h_data_set_detector_2,
          TH1D &h_data_set_detector_1_poisson_fluc, TH1D &h_data_set_detector_2_poisson_fluc,
          double n_years, std::vector<double> &lh_values_vec_rot, const int rotationZ_n, int rotationZ_0value_index){

  //////////////////////////////
  ///// Constrained Data  //////
  //////////////////////////////
  
  printf("Beginning Fake data prep--------------------------------------\n");
  printf("Fake Data:: del_21:%.9f, sin2_12:%.7f, sin2_13:%.7f\n", param_d21, param_s12, param_s13);
      
  char name[1000];
  TRandom3 *random_generator = new TRandom3();
  random_generator->SetSeed(0);

  const size_t n_detectors = detector_names.size();
  const size_t n_reactors = reactor_names.size();

  AxisCollection axes;
  std::vector<size_t> indices;
  for (size_t j = 0; j < n_detectors; j++){
    std::stringstream axis_lbl_strm;
    axis_lbl_strm<<"ev_prompt_fit_detector"<<(j+1);
    axes.AddAxis(BinAxis(axis_lbl_strm.str().c_str(), e_min, e_max, n_bins));
    indices.push_back(j);
  }
      
  ObsSet data_rep(indices);
  BinnedED data_set_pdf("data_set_pdf", axes);
  data_set_pdf.SetObservables(data_rep);
  data_set_pdf.Empty();

  ObsSet data_rep_1d(0);
  AxisCollection axes_1d;
  axes_1d.AddAxis(BinAxis("ev_prompt_fit", e_min, e_max, n_bins));
  
  
  /////////////////
  // 2d BinnedEDs//
  /////////////////
  BinnedED **reactor_osc_data_nd_pdf = new BinnedED*[n_reactors];
  // distances: { [distances to 1] ,  [distances to 2] }
  BinnedED*** reactor_unosc_data_1d_pdf = new BinnedED**[n_detectors];   
  BinnedED*** reactor_osc_data_1d_pdf = new BinnedED**[n_detectors];
  for (ULong64_t j = 0; j < n_detectors; j++){
    reactor_unosc_data_1d_pdf[j] = new BinnedED*[n_reactors];
    reactor_osc_data_1d_pdf[j] = new BinnedED*[n_reactors];
  }
  
  BinnedED **reactor_osc_total_data_axes_pdf = new BinnedED*[n_detectors];
  for (ULong64_t j = 0; j < n_detectors; j++){
    sprintf(name, "osc_total_data_e%s", detector_names[j].c_str());
    reactor_osc_total_data_axes_pdf[j] = new BinnedED(name, axes_1d);
    reactor_osc_total_data_axes_pdf[j]->SetObservables(0);
  }

  TH1D*** reactor_osc_data_1d_hist = new TH1D**[n_detectors];
  for (ULong64_t j = 0; j < n_detectors; j++)
    reactor_osc_data_1d_hist[j] = new TH1D*[n_reactors];

  bool geos_included = false;
  bool alphans_included = false;

  // read in constraint information
  std::vector<double> total_unosc_constraint_means;

  for (ULong64_t i = 0; i < n_reactors; i++){
    Double_t nd_reactor_Int_fit = 0.;
    for (ULong64_t j = 0; j < n_detectors; j++){
      sprintf(name, "%s_%s_unosc", detector_names[j].c_str(),  reactor_names[i].c_str());
      reactor_unosc_data_1d_pdf[j][i] = new BinnedED(name, axes_1d);
      reactor_unosc_data_1d_pdf[j][i]->SetObservables(0);
      reactor_osc_data_1d_pdf[j][i] = new BinnedED(reactor_names[i], axes_1d);
      reactor_osc_data_1d_pdf[j][i]->SetObservables(0);

      bool apply_oscillation = false;
      bool is_further_reactors = false;
      if ((reactor_types[i]=="PWR")||(reactor_types[i]=="BWR")){
        sprintf(name, "%s", spectrum_pwr_unosc_filepath.c_str());
	    apply_oscillation = true;
      }else if (reactor_types[i]=="PHWR"){
        sprintf(name, "%s", spectrum_phwr_unosc_filepath.c_str());
	    apply_oscillation = true;
      }else if (reactor_types[i]=="further_reactors"){
	    sprintf(name, "%s", spectrum_pwr_unosc_filepath.c_str());
	    is_further_reactors = true;
      }else if (reactor_names[i]=="uraniumthorium"){
        sprintf(name, "%s", spectrum_geo_uraniumthorium_unosc_filepath.c_str());
	    geos_included = true;
      }else if (reactor_names[i]=="alphan"){
	    sprintf(name, "%s", spectrum_bkg_alphan_unosc_filepath.c_str());
	    alphans_included = true;
      }else{
        printf("Throw: Reactor doesn't match any loaded type...\n");
        exit(0); // throw std::exception(); //continue;
      }
		  
      TFile *f_in = new TFile(name);
      file_out->cd();
      TTree *reactor_unosc_ntp = (TTree*)f_in->Get("nt");
      TNtuple *reactor_osc_ntp = new TNtuple("nt", "Oscillated Prompt Energy", "ev_fit_energy_p1");
	
      // oscillate tree
      if (apply_oscillation){
	    ntOscillate_pruned(reactor_unosc_ntp, reactor_osc_ntp, param_d21, param_s12, param_s13, distances[j][i]);
        //std::cout<<"osc: "<<detector_names[j]<<" -> "<<reactor_names[i]<<" "<<distances[j][i]<<std::endl;
      }else if (is_further_reactors)
	    ntOscillate_pruned_geo(reactor_unosc_ntp, reactor_osc_ntp, param_s12);
      else
	    ntNoOscillate_pruned(reactor_unosc_ntp, reactor_osc_ntp);

      // reset branch addresses after oscillating in function (otherwise crash before setting again below..)
      reactor_unosc_ntp->SetBranchStatus("*", 0);
      reactor_unosc_ntp->SetBranchStatus("ev_fit_energy_p1", 1); // (re-enable all branches in use)
      // fill unoscillated pdf
      Double_t ev_unosc_energy_p1;
      reactor_unosc_ntp->SetBranchAddress("ev_fit_energy_p1", &ev_unosc_energy_p1);
      for(size_t e = 0; e < reactor_unosc_ntp->GetEntries(); e++){
        reactor_unosc_ntp->GetEntry(e);
        reactor_unosc_data_1d_pdf[j][i]->Fill(ev_unosc_energy_p1);
      }

      // fill oscillated pdf
      Float_t ev_osc_energy_p1;
      reactor_osc_ntp->SetBranchAddress("ev_fit_energy_p1", &ev_osc_energy_p1);
      for(size_t e = 0; e < reactor_osc_ntp->GetEntries(); e++){
        reactor_osc_ntp->GetEntry(e);
        reactor_osc_data_1d_pdf[j][i]->Fill(ev_osc_energy_p1);
      }

      // close unoscillated reactor file
      f_in->Close();

      if (apply_oscillation || is_further_reactors){
        // work out total oscillated integral of constraints
        int normalisation_unosc = reactor_unosc_data_1d_pdf[j][i]->Integral();
        int normalisation_reactor = reactor_osc_data_1d_pdf[j][i]->Integral();
        Double_t osc_loss = normalisation_reactor/(double)normalisation_unosc;

        Double_t constraint_osc_mean = reactor_distance_to_evs(distances[j][i],P_over_D_coeff,thermal_powers[i])*osc_loss*n_years;

        if (j == 0)
          total_unosc_constraint_means.push_back(constraint_osc_mean);
        else if (j == 1)
          total_unosc_constraint_means[i] += constraint_osc_mean;

        nd_reactor_Int_fit += constraint_osc_mean;

        reactor_osc_data_1d_pdf[j][i]->Normalise(); 
        reactor_osc_data_1d_pdf[j][i]->Scale(constraint_osc_mean);

        reactor_osc_total_data_axes_pdf[j]->Add(*reactor_osc_data_1d_pdf[j][i]);

        if (n_fit == 1)
          reactor_osc_data_1d_hist[j][i] = new TH1D(DistTools::ToTH1D(*reactor_osc_data_1d_pdf[j][i]));
	  
        //data_set_pdf.Add(*reactor_osc_data_1d_pdf[j][i]);
        //printf("  added reactor %d/%d: %s, osc_survival: %.3f, norm_constraint: %.3f data_int:%.0f\n", i+1, n_reactors, reactor_names[i].c_str(), osc_loss,constraint_osc_mean, data_set_pdf.Integral());
        printf(" detector %d/%d: %s, added reactor %d/%d: %s, %u/%u -> osc_survival: %.3f, norm_constraint: %.3f dist: %.3f\n", j+1, n_detectors, detector_names[j].c_str(), i+1, n_reactors, reactor_names[i].c_str(), normalisation_reactor, normalisation_unosc, osc_loss,constraint_osc_mean,distances[j][i]);
        //std::cout<<"i: "<<i+1<<" j: "<<j+1<<" reactor: "<<reactor_names[i]<<"  reactor_osc_total_data_axes_pdf int: "<<reactor_osc_total_data_axes_pdf[j]->Integral()<<std::endl;

      }else if (geos_included){
        Double_t normalisation_unosc = reactor_unosc_data_1d_pdf[j][i]->Integral();
        Double_t normalisation_reactor = reactor_osc_data_1d_pdf[j][i]->Integral();
        Double_t osc_loss = normalisation_reactor/normalisation_unosc;

        double scale = 0.;
	  
        // geo custom flux = coefficient times 1 years predicted geos
        reactor_osc_data_1d_pdf[j][i]->Scale(0.);//geo_uth_custom_flux/(double)flux_data);

        //data_set_pdf.Add(*reactor_osc_data_1d_pdf[j][i]);
        //printf("  added geo %d/%d: %s, osc_survival: %.3f, data_int:%.0f\n", i+1, n_reactors, reactor_names[i].c_str(), osc_loss, data_set_pdf.Integral());
	
      }else if (alphans_included){
        Double_t normalisation_unosc = reactor_unosc_data_1d_pdf[j][i]->Integral();
        Double_t normalisation_reactor = reactor_osc_data_1d_pdf[j][i]->Integral();
        Double_t osc_loss = normalisation_reactor/normalisation_unosc;

        double scale = 0.;
	  
        reactor_osc_data_1d_pdf[j][i]->Normalise(); // currently custom flux number = total events in livetime
        reactor_osc_data_1d_pdf[j][i]->Scale(0.);//bkg_alphan_custom_flux);///(double)flux_data);

        //data_set_pdf.Add(*reactor_osc_data_1d_pdf[j][i]);
        //printf("  added alpha-n %d/%d: %s, osc_survival: %.3f, data_int:%.0f\n", i+1, n_reactors, reactor_names[i].c_str(), osc_loss, data_set_pdf.Integral());
      }
    }
    std::cout<<reactor_names[i]<<" -> "<<nd_reactor_Int_fit<<std::endl;
  }


  // data set integral over all detectors
  Double_t reactor_nd_Int = 0.;

  // each detector's measured spectrum
  std::vector<Double_t> data_integrals_1d(n_detectors, 0.);
  h_data_set_detector_1 = DistTools::ToTH1D(*reactor_osc_total_data_axes_pdf[0]);
  h_data_set_detector_2 = DistTools::ToTH1D(*reactor_osc_total_data_axes_pdf[1]);
  h_data_set_detector_1.SetName("h_data_set_detector_1");
  h_data_set_detector_1.SetTitle("Simulated Data at detector 1");
  h_data_set_detector_2.SetName("h_data_set_detector_2");
  h_data_set_detector_2.SetTitle("Simulated Data at detector 2");

  //poisson fluctuate data:
  double frac_poisson_change[n_detectors][n_bins];
  for (ULong64_t j = 0; j < n_detectors; j++){
    std::cout<<"det: "<<j<<" Int: "<<reactor_osc_total_data_axes_pdf[j]->Integral()<<std::endl;
    if (poisson_fluc){
      for (size_t l = 0; l < n_bins; l++){
        Double_t old_content = reactor_osc_total_data_axes_pdf[j]->GetBinContent(l);
        Double_t new_content = random_generator->Poisson(old_content);
        reactor_osc_total_data_axes_pdf[j]->SetBinContent(l, new_content);
        frac_poisson_change[j][l] = (new_content/old_content);
        //std::cout<<"det: "<<j<<" bin: "<<l<<" "<<old_content<<" -> "<<new_content<<"  "<<frac_poisson_change[j][l]<<std::endl;
      }
    }

    data_integrals_1d[j] = reactor_osc_total_data_axes_pdf[j]->Integral();
    reactor_nd_Int += data_integrals_1d[j];
  }

  h_data_set_detector_1_poisson_fluc = DistTools::ToTH1D(*reactor_osc_total_data_axes_pdf[0]);
  h_data_set_detector_2_poisson_fluc = DistTools::ToTH1D(*reactor_osc_total_data_axes_pdf[1]);
  h_data_set_detector_1_poisson_fluc.SetName("h_data_set_detector_1_poisson_fluc");
  h_data_set_detector_1_poisson_fluc.SetTitle("Simulated Data at detector 1 (Poisson fluctutated bins - Red)");
  h_data_set_detector_2_poisson_fluc.SetName("h_data_set_detector_2_poisson_fluc");
  h_data_set_detector_2_poisson_fluc.SetTitle("Simulated Data at detector 2 (Poisson fluctutated bins - Red)");
  h_data_set_detector_1_poisson_fluc.SetLineColor(2);
  h_data_set_detector_2_poisson_fluc.SetLineColor(2);

  reactor_osc_total_data_axes_pdf[0]->Normalise();
  reactor_osc_total_data_axes_pdf[1]->Normalise();

  ///  Make 2d BinnedED
  // make each reactor's 2d pdf, scale to expected value, then add to data_set_pdf
  for (ULong64_t i = 0; i < n_reactors ; i++){
    reactor_osc_data_nd_pdf[i] = new BinnedED(reactor_names[i], axes);
    reactor_osc_data_nd_pdf[i]->SetObservables(data_rep);
	
    std::vector<Double_t> e1_e2_vec(n_detectors, 0.);
    for (size_t l = 0; l < n_bins; l++){
      for (size_t k = 0; k < n_bins; k++){
        std::vector<Double_t> e1_bin_centre_vals(1, 0.);
        axes_1d.GetBinCentres(l, e1_bin_centre_vals);
        e1_e2_vec[0] = e1_bin_centre_vals[0];
        std::vector<Double_t> e2_bin_centre_vals(1, 0.);
        axes_1d.GetBinCentres(k, e2_bin_centre_vals);
        e1_e2_vec[1] = e2_bin_centre_vals[0];	      
        
        reactor_osc_data_1d_pdf[0][i]->Normalise();
        Double_t bin_content = reactor_osc_data_1d_pdf[0][i]->GetBinContent(l);
      
        for (size_t nd = 1; nd < n_detectors; nd++){
          reactor_osc_data_1d_pdf[nd][i]->Normalise();
          bin_content *= reactor_osc_data_1d_pdf[nd][i]->GetBinContent(k);
        }

        reactor_osc_data_nd_pdf[i]->Fill(e1_e2_vec, bin_content);
        //data_content[l][k] += bin_content*total_unosc_constraint_means[i];
      }
    }
    reactor_osc_data_nd_pdf[i]->Normalise();
    reactor_osc_data_nd_pdf[i]->Scale(total_unosc_constraint_means[i]);
    std::cout<<reactor_names[i]<<" total int: "<<reactor_osc_data_nd_pdf[i]->Integral()<<" contraint: "<<total_unosc_constraint_means[i]<<std::endl;
    //if (!poisson_fluc)
    data_set_pdf.Add(*reactor_osc_data_nd_pdf[i]);
  }
  std::cout<<"\n Data norm: Total number of events across detectors = "<<data_set_pdf.Integral()<<"\n"<<std::endl;  

  // apply previously calculated poisson fluctuations to data
  if (poisson_fluc){
    TH2D h2_data_temp = DistTools::ToTH2D(data_set_pdf);
    h2_data_temp.SetName("h2_data_temp");
    data_set_pdf.Empty(); // empty it then fill
    std::vector<Double_t> e1_e2_vec(n_detectors, 0.);
    for (size_t l = 0; l < n_bins; l++){
      for (size_t k = 0; k < n_bins; k++){
        std::vector<Double_t> e1_bin_centre_vals(1, 0.);
        axes_1d.GetBinCentres(l, e1_bin_centre_vals);
        e1_e2_vec[0] = e1_bin_centre_vals[0];
        std::vector<Double_t> e2_bin_centre_vals(1, 0.);
        axes_1d.GetBinCentres(k, e2_bin_centre_vals);
        e1_e2_vec[1] = e2_bin_centre_vals[0];	      
        
        double old_content = h2_data_temp.GetBinContent(l+1,k+1);
        double new_content = 0.;
        //if (frac_poisson_uncert[0][l] != -1. && frac_poisson_uncert[1][k] != -1.)
        new_content = old_content*frac_poisson_change[0][l]*frac_poisson_change[1][k];
        //std::cout<<"x: "<<l<<" y: "<<k<<" "<<frac_poisson_change[0][l]<<" "<<frac_poisson_change[1][k]<<std::endl;
        //std::cout<<"x: "<<l<<" y: "<<k<<" "<<old_content<<" -> "<<new_content<<std::endl;
        data_set_pdf.Fill(e1_e2_vec,new_content);
      }
    }
  }
  data_set_pdf.Normalise();
  
  std::cout<<"\n Data norm: Total number of events across detectors = "<<reactor_nd_Int<<"\n"<<std::endl;
  data_set_pdf.Scale(reactor_nd_Int);
  
  std::cout<<"\n Data norm: Total number of events across detectors = "<<data_set_pdf.Integral()<<"\n"<<std::endl;  
  

  // Data method old
  /*std::vector<Double_t> e1_e2_vec(n_detectors, 0.);
  for (size_t l = 0; l < n_bins; l++){
    for (size_t k = 0; k < n_bins; k++){
      std::vector<Double_t> e1_bin_centre_vals(1, 0.);
      axes_1d.GetBinCentres(l, e1_bin_centre_vals);
      e1_e2_vec[0] = e1_bin_centre_vals[0];
      std::vector<Double_t> e2_bin_centre_vals(1, 0.);
      axes_1d.GetBinCentres(k, e2_bin_centre_vals);
      e1_e2_vec[1] = e2_bin_centre_vals[0];	      

      Double_t bin_content = reactor_osc_total_data_axes_pdf[0]->GetBinContent(l);
      
      for (size_t nd = 1; nd < n_detectors; nd++)
        bin_content *= reactor_osc_total_data_axes_pdf[nd]->GetBinContent(k);

      //reactor_nd_Int += bin_content;
      //std::cout<<"e1: "<<e1_e2_vec[0]<<" e2: "<<e1_e2_vec[1]<<"  ->  "<<bin_content<<std::endl;
      data_set_pdf.Fill(e1_e2_vec, bin_content);
    }
  }
  data_set_pdf.Normalise();
  
  std::cout<<"\n Data norm: Total number of events across detectors = "<<reactor_nd_Int<<"\n"<<std::endl;
  data_set_pdf.Scale(reactor_nd_Int);
  */

  
  h2_data = DistTools::ToTH2D(data_set_pdf);
  h2_data.SetName("h2_data");

  /*TFile *file_data_out = new TFile(data_path.c_str(), "RECREATE");
  TH2D DataDist = DistTools::ToTH2D(data_set_pdf);
  DataDist.SetName("h_data");
  DataDist.Write();
  TH1D DataDist_e1 = DistTools::ToTH1D(data_set_pdf_e1);
  DataDist_e1.SetName("h_data_e1");
  DataDist_e1.Write();
  TH1D DataDist_e2 = DistTools::ToTH1D(data_set_pdf_e2);
  DataDist_e2.SetName("h_data_e2");
  DataDist_e2.Write();
  Double_t test_2d_Int = test_2d.Integral();
  std::cout<<"test_2d Int: "<<test_2d_Int<<std::endl;
  test_2d.Write();
  file_data_out->Close();
  */


  ////////////////////////////////////
  ////////////  Fitting  /////////////
  ////////////////////////////////////
  //fit_validity = true;
  std::vector<Double_t> reactor_power_resolutions(n_reactors, 0.);
  for (int rot_i = 0; rot_i < rotationZ_n ; rot_i++){
    if (rotationZ_n > 1) std::cout<<"\n Rotation: "<<rot_i<<std::endl;

    printf("\n\nBegin fit--------------------------------------\n");
    printf("LHFit_fit:: del_21:%.9f, sin2_12:%.7f, sin2_13:%.7f\n", param_d21, param_s12, param_s13);

    //ObsSet fit_rep(0);
    // create LH function
    BinnedNLLH lh_function;
    lh_function.SetBufferAsOverflow(true);
    int buff = 1;
    lh_function.SetBuffer(0, buff, buff);
    lh_function.SetDataDist(data_set_pdf); // initialise withe the data set

    // setup max and min ranges
    ParameterDict minima;
    ParameterDict maxima;
    ParameterDict initial_val;
    ParameterDict initial_err;

    BinnedED reactor_osc_pdf_fitosc_sum("reactor_osc_pdf_fitosc_sum",axes);
    reactor_osc_pdf_fitosc_sum.SetObservables(data_rep);

    BinnedED **reactor_osc_nd_pdf = new BinnedED*[n_reactors];
    // distances: { [distances to 1] ,  [distances to 2] }
    BinnedED*** reactor_unosc_1d_pdf = new BinnedED**[n_detectors];   
    BinnedED*** reactor_osc_1d_pdf = new BinnedED**[n_detectors];
    for (ULong64_t j = 0; j < n_detectors; j++){
      reactor_unosc_1d_pdf[j] = new BinnedED*[n_reactors];
      reactor_osc_1d_pdf[j] = new BinnedED*[n_reactors];
    }

    Double_t data_set_pdf_integral = data_set_pdf.Integral();
    std::cout<<" Data Int: "<<data_set_pdf_integral<<std::endl;

    geos_included = false;

    std::vector<std::vector<Double_t> >osc_losses;   // losses: { [reactor A losses to detectors] ,  [reactor B losses to detectors] ..}

    for (ULong64_t i = 0; i < n_reactors; i++){
      std::vector<Double_t> reactor_osc_losses;
      osc_losses.push_back(reactor_osc_losses);

      Double_t nd_reactor_Int_fit = 0.;
      for (ULong64_t j = 0; j < n_detectors; j++){
        sprintf(name, "%s_%s_unosc", detector_names[j].c_str(),  reactor_names[i].c_str());
        reactor_unosc_1d_pdf[j][i] = new BinnedED(name, axes_1d);
        reactor_unosc_1d_pdf[j][i]->SetObservables(0);
        reactor_osc_1d_pdf[j][i] = new BinnedED(reactor_names[i], axes_1d);
        reactor_osc_1d_pdf[j][i]->SetObservables(0);

        bool apply_oscillation = false;
        bool is_further_reactors = false;
        if ((reactor_types[i]=="PWR")||(reactor_types[i]=="BWR")){
          sprintf(name, "%s", spectrum_pwr_unosc_filepath.c_str());
          apply_oscillation = true;
        }else if (reactor_types[i]=="PHWR"){
          sprintf(name, "%s", spectrum_phwr_unosc_filepath.c_str());
          apply_oscillation = true;
        }else if (reactor_types[i]=="further_reactors"){
          sprintf(name, "%s", spectrum_pwr_unosc_filepath.c_str());
          is_further_reactors = true;
        }else if (reactor_names[i]=="uraniumthorium"){
          sprintf(name, "%s", spectrum_geo_uraniumthorium_unosc_filepath.c_str());
          geos_included = true;
        }else if (reactor_names[i]=="alphan"){
          sprintf(name, "%s", spectrum_bkg_alphan_unosc_filepath.c_str());
          alphans_included = true;
        }else{
          printf("Throw: Reactor doesn't match any loaded type...\n");
          exit(0);
        }
      
        TFile *f_in = new TFile(name);
        file_out->cd();
        TTree *reactor_unosc_ntp = (TTree*)f_in->Get("nt");
        TNtuple *reactor_osc_ntp = new TNtuple("nt", "Oscillated Prompt Energy", "ev_fit_energy_p1");

        // oscillate tree
        if (apply_oscillation)
          ntOscillate_pruned(reactor_unosc_ntp, reactor_osc_ntp, param_d21, param_s12, param_s13, distances_fit[rot_i][j][i]);
        else if (is_further_reactors)
          ntOscillate_pruned_geo(reactor_unosc_ntp, reactor_osc_ntp, param_s12);
        else
          ntNoOscillate_pruned(reactor_unosc_ntp, reactor_osc_ntp);

        // reset branch addresses after oscillating in function (otherwise crash before setting again below..)
        reactor_unosc_ntp->SetBranchStatus("*", 0);
        reactor_unosc_ntp->SetBranchStatus("ev_fit_energy_p1", 1); // (re-enable all branches in use)

        // fill unoscillated pdf
        Double_t ev_unosc_energy_p1;
        reactor_unosc_ntp->SetBranchAddress("ev_fit_energy_p1", &ev_unosc_energy_p1);
        for(size_t e = 0; e < reactor_unosc_ntp->GetEntries(); e++){
          reactor_unosc_ntp->GetEntry(e);
          reactor_unosc_1d_pdf[j][i]->Fill(ev_unosc_energy_p1);
        }
      
        // fill oscillated pdf
        Float_t ev_osc_energy_p1;
        reactor_osc_ntp->SetBranchAddress("ev_fit_energy_p1", &ev_osc_energy_p1);
        for(size_t e = 0; e < reactor_osc_ntp->GetEntries(); e++){
          reactor_osc_ntp->GetEntry(e);
          reactor_osc_1d_pdf[j][i]->Fill(ev_osc_energy_p1);
        }

        // close unoscillated reactor file
        f_in->Close();

        if (apply_oscillation || is_further_reactors){
          // work out total oscillated integral of constraints
          int normalisation_unosc = reactor_unosc_1d_pdf[j][i]->Integral();
          int normalisation_reactor = reactor_osc_1d_pdf[j][i]->Integral();
          Double_t osc_loss = normalisation_reactor/(double)normalisation_unosc;
	
          osc_losses[i].push_back(osc_loss);
	  
          Double_t constraint_osc_mean = reactor_distance_to_evs(distances_fit[rot_i][j][i],P_over_D_coeff,thermal_powers[i])*osc_loss*n_years;
	
          reactor_osc_1d_pdf[j][i]->Normalise(); 
          //reactor_osc_1d_pdf[j][i]->Scale(constraint_osc_mean);

          nd_reactor_Int_fit += constraint_osc_mean;
          //data_set_pdf.Add(*reactor_osc_1d_pdf[j][i]);
          //printf("  added reactor %d/%d: %s, osc_survival: %.3f, norm_constraint: %.3f data_int:%.0f\n", i+1, n_reactors, reactor_names[i].c_str(), osc_loss,constraint_osc_mean, data_set_pdf.Integral());
          printf(" detector %d/%d: %s, added reactor %d/%d: %s, %i/%i -> osc_survival: %.3f, norm_constraint: %.3f \n", j+1, n_detectors, detector_names[j].c_str(), i+1, n_reactors, reactor_names[i].c_str(), normalisation_reactor, normalisation_unosc, osc_loss,constraint_osc_mean);
          //printf(" detector %d/%d: %s, added reactor %d/%d: %s, osc_survival: %.3f, norm_constraint: %.3f \n", j+1, n_detectors, detector_names[j].c_str(), i+1, n_reactors, reactor_names[i].c_str(), osc_loss,constraint_osc_mean);
	
        }else if (geos_included){
          Double_t normalisation_unosc = reactor_unosc_1d_pdf[j][i]->Integral();
          Double_t normalisation_reactor = reactor_osc_1d_pdf[j][i]->Integral();
          Double_t osc_loss = normalisation_reactor/normalisation_unosc;
	
          double scale = 0.;
	
          // geo custom flux = coefficient times 1 years predicted geos
          reactor_osc_1d_pdf[j][i]->Scale(0.);//geo_uth_custom_flux/(double)flux_data);
	
          //data_set_pdf.Add(*reactor_osc_1d_pdf[j][i]);
          //printf("  added geo %d/%d: %s, osc_survival: %.3f, data_int:%.0f\n", i+1, n_reactors, reactor_names[i].c_str(), osc_loss, data_set_pdf.Integral());
	
        }else if (alphans_included){
          Double_t normalisation_unosc = reactor_unosc_1d_pdf[j][i]->Integral();
          Double_t normalisation_reactor = reactor_osc_1d_pdf[j][i]->Integral();
          Double_t osc_loss = normalisation_reactor/normalisation_unosc;
	
          double scale = 0.;
	  
          reactor_osc_1d_pdf[j][i]->Normalise(); // currently custom flux number = total events in livetime
          reactor_osc_1d_pdf[j][i]->Scale(0.);//bkg_alphan_custom_flux);///(double)flux_data);
	
          //data_set_pdf.Add(*reactor_osc_1d_pdf[j][i]);
          //printf("  added alpha-n %d/%d: %s, osc_survival: %.3f, data_int:%.0f\n", i+1, n_reactors, reactor_names[i].c_str(), osc_loss, data_set_pdf.Integral());
        }
      
        /////////////// end of detector loop for given reactor ^ /////////////
        if (j == n_detectors-1){
          reactor_osc_nd_pdf[i] = new BinnedED(reactor_names[i], axes);
          reactor_osc_nd_pdf[i]->SetObservables(data_rep);
	
          std::vector<Double_t> e1_e2_vec(n_detectors, 0.);
          for (size_t l = 0; l < n_bins; l++){
            for (size_t k = 0; k < n_bins; k++){
              std::vector<Double_t> e1_bin_centre_vals(1, 0.);
              axes_1d.GetBinCentres(l, e1_bin_centre_vals);
              e1_e2_vec[0] = e1_bin_centre_vals[0];
              std::vector<Double_t> e2_bin_centre_vals(1, 0.);
              axes_1d.GetBinCentres(k, e2_bin_centre_vals);
              e1_e2_vec[1] = e2_bin_centre_vals[0];
	      
              Double_t bin_content = reactor_osc_1d_pdf[0][i]->GetBinContent(l);
              for (size_t nd = 1; nd < n_detectors; nd++)
                bin_content *= reactor_osc_1d_pdf[nd][i]->GetBinContent(k);
	      
              reactor_osc_nd_pdf[i]->Fill(e1_e2_vec, bin_content);
            }
          }
          reactor_osc_nd_pdf[i]->Normalise();

          // Setting optimisation limits
          sprintf(name, "%s_norm", reactor_names[i].c_str());
          Double_t min = 0.;
          Double_t max = 1.5*nd_reactor_Int_fit;
          if (reactor_names[i] == "Yongbyon")
            max = 2.5*nd_reactor_Int_fit;
          if (min < 0) min = 0;
          minima[name] = min;
          maxima[name] = max;
          Double_t random = random_generator->Uniform(0.8,1.2);
          std::cout<<name<<" -> "<<nd_reactor_Int_fit<<std::endl;
          initial_val[name] = nd_reactor_Int_fit*random;
          initial_err[name] = 0.5*nd_reactor_Int_fit;

          lh_function.AddPdf(*reactor_osc_nd_pdf[i]);
          if (known_reactor_constraints[i] > 0.){
            Double_t sigma =  nd_reactor_Int_fit*known_reactor_constraints[i];
            lh_function.SetConstraint(name, nd_reactor_Int_fit, sigma);
            printf("  added CONSTRAINED reactor %d/%d: %s (min:%.3f max:%.3f) constraint: %.3f , sigma: %.3f data_int:%.0f\n", i+1, n_reactors, reactor_names[i].c_str(), min, max, nd_reactor_Int_fit, sigma, data_set_pdf_integral);
          }
          else
            printf("  added reactor %d/%d: %s (min:%.3f max:%.3f) constraint: %.3f , sigma: %.3f data_int:%.0f\n", i+1, n_reactors, reactor_names[i].c_str(), min, max, nd_reactor_Int_fit, 0., data_set_pdf_integral);

        }
      }
    }

    // fit
    printf("Built LH function, fitting...\n");
    Minuit min;
    min.SetMethod("Migrad");
    min.SetMaxCalls(100000);
    min.SetTolerance(0.01);
    min.SetMinima(minima);
    min.SetMaxima(maxima);
    min.SetInitialValues(initial_val);
    min.SetInitialErrors(initial_err);
  
    FitResult fit_result = min.Optimise(&lh_function);
    fit_result.SetPrintPrecision(9);
    ParameterDict best_fit = fit_result.GetBestFit();
    fit_result.Print();

    //if (rot_i == 0)
    if (fit_validity == 1)
      fit_validity *= fit_result.GetValid();
    else
      fit_validity = fit_result.GetValid();
    

    /////////////
    // checking only for zero rotation!!!
    if (n_fit == 1 && fit_validity && rot_i == 0){
      for (ULong64_t i = 0; i < n_reactors; i++){
        for (ULong64_t j = 0; j < n_detectors; j++){
          if (reactor_names[i] == "Yongbyon"){
            //std::cout<<" FOUND YB i:"<<i<<" j:"<<j<<std::endl;
            if (j == 0){
              reactor_osc_data_1d_hist[j][i]->SetFillColor(8);
              reactor_osc_data_1d_hist[j][i]->SetLineColor(1);
              hs_e1->Add(reactor_osc_data_1d_hist[j][i]);
            }else if (j == 1){
              reactor_osc_data_1d_hist[j][i]->SetFillColor(3);
              reactor_osc_data_1d_hist[j][i]->SetLineColor(1);
              hs_e2->Add(reactor_osc_data_1d_hist[j][i]);
            }
          }
        }
      }
    
      for (ULong64_t i = 0; i < n_reactors; i++){
        for (ULong64_t j = 0; j < n_detectors; j++){
          if (reactor_names[i] != "Yongbyon"){
            //std::cout<<" ADDING OTHER i:"<<i<<" j:"<<j<<std::endl;
            if (j == 0){
              reactor_osc_data_1d_hist[j][i]->SetFillColor(9);
              reactor_osc_data_1d_hist[j][i]->SetLineColor(1);
              hs_e1->Add(reactor_osc_data_1d_hist[j][i]);
            }else if (j == 1){ 
              reactor_osc_data_1d_hist[j][i]->SetFillColor(4);
              reactor_osc_data_1d_hist[j][i]->SetLineColor(1);
              hs_e2->Add(reactor_osc_data_1d_hist[j][i]);
            }
          }
        }
      }

    }

    Double_t lh_val = 99999; // positive non-sensical value to return if fit is not valid
    if (fit_validity == true)
      lh_val = (-1.)*lh_function.Evaluate();

    //Double_t total_power_resolution = 0.;
    //std::vector<Double_t> reactor_power_resolutions(n_reactors, 0.);
    /////////////////////////
    if (rotationZ_0value_index == rot_i){
      for (ULong64_t i = 0; i < n_reactors; i++){
        sprintf(name, "%s_norm", reactor_names[i].c_str());
        std::vector<Double_t> reactor_distances;
        std::vector<Double_t> reactor_osc_losses;
        for (ULong64_t j = 0; j < n_detectors; j++){
          reactor_distances.push_back(distances[j][i]);
          reactor_osc_losses.push_back(osc_losses[i][j]); // losses: { [reactor A losses to detectors] ,  [reactor B losses to detectors] ..}
        }
        Double_t reactor_therm_power = nd_evs_integral_to_power(best_fit.at(name), n_detectors, reactor_distances, reactor_osc_losses, P_over_D_coeff);
        Double_t reactor_therm_power_uncert = (reactor_therm_power - thermal_powers[i])/thermal_powers[i];
        //total_power_resolution += reactor_therm_power_uncert;
        Double_t true_data_nd_int = 0.;
        for (ULong64_t j = 0; j < n_detectors; j++)
          true_data_nd_int += reactor_osc_data_1d_pdf[j][i]->Integral();
        std::cout<<name<<" best fit nd Integral: "<<best_fit.at(name)<<"  true nd Int: "<<true_data_nd_int<<" -> Fit therm power: "<<reactor_therm_power<<"  true: therm power: "<<thermal_powers[i]<<" resolution: "<<reactor_therm_power_uncert<<std::endl;
        reactor_power_resolutions[i] = reactor_therm_power_uncert;

        if (reactor_names[i] == "Yongbyon"){
          if (best_fit.at(name) < 5.)
            fit_validity = 0;
        }
      }
    }
    //total_power_resolution = total_power_resolution/(double)n_reactors;
  
    /*for (ULong64_t j = 0; j < n_reactors; j++){
      sprintf(name, "%s_norm", reactor_names[j].c_str());
      reactor_osc_pdf[j]->Normalise();
      reactor_osc_pdf[j]->Scale(best_fit.at(name));
      reactor_osc_pdf_fitosc_sum.Add(*reactor_osc_pdf[j]);
      }
  
      if (geos_included)
      fit_geo_uth_norm = best_fit.at("uraniumthorium_norm");
      else
      fit_geo_uth_norm = 0.;
  
      // write plots to file (only 'good' plots - those with the best fit values)
      if (param_d21>=param_d21_plot_min && param_d21<=param_d21_plot_max && param_s12>=param_s12_plot_min && param_s12<=param_s12_plot_max){
      file_out->cd();
      TCanvas* c_data_fit = new TCanvas("c_data_fit");
      c_data_fit->cd();

      // and their sum
      TH1D reactor_hist_fitosc_sum = DistTools::ToTH1D(reactor_osc_pdf_fitosc_sum);
      sprintf(name, "reactor_hist_fitosc_sum_d21%.9f_s12%.7f_s13%.7f", param_d21, param_s12, param_s13);
      reactor_hist_fitosc_sum.SetName(name);
      reactor_hist_fitosc_sum.GetXaxis()->SetTitle("Energy (MeV)");
      reactor_hist_fitosc_sum.GetYaxis()->SetTitle("Counts");
      reactor_hist_fitosc_sum.SetLineColor(kRed);
      reactor_hist_fitosc_sum.Write();
      reactor_hist_fitosc_sum.Draw("same");
    
      // data set
      TH1D data_set_hist = DistTools::ToTH1D(data_set_pdf);
      // data_hist.Sumw2();
      sprintf(name, "data_set_hist");
      data_set_hist.SetName(name);
      data_set_hist.GetYaxis()->SetTitle("Counts");
      data_set_hist.GetXaxis()->SetTitle("Energy (MeV)");
      data_set_hist.SetLineColor(kBlue);
      data_set_hist.Write();
      data_set_hist.Draw("same");
    
      file_out->cd();
      c_data_fit->Write();
    
      // reactor pdfs
      TH1D reactor_osc_hist;
      for (ULong64_t j = 0; j < n_reactors; j++){
      reactor_osc_hist = DistTools::ToTH1D(*reactor_osc_pdf[j]);
      // data_hist.Sumw2();
      sprintf(name, "reactor_osc_pdf_%s_d21%.9f_s12%.7f_s13%.7f", reactor_names[j].c_str(), param_d21, param_s12, param_s13);
      reactor_osc_hist.SetName(name);
      reactor_osc_hist.GetYaxis()->SetTitle("Counts");
      reactor_osc_hist.GetXaxis()->SetTitle("Energy (MeV)");
      reactor_osc_hist.Write();
      }
    
      // pdfs of spectra
      TH1D reactor_unosc_hist;
      for (ULong64_t j = 0; j < n_reactors; j++){
      reactor_unosc_hist = DistTools::ToTH1D(*reactor_unosc_pdf[j]);
      // reactor_unosc_hist.Sumw2();
      sprintf(name, "reactor_unosc_pdf_%s_d21%.9f_s12%.7f_s13%.7f", reactor_names[j].c_str(), param_d21, param_s12, param_s13);
      reactor_unosc_hist.SetName(name);
      reactor_unosc_hist.GetYaxis()->SetTitle("Counts");
      reactor_unosc_hist.GetXaxis()->SetTitle("Energy (MeV)");
      reactor_unosc_hist.Write();
      }
      }
    */
    printf("fit valid: %d, lh_value:%.9f\n", fit_validity, lh_val);
    printf("End fit--------------------------------------\n");
    //lh_value = lh_val;
    lh_values_vec_rot[rot_i] = lh_val;
  }
  //return lh_val;
  return reactor_power_resolutions;//total_power_resolution;
  //return 0;
}

int main(int argc, char *argv[]) {

  if (argc < 12){
      std::cout<<"Error: 11(12) argument expected."<<std::endl;
      return 1; // return>0 indicates error code
  }
  else{
    const std::string reactor_monitoring_folder = argv[1];
    const Double_t reactor_separation_d = atof(argv[2]);

    const Double_t number_of_ktonnes = atof(argv[3]);
    const size_t number_of_fits = atoi(argv[4]);

    const std::string detector_names_csv_file = argv[5];
    //"/data/snoplus/blakei/antinu/reactor_monitoring/csv_files/detector_names.csv";

    const std::string reactor_names_csv_file = argv[6];
      //"/data/snoplus/blakei/antinu/reactor_monitoring/csv_files/reactor_names.csv";

    const bool poisson_fluc = atoi(argv[7]);

    const std::string out_file_name = argv[8];
    
    Double_t rotationZ_min = atof(argv[9]);
    Double_t rotationZ_max = atof(argv[10]);
    Double_t rotationZ_n = atof(argv[11]);

    size_t job_id = 0;
    if(argc == 13)
      job_id = atoi(argv[12]);
    
    const std::string &info_file = "";
    const std::string &spectrum_phwr_unosc_filepath = "";
    const std::string &spectrum_pwr_unosc_filepath = "/data/snoplus2/blakei/Monitoring/pruned_ntuples/charlie_scaling/reactors/latest_osc/pdfs/pwr_pdf_flux30000_combined_cleanround299.ntuple.root";
    const std::string &spectrum_geo_uraniumthorium_unosc_filepath = "";
    const std::string &spectrum_bkg_alphan_unosc_filepath = "";
    const std::string &constraints_info_file = "";
    const double s12 = 0.306;
    const double d21 = 7.51e-5;
    const double s13 = 0.0215;
    const size_t x_bin = 0;
    const size_t y_bin = 0;
    const double flux_data = 1000.;
    const double mc_scale_factor = 1.;
    const double e_min = 0.9;
    const double e_max = 7.7;
    const size_t n_bins = 68;//51; //40;
    const double param_d21_plot_min = 9e-5;
    const double param_d21_plot_max = 9e-5;
    const double param_s12_plot_min = 0.29;
    const double param_s12_plot_max = 0.3;
    const bool constrained_data = true;

    std::stringstream _dstream;
    _dstream<<"_d"<<reactor_separation_d<<"km";
    std::stringstream _rotstream;
    _rotstream<<"_rot"<<rotationZ_min<<"_"<<rotationZ_max<<"_"<<rotationZ_n;
    //if (rotationZ != 0.)
    //_dstream<<"_rot"<<rotationZ<<"deg";
    std::stringstream job_id_stream;
    job_id_stream<<"_r"<<job_id;
    const std::string &data_path = reactor_monitoring_folder+"sim_2d_data"+_dstream.str()+job_id_stream.str()+".root";
    const std::string &out_filename_csv = reactor_monitoring_folder+"output"+_dstream.str()+job_id_stream.str()+".csv";
    const std::string &out_filename_plots = reactor_monitoring_folder+"plot_2d_output"+_dstream.str()+job_id_stream.str()+".root";
    const std::string &out_filename_res_plots = reactor_monitoring_folder+out_file_name+"_"+_dstream.str()+job_id_stream.str()+".root";

    
    //double P_over_D_coeff = 219.*number_of_ktonnes;
    double cut_effic = 0.68; //<- ~clean300
    double P_over_D_coeff = 408.282*cut_effic*number_of_ktonnes;

    double n_years = 1.;
    
    std::ifstream in_csv_reactor;
    in_csv_reactor.open(reactor_names_csv_file.c_str());
    
    std::ifstream in_csv_detector;
    in_csv_detector.open(detector_names_csv_file.c_str());


        ///////////////////////////////////
    //// detector positions, names ////
    std::vector<std::string> detector_names;
    std::vector<Double_t> detector_latitudes;
    std::vector<Double_t> detector_longitudes;
    std::vector<Double_t> detector_altitudes;
    
    std::string detector_name,detector_latitude,detector_longitude,detector_altitude;
    ULong64_t line_no = 0;

    while(in_csv_detector.good()){
      std::getline(in_csv_detector,detector_name,',');
      std::getline(in_csv_detector,detector_latitude,',');
      std::getline(in_csv_detector,detector_longitude,',');
      std::getline(in_csv_detector,detector_altitude,'\n');

      if (line_no>0){ //skip csv header
	    if (strcmp(detector_name.c_str(),"")!=0) {
          detector_names.push_back(detector_name);
          detector_latitudes.push_back(atof(detector_latitude.c_str()));
          detector_longitudes.push_back(atof(detector_longitude.c_str()));
          detector_altitudes.push_back(atof(detector_altitude.c_str()));
	    }
      }
      line_no++;
    }
    in_csv_detector.close();
    
    const size_t n_detectors = detector_names.size();

    
    std::vector<TVector3> detector_positions;
    for (size_t i = 0; i < detector_names.size(); i++)
      detector_positions.push_back(LLAtoECEF(detector_latitudes[i], detector_longitudes[i], detector_altitudes[i]));
    
    TVector3 mean_detector_pos = TVector3(0., 0., 0.);
    for (size_t i = 0; i < detector_positions.size(); i++)
      mean_detector_pos += detector_positions[i];
    mean_detector_pos = mean_detector_pos*(1./(double)detector_positions.size());
    
    //This assumes flat earth lols
    for (size_t i = 0; i < detector_positions.size(); i++){
      TVector3 mid_point_to_detector_dir = (detector_positions[i] - mean_detector_pos)*(1./((detector_positions[i] - mean_detector_pos).Mag()));
      detector_positions[i] = (mid_point_to_detector_dir*reactor_separation_d) + mean_detector_pos;
    }

    if (detector_names.size() != detector_positions.size()){
      std::cout<<"detector name and pos not equal!"<<std::endl;
      std::cout<<"detector names: "<<detector_names.size()<<" positions: "<<detector_positions.size()<<std::endl;
      exit(0);
    }

    //const double avg_total_therm_power = 5000.;
    /////////////////////////////////////////
    //// reactor names, positions types  ////    
    std::vector<std::string> reactor_names;
    std::vector<Double_t> reactor_latitude;
    std::vector<Double_t> reactor_longitude;
    std::vector<Double_t> reactor_altitude;
    std::vector<Double_t> thermal_powers;
    std::vector<Double_t> known_reactor_constraints;
   
    
    std::string reactor_name,latitude,longitude,altitude,power,constrain;
    line_no = 0;
    
    while(in_csv_reactor.good()){
      std::getline(in_csv_reactor,reactor_name,',');
      std::getline(in_csv_reactor,latitude,',');
      std::getline(in_csv_reactor,longitude,',');
      std::getline(in_csv_reactor,altitude,',');
      std::getline(in_csv_reactor,power,',');
      std::getline(in_csv_reactor,constrain,'\n');

      if (line_no>0){ //skip csv header
	    if (strcmp(reactor_name.c_str(),"")!=0) {
          reactor_names.push_back(reactor_name);
          reactor_latitude.push_back(atof(latitude.c_str()));
          reactor_longitude.push_back(atof(longitude.c_str()));
          reactor_altitude.push_back(atof(altitude.c_str()));
          thermal_powers.push_back(atof(power.c_str()));
          known_reactor_constraints.push_back(atof(constrain.c_str()));
	    }
      }
      line_no++;
    }
    in_csv_reactor.close();
    
    const size_t n_reactors = reactor_names.size();

    std::vector<TVector3> reactor_positions;
    std::vector<std::vector<TVector3> >reactor_positions_fit;
    for (size_t i = 0; i < reactor_latitude.size(); i++){
      TVector3 reactor_pos = LLAtoECEF(reactor_latitude[i], reactor_longitude[i], reactor_altitude[i]);
      reactor_positions.push_back(reactor_pos);
      /*if (rotationZ != 0){
       if (reactor_names[i] == "Yongbyon"){
          std::cout<<"Rot: Yongbyon: "<<reactor_pos.X()<<","<<reactor_pos.Y()<<","<<reactor_pos.Z();
          rotationZ = rotationZ*3.14159265359/180.;  //deg to rad
          reactor_pos.Rotate(rotationZ,mean_detector_pos);
          std::cout<<" -> "<<reactor_pos.X()<<","<<reactor_pos.Y()<<","<<reactor_pos.Z()<<std::endl;
        }
      }
      reactor_positions_fit.push_back(reactor_pos);*/
    }
    double del_rotation = (rotationZ_max - rotationZ_min)/(double)rotationZ_n;
    int rotationZ_0value_index;
    for (int rot_i = 0; rot_i < rotationZ_n ; rot_i++){
      double rotationZ = rotationZ_min + rot_i*del_rotation;
      if (rotationZ == 0.){
        rotationZ_0value_index = rot_i;
        std::cout<<"\n\n 0 rotation index == "<<rotationZ_0value_index<<"\n\n"<<std::endl;
      }
      std::cout<<"Rotation: "<<rotationZ;
      std::vector<TVector3> reactor_positions_fit_rotation;
      for (size_t i = 0; i < reactor_latitude.size(); i++){
        TVector3 reactor_pos = LLAtoECEF(reactor_latitude[i], reactor_longitude[i], reactor_altitude[i]);
        if (reactor_names[i] == "Yongbyon"){
          std::cout<<"Rot: Yongbyon: "<<reactor_pos.X()<<","<<reactor_pos.Y()<<","<<reactor_pos.Z();
          rotationZ = rotationZ*3.14159265359/180.;  //deg to rad
          reactor_pos.Rotate(rotationZ,mean_detector_pos);
          std::cout<<" -> "<<reactor_pos.X()<<","<<reactor_pos.Y()<<","<<reactor_pos.Z()<<std::endl;
        }
        reactor_positions_fit_rotation.push_back(reactor_pos);
      }
      reactor_positions_fit.push_back(reactor_positions_fit_rotation);
    }
    
    std::vector<std::string> reactor_types;
    for (size_t i = 0; i < reactor_names.size(); i++)
      reactor_types.push_back("PWR");

    if (reactor_names.size() != reactor_positions.size() || reactor_names.size() != thermal_powers.size()){
      std::cout<<"reactor name, powers and pos not equal!"<<std::endl;
      exit(0);
    }
          
    std::vector<std::vector<Double_t> > distances; // { [distances to 1] ,  [distances to 2] }
    for (size_t i = 0; i < n_detectors; i++){
      std::cout<<"\n reactor distances to detector "<<detector_names[i]<<": "<<std::endl;
      std::vector<Double_t> distances_to_detector;
      distances.push_back(distances_to_detector);
      //std::vector<Double_t> distances_fit_to_detector;
      //distances_fit.push_back(distances_fit_to_detector);
      for (size_t j = 0; j < n_reactors; j++){
        distances[i].push_back((detector_positions[i] - reactor_positions[j]).Mag());
        //distances_fit[i].push_back((detector_positions[i] - reactor_positions_fit[j]).Mag());
        std::cout<<detector_names[i]<<"<->"<<reactor_names[j]<<" = "<<distances[i][j]<<"km"<<std::endl;
      }
    }

    std::vector<std::vector<std::vector<Double_t> > >distances_fit; // { [distances to 1] ,  [distances to 2] }
    for (int rot_i = 0; rot_i < rotationZ_n ; rot_i++){
      std::vector<std::vector<Double_t> > distances_to_detector_rot;
      distances_fit.push_back(distances_to_detector_rot);
      for (size_t i = 0; i < n_detectors; i++){
        std::vector<Double_t> distances_fit_to_detector;
        distances_fit[rot_i].push_back(distances_fit_to_detector);
        for (size_t j = 0; j < n_reactors; j++){
          distances_fit[rot_i][i].push_back((detector_positions[i] - reactor_positions[j]).Mag());
        }
      }
    }

    /*if (rotationZ != 0.){
      for (size_t i = 0; i < n_detectors; i++){
        std::cout<<"\n ROTATED reactor distances to detector "<<detector_names[i]<<": "<<std::endl;
        for (size_t j = 0; j < n_reactors; j++){
          std::cout<<detector_names[i]<<"<->"<<reactor_names[j]<<" = "<<distances_fit[i][j]<<"km"<<std::endl;
        }
      }
      }*/
    std::cout<<"Centre of Detector positions: ("<<mean_detector_pos[0]<<", "<<mean_detector_pos[1]<<", "<<mean_detector_pos[2]<<")"<<std::endl;
    
    /*const size_t n_inter_ds = 10;
    const Double_t min_inter_d = 0.;
    const Double_t max_inter_d = 100.;
    Double_t del_inter_d = max_inter_d/(double)n_inter_ds;
    */
    //std::vector<std::vector<Double_t> > reactor_power_res_vs_d_vec;  //n_reactors
    //std::vector<Double_t> d_vals;
    //for (size_t i = 0; i < n_reactors; i++){
    /*std::stringstream reactor_vsd_name_strm;
      reactor_vsd_name_strm<<"reactor_"<<reactor_names[i]<<"_total_power_resolution_vs_d";
      TH1D h_power_res_vs_d(reactor_vsd_name_strm.str().c_str(), reactor_vsd_name_strm.str().c_str(), n_inter_ds, min_inter_d, max_inter_d);
    */
    //std::vector<Double_t> power_res_vals;
    //reactor_power_res_vs_d_vec.push_back(power_res_vals);
    //}
    
    std::vector<std::vector<TH1D> > reactor_power_uncert_vs_d_vec; //{ [dist1: {reac1, reac2, ..}]  [dist2: {reac1, reac2, ..}]}
    size_t n_ds = 0;
    //for (Double_t reactor_separation_d = min_inter_d; reactor_separation_d <= max_inter_d; reactor_separation_d += del_inter_d){
    //d_vals.push_back(reactor_separation_d);
    n_ds += 1;
    std::vector<TH1D> reactor_power_uncert_at_d_vec;  //n_reactors
    
    for (size_t i = 0; i < n_reactors; i++){
      std::stringstream reactor_power_uncert_strm;
      reactor_power_uncert_strm<<reactor_names[i]<<"_power_uncert_d_"<<reactor_separation_d;
      TH1D h_reactor_total_power_uncert_d(reactor_power_uncert_strm.str().c_str(), reactor_power_uncert_strm.str().c_str(), 100, -1., 1.);
	
      reactor_power_uncert_at_d_vec.push_back(h_reactor_total_power_uncert_d);
    }
      
    printf("Begin--------------------------------------\n");
  
    //readInfoFile(info_file, reactor_names, distances, reactor_types, n_cores, powers, power_errs);
      
    //reactor_names.insert(reactor_names.begin(), reactor_names_arr, reactor_names_arr + n_reactors);
    //distances.insert(distances.begin(), distances_arr, distances_arr + n_reactors);
            
    // read constraint info for each reactor in the info file (one at time to ensure they match correctly)
    /*for (size_t i=0; i<(size_t)reactor_names.size(); i++){
      double fit_mean, fit_mean_err, fit_sigma, fit_sigma_err;
      readConstraintsInfoFile(constraints_info_file, reactor_names[i].c_str(), fit_mean, fit_mean_err, fit_sigma, fit_sigma_err);
      constraint_means.push_back(fit_mean);
      constraint_mean_errs.push_back(fit_mean_err);
      constraint_sigmas.push_back(fit_sigma);
      constraint_sigma_errs.push_back(fit_sigma_err);
      }
    */
    for (size_t i = 0; i < n_detectors; i++){
      for (size_t j = 0; j < n_reactors; j++){
        printf("detector_name:%s, reactor_name:%s, fit_mean: %.3f\n", detector_names[i].c_str(),reactor_names[j].c_str(), reactor_distance_to_evs(distances[i][j],P_over_D_coeff,thermal_powers[j]));
      }
    }
      
    // initialise data
    /*if (constrained_data == false)
      LHFit_initialise(data_set_pdf, data_path, flux_data, e_min, e_max, n_bins);
      else
      LHFit_load_fake_data(data_set_pdf, data_path, flux_data, e_min, e_max, n_bins);
    */
    ////save objects to file
    printf("Save objects to file...\n");
    TFile *file_out = new TFile(out_filename_plots.c_str(), "RECREATE");
    bool fit_validity = 0;
    ULong64_t fit_try_max = 30;
    ULong64_t print_plots = 0;
      
    printf("running: d_21:%.9f(%.9f-%.9f) s_12:%.7f(%.7f-%.7f)\n", d21, param_d21_plot_min, param_d21_plot_max, s12, param_s12_plot_min, param_s12_plot_max);
    if (d21 >= param_d21_plot_min && d21 <= param_d21_plot_max && s12 >= param_s12_plot_min && s12<=param_s12_plot_max){
      printf("writing plots to: %s\n", out_filename_plots.c_str());
      print_plots++;
    }
      
    printf("print plots: %d\n\n ", print_plots);
    //printf("Fit number: %llu of %llu\n", i+1, n_parameter_sets);
      
    //double lh_value = 999999;
    //double total_power_resolution = 9999;
    std::vector<Double_t> reactor_power_uncert(n_reactors, 0.);
      
    double fit_geo_uth_norm = -999999;
      
    fit_validity = 0;

    THStack * hs_e1 = new THStack("hs_e1", "Simulated Data: YongByon (green) reactor + rest at detector 1");
    THStack * hs_e2 = new THStack("hs_e2", "Simulated Data: YongByon (green) reactor + rest at detector 2");

    TH1D h_data_set_detector_1;
    TH1D h_data_set_detector_2;
    TH1D h_data_set_detector_1_poisson_fluc;
    TH1D h_data_set_detector_2_poisson_fluc;

    TH2D h2_data;
    TH2D * h2_best_fit = new TH2D("h2_best_fit","h2_best_fit",n_bins, e_min, e_max, n_bins, e_min, e_max);

    std::vector<std::vector<double> >lh_values_vec;

    for (ULong64_t n_fit=1; n_fit<=number_of_fits; n_fit++) {
      std::cout<<"\n\nd: "<<reactor_separation_d<<" n_ktonnes:  "<<number_of_ktonnes<<"  fit: "<<n_fit<<"/"<<number_of_fits<<std::endl;

      for (ULong64_t fit_try=1; fit_try<=fit_try_max; fit_try++) {
        std::vector<double> lh_values_vec_rot; //double lh_value;
        for (int rot_i = 0; rot_i < rotationZ_n ; rot_i++)
          lh_values_vec_rot.push_back(999999);
        reactor_power_uncert = LHFit_fit(spectrum_phwr_unosc_filepath,
                                         spectrum_pwr_unosc_filepath,
                                         spectrum_geo_uraniumthorium_unosc_filepath,
                                         spectrum_bkg_alphan_unosc_filepath,
                                         detector_names,
                                         reactor_names, reactor_types,
                                         distances,
                                         distances_fit,
                                         known_reactor_constraints,/*constraint_sigmas,*/
                                         file_out,
                                         d21, s12, s13,
                                         fit_validity, e_min, e_max, n_bins,
                                         flux_data, mc_scale_factor,
                                         fit_geo_uth_norm,
                                         param_d21_plot_min, param_d21_plot_max,
                                         param_s12_plot_min, param_s12_plot_max,
                                         P_over_D_coeff, thermal_powers,
                                         data_path, hs_e1, hs_e2, n_fit,
                                         h2_data, h2_best_fit, poisson_fluc,
                                         h_data_set_detector_1, h_data_set_detector_2,
                                         h_data_set_detector_1_poisson_fluc,
                                         h_data_set_detector_2_poisson_fluc,
                                         n_years, lh_values_vec_rot, rotationZ_n,
                                         rotationZ_0value_index);
	
        if (fit_validity==0)
          printf("Fit invalid... retrying (attempt no: %llu)\n", fit_try);
        else{
          printf("Fit valid. (attempt no: %llu)\n", fit_try);
          lh_values_vec.push_back(lh_values_vec_rot);//lh_value);
          for (ULong64_t i = 0; i < n_reactors; i++){
            reactor_power_uncert_at_d_vec[i].Fill(reactor_power_uncert[i]);
          }
          fit_try = fit_try_max+1;
        }
      }
    }
      
    // close output file
    file_out->Close();
      
    if (print_plots==0) { //if no plots passed the plot cuts, then delete the then empty output file.
      usleep(10000); // wait for the file to finish writing
      if (remove(out_filename_plots.c_str()) != 0)
	printf("Error: deletetion of temporary file not successful...\n");
    }
    //Write fit coefficients to txt file
    //printf("writing to: %s\n", out_filename_csv.c_str());
    //ofstream outcsvfile;
    //outcsvfile.open(out_filename_csv.c_str(), std::ios_base::app);
    //outcsvfile<<s12<<","<<d21<<","<<x_bin<<","<<y_bin<<","<<lh_value<<","<<fit_validity<<","<<fit_geo_uth_norm<<"\n";
    //outcsvfile<<s12<<","<<d21<<","<<x_bin<<","<<y_bin<<","<<"\n";
    //outcsvfile.close();
    
    reactor_power_uncert_vs_d_vec.push_back(reactor_power_uncert_at_d_vec);
    //}
  
    //std::vector<std::vector<TH1D> > reactor_power_uncert_vs_d_vec; //{ [dist1: {reac1, reac2, ..}]  [dist2: {reac1, reac2, ..}]}
    /*
      std::vector<TH1D> reactor_power_res_vs_d_vec;  //n_reactors
      for (size_t i = 0; i < n_reactors; i++){
      std::stringstream reactor_vsd_name_strm;
      reactor_vsd_name_strm<<"reactor_"<<reactor_names[i]<<"_total_power_uncert_d_"<<reactor_separation_d;
      TH1D h_power_res_vs_d(reactor_vsd_name_strm.c_str(), reactor_vsd_name_strm.c_str(), n_inter_ds, 0, max_inter_d);
      }
    */

    TFile *file_res_plots_out = new TFile(out_filename_res_plots.c_str(), "RECREATE");
    if (job_id == 0){
      //THStack * hs_e1 = new THStack("hs_e1", "Simulated Data: YongByon (green) reactor + rest at detector 1");
      /*TH1D h_data_set_YB_pdf_e1 = DistTools::ToTH1D(data_set_YB_pdf_e1);
      h_data_set_YB_pdf_e1.SetName("h_data_set_YB_pdf_e1");
      h_data_set_YB_pdf_e1.SetFillColor(8);
      hs_e1->Add(&h_data_set_YB_pdf_e1);
      TH1D h_data_set_sans_YB_pdf_e1 = DistTools::ToTH1D(data_set_sans_YB_pdf_e1);
      h_data_set_sans_YB_pdf_e1.SetName("h_data_set_sans_YB_pdf_e1");
      h_data_set_sans_YB_pdf_e1.SetFillColor(9);
      hs_e1->Add(&h_data_set_sans_YB_pdf_e1);

      THStack * hs_e2 = new THStack("hs_e2", "Simulated Data: YongByon (green) reactor + rest at detector 2");
      TH1D h_data_set_YB_pdf_e2 = DistTools::ToTH1D(data_set_YB_pdf_e2);
      h_data_set_YB_pdf_e2.SetName("h_data_set_YB_pdf_e2");
      h_data_set_YB_pdf_e2.SetFillColor(3);
      hs_e2->Add(&h_data_set_YB_pdf_e2);
      TH1D h_data_set_sans_YB_pdf_e2 = DistTools::ToTH1D(data_set_sans_YB_pdf_e2);
      h_data_set_sans_YB_pdf_e2.SetName("h_data_set_sans_YB_pdf_e2");
      h_data_set_sans_YB_pdf_e2.SetFillColor(4);
      hs_e2->Add(&h_data_set_sans_YB_pdf_e2);
      */
      TCanvas *c_poisson_v_no_poisson_data_detector_1 = new TCanvas("c_poisson_v_no_poisson_data_detector_1");
      h_data_set_detector_1_poisson_fluc.Draw("same");
      h_data_set_detector_1.Draw("same");
      c_poisson_v_no_poisson_data_detector_1->Write();
      TCanvas *c_poisson_v_no_poisson_data_detector_2 = new TCanvas("c_poisson_v_no_poisson_data_detector_2");
      h_data_set_detector_2_poisson_fluc.Draw("same");
      h_data_set_detector_2.Draw("same");
      c_poisson_v_no_poisson_data_detector_2->Write();
      
      hs_e1->Write();
      hs_e2->Write();
      h2_data.Write();
    }

    for (size_t j = 0; j < n_ds; j++){
      for (size_t i = 0; i < n_reactors; i++){
        reactor_power_uncert_vs_d_vec[j][i].Write();
        //Double_t hist_mean = reactor_power_uncert_vs_d_vec[j][i].GetMean();
        //Double_t hist_sigma = reactor_power_uncert_vs_d_vec[j][i].GetStdDev();
        //reactor_power_res_vs_d_vec[i].SetBinContent(j+1, hist_sigma);
        //reactor_power_res_vs_d_vec[i].push_back(hist_sigma);
      }
    }
    /*for (size_t i = 0; i < n_reactors; i++){
      std::stringstream reactor_vsd_name_strm;
      reactor_vsd_name_strm<<"reactor_"<<reactor_names[i]<<"_total_power_resolution_vs_d";
      std::vector<Double_t> power_res_vals = reactor_power_res_vs_d_vec[i];
      TGraph g_power_res_vs_d(d_vals.size(), &d_vals[0], &power_res_vals[0]);
      g_power_res_vs_d.SetName(reactor_vsd_name_strm.str().c_str());
      g_power_res_vs_d.Write();
      }
    */

    TCanvas* c_detector_reactor_pos = new TCanvas("c_detector_reactor_pos");

    // plotting detector and reactor positions
    TMultiGraph *mg = new TMultiGraph();
    double detector_posx[n_detectors];
    double detector_posy[n_detectors];
    double reactor_posx[n_reactors];
    double reactor_posy[n_reactors];
    for (size_t i = 0; i < n_detectors; i++){
      //detector_posx[i] = (mean_detector_pos - detector_positions[i]).X();
      //detector_posy[i] = (mean_detector_pos - detector_positions[i]).Y();
      //detector_posx[i] = (180 - detector_latitudes[i]);
      //detector_posy[i] = (180 - detector_longitudes[i]);
      detector_posx[i] = (detector_longitudes[i]);
      detector_posy[i] = (detector_latitudes[i]);
    }
    for (size_t i = 0; i < n_reactors; i++){
      //reactor_posx[i] = (mean_detector_pos - reactor_positions[i]).X();
      //reactor_posy[i] = (mean_detector_pos - reactor_positions[i]).Y();
      //reactor_posx[i] = (180 - reactor_latitude[i]);
      //reactor_posy[i] = (180 - reactor_longitude[i]);
      reactor_posx[i] = (reactor_longitude[i]);
      reactor_posy[i] = (reactor_latitude[i]);
    }
    
    TGraph *gr_detpos = new TGraph(n_detectors, detector_posx, detector_posy);
    //TGraph *gr_detpos = new TGraph(n_detectors, &detector_latitudes[0], &detector_longitudes[0]);
    gr_detpos->SetMarkerColor(2);
    gr_detpos->SetMarkerStyle(7);
    for (size_t i = 0; i < n_detectors; i++){
      std::stringstream name_str;
      name_str<<detector_names[i];//<<"_"<<detector_posx[i]<<"_"<<detector_posy[i];
      TLatex *latex = new TLatex(gr_detpos->GetX()[i], gr_detpos->GetY()[i],name_str.str().c_str());
      latex->SetTextSize(0.03);
      gr_detpos->GetListOfFunctions()->Add(latex); gr_detpos->Draw("alp");
    }
    TGraph *gr_reacpos = new TGraph(n_reactors, reactor_posx, reactor_posy);
    //TGraph *gr_reacpos = new TGraph(n_reactors, &reactor_latitude[0], &reactor_longitude[0]);
    gr_reacpos->SetMarkerColor(4);
    gr_reacpos->SetMarkerStyle(7);
    for (size_t i = 0; i < n_reactors; i++){
      std::stringstream name_str;
      name_str<<reactor_names[i];//<<"_"<<reactor_posx[i]<<"_"<<reactor_posy[i];
      TLatex *latex = new TLatex(gr_reacpos->GetX()[i], gr_reacpos->GetY()[i],name_str.str().c_str());
      latex->SetTextSize(0.03);
      gr_reacpos->GetListOfFunctions()->Add(latex); gr_reacpos->Draw("alp");
    }

    mg->Add(gr_detpos);
    mg->Add(gr_reacpos);

    mg->Draw("AP");
    c_detector_reactor_pos->Write();

    file_res_plots_out->Close();

    ///////////////
    // record lh_values in ntuple
    std::string out_LHtree_filename_res_plots = reactor_monitoring_folder+"LHTree_"+out_file_name+"_"+_dstream.str()+_rotstream.str()+job_id_stream.str()+".root";

    TFile *file_LH_plots_out = new TFile(out_LHtree_filename_res_plots.c_str(), "RECREATE");
    
    ///////////////////
    // filling every rot_i^th entry to be fit value for that rotation
    TNtuple *LH_vals_ntp = new TNtuple("nt", "LH_values", "lh_value:rot_value");
    double lh_val_ntp;
    double rot_value;
    //for (int rot_i = 0; rot_i < rotationZ_n ; rot_i++){
    //lh_val_ntp.push_back(99999);
    for (size_t i = 0; i < lh_values_vec.size(); i++){
      for (int rot_i = 0; rot_i < rotationZ_n ; rot_i++){
        double rotationZ = rotationZ_min + rot_i*del_rotation;
        lh_val_ntp = lh_values_vec[i][rot_i];
        rot_value = rotationZ;
        LH_vals_ntp->Fill(lh_val_ntp, rot_value);
      }
    }
    //////////////

    LH_vals_ntp->Write();

    file_LH_plots_out->Close();
    
    printf("End--------------------------------------\n");
      
    return 0; // completed successfully
  }
}
