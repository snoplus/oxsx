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
#include <TGraph.h>
#include <THStack.h>

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

Double_t evs_to_power(double int_evs, double distance, double osc_loss, const double P_over_D_coeff){
  //calc in steps
  double avg_total_therm_power = int_evs/(double)P_over_D_coeff;
  avg_total_therm_power *= (distance*distance/osc_loss);
  avg_total_therm_power = avg_total_therm_power;
  return avg_total_therm_power;
}


std::vector<Double_t>
LHFit_fit(const std::string &spectrum_phwr_unosc_filepath,
	  const std::string &spectrum_pwr_unosc_filepath,
	  const std::string &spectrum_geo_uraniumthorium_unosc_filepath,
	  const std::string &spectrum_bkg_alphan_unosc_filepath,
	  const std::vector<std::string> &detector_names,
	  const std::vector<std::string> &reactor_names, std::vector<std::string> &reactor_types,
	  std::vector< Double_t > &distances,
	  std::vector<Double_t> &known_reactor_constraints, std::vector<Double_t> &constraint_sigmas,
	  TFile *file_out,
	  Double_t param_d21, Double_t param_s12, Double_t param_s13,
	  bool &fit_validity,
	  const double e_min, const double e_max, const size_t n_bins,
	  const double flux_data, const double mc_scale_factor,
	  double &fit_geo_uth_norm,
	  const double param_d21_plot_min, const double param_d21_plot_max, const double param_s12_plot_min, 
	  const double param_s12_plot_max, const double P_over_D_coeff, const std::vector<double> thermal_powers,
          const std::string &data_path, THStack *hs, size_t n_fit, std::string reactor_to_monitor, double n_years){


  //////////////////////////////
  ///// Constrained Data  //////
  //////////////////////////////
  
  printf("Beginning Fake data prep--------------------------------------\n");
  printf("Fake Data:: del_21:%.9f, sin2_12:%.7f, sin2_13:%.7f\n", param_d21, param_s12, param_s13);

  char name[1000];
  TRandom3 *random_generator = new TRandom3();
  random_generator->SetSeed(0);
  
  const ULong64_t n_reactors = reactor_names.size();

  // set up binning
  ObsSet data_rep(0);
  AxisCollection axes;
  axes.AddAxis(BinAxis("ev_prompt_fit", e_min, e_max, n_bins));

  BinnedED data_set_pdf("data_set_pdf", axes);
  data_set_pdf.SetObservables(data_rep);
  
  BinnedED **reactor_unosc_data_pdf = new BinnedED*[n_reactors];
  BinnedED **reactor_osc_data_pdf = new BinnedED*[n_reactors];

  TH1D **reactor_osc_data_hist = new TH1D*[n_reactors];

  bool geos_included = false;
  bool alphans_included = false;
  
  for (ULong64_t i = 0; i < n_reactors; i++){
    sprintf(name, "%s_unosc", reactor_names[i].c_str());
    reactor_unosc_data_pdf[i] = new BinnedED(name, axes);
    reactor_unosc_data_pdf[i]->SetObservables(0);
    reactor_osc_data_pdf[i] = new BinnedED(reactor_names[i], axes);
    reactor_osc_data_pdf[i]->SetObservables(0);

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
    if (apply_oscillation)
      ntOscillate_pruned(reactor_unosc_ntp, reactor_osc_ntp, param_d21, param_s12, param_s13, distances[i]);
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
    for(size_t j = 0; j < reactor_unosc_ntp->GetEntries(); j++){
      reactor_unosc_ntp->GetEntry(j);
      reactor_unosc_data_pdf[i]->Fill(ev_unosc_energy_p1);
    }

    // fill oscillated pdf
    Float_t ev_osc_energy_p1;
    reactor_osc_ntp->SetBranchAddress("ev_fit_energy_p1", &ev_osc_energy_p1);
    for(size_t j = 0; j < reactor_osc_ntp->GetEntries(); j++){
      reactor_osc_ntp->GetEntry(j);
      reactor_osc_data_pdf[i]->Fill(ev_osc_energy_p1);
    }

    // close unoscillated reactor file
    f_in->Close();

	if (apply_oscillation || is_further_reactors){
	  // work out total oscillated integral of constraints
	  int normalisation_unosc = reactor_unosc_data_pdf[i]->Integral();
	  int normalisation_reactor = reactor_osc_data_pdf[i]->Integral();
	  double osc_loss = normalisation_reactor/(double)normalisation_unosc;

	  Double_t constraint_osc_mean = reactor_distance_to_evs(distances[i],P_over_D_coeff,thermal_powers[i])*osc_loss*mc_scale_factor*n_years;
	  
	  reactor_osc_data_pdf[i]->Normalise(); 


	  /*TH1D *hist_for_smoothing = new TH1D(DistTools::ToTH1D(*reactor_osc_data_pdf[i]));
        hist_for_smoothing->Smooth();
        std::vector<Double_t> smooth_bin_contents;
        for (int bin = 1; bin < hist_for_smoothing->GetNbinsX()+1; bin++){
	    Double_t content = hist_for_smoothing->GetBinContent(bin);
	    smooth_bin_contents.push_back(content);
        }
        reactor_osc_data_pdf[i]->SetBinContents(smooth_bin_contents);
	  */
	  
	  reactor_osc_data_pdf[i]->Scale(constraint_osc_mean);

	  data_set_pdf.Add(*reactor_osc_data_pdf[i]);

	  if (n_fit == 1)
	    reactor_osc_data_hist[i] = new TH1D(DistTools::ToTH1D(*reactor_osc_data_pdf[i]));

      printf("added reactor %d/%d: %s, %u/%u -> osc_survival: %.3f, norm_constraint: %.3f data_int:%.0f \n", i+1, n_reactors, reactor_names[i].c_str(), normalisation_reactor, normalisation_unosc, osc_loss,constraint_osc_mean, data_set_pdf.Integral());
	  
	}else if (geos_included){
	  Double_t normalisation_unosc = reactor_unosc_data_pdf[i]->Integral();
	  Double_t normalisation_reactor = reactor_osc_data_pdf[i]->Integral();
	  Double_t osc_loss = normalisation_reactor/normalisation_unosc;

	  double scale = 0.;
	  
	  // geo custom flux = coefficient times 1 years predicted geos
	  reactor_osc_data_pdf[i]->Scale(0.);//geo_uth_custom_flux/(double)flux_data);

	  data_set_pdf.Add(*reactor_osc_data_pdf[i]);

	  printf("  added geo %d/%d: %s, osc_survival: %.3f, data_int:%.0f\n", i+1, n_reactors, reactor_names[i].c_str(), osc_loss, data_set_pdf.Integral());
	
	}else if (alphans_included){
	  Double_t normalisation_unosc = reactor_unosc_data_pdf[i]->Integral();
	  Double_t normalisation_reactor = reactor_osc_data_pdf[i]->Integral();
	  Double_t osc_loss = normalisation_reactor/normalisation_unosc;

	  double scale = 0.;
	  
	  reactor_osc_data_pdf[i]->Normalise(); // currently custom flux number = total events in livetime
	  reactor_osc_data_pdf[i]->Scale(0.);//bkg_alphan_custom_flux);///(double)flux_data);

	  data_set_pdf.Add(*reactor_osc_data_pdf[i]);

	  printf("  added alpha-n %d/%d: %s, osc_survival: %.3f, data_int:%.0f\n", i+1, n_reactors, reactor_names[i].c_str(), osc_loss, data_set_pdf.Integral());
	
	}
	
  }
  
  // Poisson_fluctuations
  for(size_t i = 0; i < data_set_pdf.GetNBins(); i++){
    Double_t old_content = data_set_pdf.GetBinContent(i);
    Double_t new_content = random_generator->Poisson(old_content);
    data_set_pdf.SetBinContent(i, new_content);
  }
  /*TFile *file_data_out = new TFile(data_path.c_str(), "RECREATE");
  TH1D DataDist = DistTools::ToTH1D(data_set_pdf);
  DataDist.Write();
  file_data_out->Close();
  */

  ////////////////////////////////////
  ////////////  Fitting  /////////////
  ////////////////////////////////////
  
  
  printf("Begin fit--------------------------------------\n");
  printf("LHFit_fit:: del_21:%.9f, sin2_12:%.7f, sin2_13:%.7f\n", param_d21, param_s12, param_s13);

  ObsSet fit_rep(0);
  
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

  TH1D *reactor_hist = new TH1D[n_reactors];

  BinnedED reactor_osc_pdf_fitosc_sum("reactor_osc_pdf_fitosc_sum",axes);
  reactor_osc_pdf_fitosc_sum.SetObservables(data_rep);

  BinnedED **reactor_unosc_pdf = new BinnedED*[n_reactors];
  BinnedED **reactor_osc_pdf = new BinnedED*[n_reactors];

  Double_t data_set_pdf_integral = data_set_pdf.Integral();

  geos_included = false;

  std::vector<Double_t> osc_losses;
  for (ULong64_t i = 0; i < n_reactors; i++){
    // for each reactor, load spectrum pdf for reactor type
    sprintf(name, "%s_unosc", reactor_names[i].c_str());
    reactor_unosc_pdf[i] = new BinnedED(name, axes);
    reactor_unosc_pdf[i]->SetObservables(0);
    reactor_osc_pdf[i] = new BinnedED(reactor_names[i], axes);
    reactor_osc_pdf[i]->SetObservables(0);

    bool apply_oscillation = false;
    bool is_further_reactors = false;
    if ((reactor_types[i]=="PWR")||(reactor_types[i]=="BWR")){
      sprintf(name, "%s", spectrum_pwr_unosc_filepath.c_str());
      apply_oscillation = true;
      std::cout<<"reactor: "<<reactor_names[i];
    }else if (reactor_types[i]=="PHWR"){
      sprintf(name, "%s", spectrum_phwr_unosc_filepath.c_str());
      std::cout<<"reactor: "<<reactor_names[i];
      apply_oscillation = true;
    }else if (reactor_types[i]=="further_reactors"){
      sprintf(name, "%s", spectrum_pwr_unosc_filepath.c_str());
      is_further_reactors = true;
      std::cout<<"further reactor";
    }else if (reactor_names[i]=="uraniumthorium"){
      sprintf(name, "%s", spectrum_geo_uraniumthorium_unosc_filepath.c_str());
      geos_included = true;	
    }else if (reactor_names[i]=="alphan"){
      sprintf(name, "%s", spectrum_bkg_alphan_unosc_filepath.c_str());
    }else{
      printf("Throw: Reactor doesn't match any loaded type...\n");
      exit(0); // throw std::exception(); //continue;
    }

    TFile *f_in = new TFile(name);
    file_out->cd(); // switch to output file (for ntuple to use)
    TTree *reactor_unosc_ntp = (TTree*)f_in->Get("nt");
    TNtuple *reactor_osc_ntp = new TNtuple("nt", "Oscillated Prompt Energy", "ev_fit_energy_p1");

    // oscillate tree
    if (apply_oscillation){
      ntOscillate_pruned(reactor_unosc_ntp, reactor_osc_ntp, param_d21, param_s12, param_s13, distances[i]);
      std::cout<<" oscillate ";
    }else if (is_further_reactors){
      ntOscillate_pruned_geo(reactor_unosc_ntp, reactor_osc_ntp, param_s12);
      std::cout<<" oscillate_geo";
    }else{
      ntNoOscillate_pruned(reactor_unosc_ntp, reactor_osc_ntp);
    }

    // reset branch addresses after oscillating in function (otherwise crash before setting again below..)
    reactor_unosc_ntp->SetBranchStatus("*", 0);
    reactor_unosc_ntp->SetBranchStatus("ev_fit_energy_p1", 1); // (re-enable all branches in use)

    // fill unoscillated pdf
    Double_t ev_unosc_energy_p1;
    reactor_unosc_ntp->SetBranchAddress("ev_fit_energy_p1", &ev_unosc_energy_p1);
    for(size_t j = 0; j < reactor_unosc_ntp->GetEntries(); j++){
      reactor_unosc_ntp->GetEntry(j);
      reactor_unosc_pdf[i]->Fill(ev_unosc_energy_p1);
    }

    // fill oscillated pdf
    Float_t ev_osc_energy_p1;
    reactor_osc_ntp->SetBranchAddress("ev_fit_energy_p1", &ev_osc_energy_p1);
    for(size_t j = 0; j < reactor_osc_ntp->GetEntries(); j++){
      reactor_osc_ntp->GetEntry(j);
      reactor_osc_pdf[i]->Fill(ev_osc_energy_p1);
    }

    // close unoscillated reactor file
    f_in->Close();

    if (apply_oscillation || is_further_reactors){
      // work out total oscillated integral of constraints
      int normalisation_unosc = reactor_unosc_pdf[i]->Integral();
      int normalisation_reactor = reactor_osc_pdf[i]->Integral();
      double osc_loss = normalisation_reactor/(double)normalisation_unosc;

      osc_losses.push_back(osc_loss);
      
      Double_t constraint_osc_mean = reactor_distance_to_evs(distances[i],P_over_D_coeff,thermal_powers[i])*osc_loss*mc_scale_factor;
      //Double_t constraint_osc_sigma = constraint_osc_mean*0.1;
      reactor_osc_pdf[i]->Normalise(); //remove number of events from mc
      reactor_unosc_pdf[i]->Scale(1./flux_data); // osc pdf gets fitted, the unosc doesn't, scale it simply for plotting..
      
      // Setting optimisation limits
      sprintf(name, "%s_norm", reactor_names[i].c_str());
      Double_t min = 0.;//constraint_osc_mean-2.*constraint_osc_sigma; // let min and max float within 2 sigma
      Double_t max = 3.*constraint_osc_mean;//constraint_osc_mean+2.*constraint_osc_sigma;
      if (min < 0) min = 0;
        minima[name] = min;
      maxima[name] = max;
      //printf("  added reactor %d/%d: %s, osc_survival: %.3f, norm_constraint: %.3f (min:%.3f max:%.3f) err: %.3f data_int:%.0f\n", i+1, n_reactors, reactor_names[i].c_str(), osc_loss,constraint_osc_mean, min, max, constraint_osc_sigma, data_set_pdf_integral);
      Double_t random = random_generator->Uniform(0.5,1.5);
      initial_val[name] = constraint_osc_mean*random;
      initial_err[name] = 0.5*constraint_osc_mean; //constraint_osc_sigma;

      lh_function.AddPdf(*reactor_osc_pdf[i]);
      //lh_function.SetConstraint(name, constraint_osc_mean, constraint_osc_sigma);
      if (known_reactor_constraints[i] > 0.){
	  Double_t sigma =  constraint_osc_mean*known_reactor_constraints[i];
	  lh_function.SetConstraint(name, constraint_osc_mean,sigma);
	  printf("  added CONSTRAINED reactor %d/%d: %s, %u/%u -> osc_survival: %.3f, norm_constraint: %.3f (min:%.3f max:%.3f) data_int:%.0f\n", i+1, n_reactors, reactor_names[i].c_str(),  normalisation_reactor, normalisation_unosc, osc_loss,constraint_osc_mean, min, max, data_set_pdf_integral);
	}
	else
	  printf("  added reactor %d/%d: %s, osc_survival: %.3f, norm_constraint: %.3f (min:%.3f max:%.3f) data_int:%.0f\n", i+1, n_reactors, reactor_names[i].c_str(), osc_loss,constraint_osc_mean, min, max, data_set_pdf_integral);
    }
    else{
      // Setting optimisation limits
      sprintf(name, "%s_norm", reactor_names[i].c_str());
      Double_t min = 0; // let min and max float within 2 sigma
      Double_t max = 500;
      if (min < 0) min = 0;
      minima[name] = min;
      maxima[name] = max;
      printf("  added reactor %d/%d: %s, norm: %.3f (min:%.3f max:%.3f) sigma: %.3f data_int:%.0f\n", i+1, n_reactors, reactor_names[i].c_str(), 0 , min, max, 0, data_set_pdf_integral);
      Double_t random = random_generator->Uniform(0.5,1.5);
      initial_val[name] = min + (max-min)*random;
      initial_err[name] = min + (max-min)*random;

      lh_function.AddPdf(*reactor_osc_pdf[i]);
      //lh_function.SetConstraint(name, constraint_osc_mean, constraint_osc_sigma);
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
  fit_validity = fit_result.GetValid();


  if (n_fit == 1 && fit_validity){
    for (ULong64_t i = 0; i < n_reactors; i++){
      if (reactor_names[i] == reactor_to_monitor){
	reactor_osc_data_hist[i]->SetFillColor(8);
	reactor_osc_data_hist[i]->SetLineColor(1);
	hs->Add(reactor_osc_data_hist[i]);
      }
    }
    
    for (ULong64_t i = 0; i < n_reactors; i++){
      if (reactor_names[i] != reactor_to_monitor){
	reactor_osc_data_hist[i]->SetFillColor(9);
	reactor_osc_data_hist[i]->SetLineColor(1);
	hs->Add(reactor_osc_data_hist[i]);
      }
    }
  }

  Double_t lh_val = 99999; // positive non-sensical value to return if fit is not valid
  if (fit_validity == true)
    lh_val = (-1.)*lh_function.Evaluate(); //lh_function.Evaluate();

  std::vector<Double_t> reactor_power_uncerts(n_reactors, 0.);
  for (ULong64_t i = 0; i < n_reactors; i++){
    sprintf(name, "%s_norm", reactor_names[i].c_str());
    Double_t reactor_therm_power = evs_to_power(best_fit.at(name), distances[i], osc_losses[i], P_over_D_coeff);
    Double_t reactor_therm_power_uncert = (reactor_therm_power - thermal_powers[i])/thermal_powers[i];
    //total_power_resolution += reactor_therm_power_uncert;
    std::cout<<name<<" best fit nd Integral: "<<best_fit.at(name)<<" -> Fit therm power: "<<reactor_therm_power<<"  true: therm power: "<<thermal_powers[i]<<" resolution: "<<reactor_therm_power_uncert<<std::endl;
    reactor_power_uncerts[i] = reactor_therm_power_uncert;
  }
  

  /*
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
  //return lh_val;
  return reactor_power_uncerts;
}

int main(int argc, char *argv[]) {

  if (argc < 4){
      std::cout<<"Error: 3(4) argument expected."<<std::endl;
      return 1; // return>0 indicates error code
  }
  else{
    const std::string reactor_monitoring_folder = argv[1];
    const Double_t number_of_ktonnes = atof(argv[2]);
    const size_t number_of_fits = atoi(argv[3]);

    size_t job_id = 0;
    if(argc == 5)
      job_id = atoi(argv[4]);

    const std::string &info_file = "";
    const std::string &spectrum_phwr_unosc_filepath = "";
    //"/data/snoplus/blakei/antinu/reactor_monitoring/pdfs/pwr_pdf_flux1_day365_passcombined5000_cleanround4.ntuple.root";//"/data/snoplusmc/antinu/ntuples/snop_rat6176_2017/reactors/scintFitter/flux1/pdfs/pwr_pdf_flux1_day365_passcombined3000_cleanround4.ntuple.root";
    const std::string &spectrum_pwr_unosc_filepath = "/data/snoplus2/blakei/Monitoring/pruned_ntuples/charlie_scaling/reactors/latest_osc/pdfs/pwr_pdf_flux30000_combined_cleanround299.ntuple.root";
    //"/data/snoplus/blakei/antinu/reactor_monitoring/pdfs/pwr_pdf_flux1_day365_passcombined5000_cleanround4.ntuple.root";//"/data/snoplusmc/antinu/ntuples/snop_rat6176_2017/reactors/scintFitter/flux1/pdfs/pwr_pdf_flux1_day365_passcombined3000_cleanround4.ntuple.root";
    const std::string &spectrum_geo_uraniumthorium_unosc_filepath = "";
    //"/data/snoplusmc/antinu/ntuples/snop_rat6169_2017/reactors/scintFitter/flux1/pdfs/pwr_pdf_flux1_day365_passcombined1000_cleanround4.ntuple.root";
    const std::string &spectrum_bkg_alphan_unosc_filepath = "";
    //"/data/snoplusmc/antinu/ntuples/snop_rat6169_2017/reactors/scintFitter/flux1/pdfs/pwr_pdf_flux1_day365_passcombined1000_cleanround4.ntuple.root";

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
    const size_t n_bins = 68;
    const double param_d21_plot_min = 9e-5;
    const double param_d21_plot_max = 9e-5;
    const double param_s12_plot_min = 0.29;
    const double param_s12_plot_max = 0.3;
    const bool constrained_data = true;

    //std::stringstream _dstream;
    //_dstream<<"_d"<<reactor_separation_d<<"km";
    std::stringstream job_id_stream;
    job_id_stream<<"_r"<<job_id;
    const std::string &data_path = reactor_monitoring_folder+"sim_1d_data_1_"+job_id_stream.str()+".root";
    const std::string &out_filename_csv = reactor_monitoring_folder+"output"+job_id_stream.str()+".csv";
    const std::string &out_filename_plots = reactor_monitoring_folder+"plot_2d_output"+job_id_stream.str()+".root";
    const std::string &out_filename_res_plots = reactor_monitoring_folder+"plot_e_uncert"+job_id_stream.str()+".root";
    

    printf("Begin--------------------------------------\n");

    //double P_over_D_coeff = 219.*number_of_ktonnes;
    double cut_effic = 0.68; //<- ~clean300
    double P_over_D_coeff = 408.282*cut_effic*number_of_ktonnes;

    double n_years = 1.;
    
    const std::string reactor_to_monitor = "Yongbyon";
    std::vector<std::string> reactor_names;
    reactor_names.push_back("Hanul"); // Reactor A
    reactor_names.push_back("Hanbit"); // Reactor B
    reactor_names.push_back("Wolsong"); // Reactor C
    reactor_names.push_back("Kori_ShinKori"); // Reactor D

    reactor_names.push_back("Yongbyon"); // Reactor E
    //reactor_names.push_back("Taechon"); // Reactor F

    reactor_names.push_back("Hongyanhe"); // Reactor F
    reactor_names.push_back("Tianwan"); // Reactor F
    const size_t n_reactors = reactor_names.size();

    std::vector<Double_t> thermal_powers;
    thermal_powers.push_back(6*2800.);
    thermal_powers.push_back(6*2800.);
    thermal_powers.push_back(4*2060.);
    thermal_powers.push_back(1730+1880+2900+2900+(4.*2825));

    thermal_powers.push_back(150.);
    //thermal_powers.push_back(600.);

    thermal_powers.push_back(3.*2905);
    thermal_powers.push_back(2.*3000);
    
    std::vector<Double_t> known_reactor_constraints;
    known_reactor_constraints.push_back(0.06);
    known_reactor_constraints.push_back(0.06);
    known_reactor_constraints.push_back(0.06);
    known_reactor_constraints.push_back(0.06);

    known_reactor_constraints.push_back(0.99);
    //known_reactor_constraints.push_back(-1.);

    known_reactor_constraints.push_back(0.06);
    known_reactor_constraints.push_back(0.06);
    
    std::vector<TVector3> reactor_locations;
    reactor_locations.push_back(TVector3(36.979203, 129.385462, 0.) );
    reactor_locations.push_back(TVector3(35.4105, 126.4175, 0.) );
    reactor_locations.push_back(TVector3(35.713450, 129.475759, 0.) );
    reactor_locations.push_back(TVector3(35.320423, 129.290547, 0.) );

    reactor_locations.push_back(TVector3(39.797325, 125.755063, 0.) );
    //reactor_locations.push_back(TVector3(39.932491, 125.485561, 0.) );

    reactor_locations.push_back(TVector3(39.795132, 121.481271, 0.) );
    reactor_locations.push_back(TVector3(34.686477, 119.459362, 0.) );
    

    // Iran //
    /*const std::string reactor_to_monitor = "Bushehr";

    std::vector<std::string> reactor_names;
    reactor_names.push_back("Bushehr"); // Reactor A
    reactor_names.push_back("Barakah"); // Reactor A
    
    const size_t n_reactors = reactor_names.size();

    std::vector<Double_t> thermal_powers;
    thermal_powers.push_back(3000.);
    thermal_powers.push_back(4*4000.);
    
    std::vector<Double_t> known_reactor_constraints;
    known_reactor_constraints.push_back(-1.);
    known_reactor_constraints.push_back(-1.);

    std::vector<TVector3> reactor_locations;
    reactor_locations.push_back(TVector3(28.829717, 50.886004, 0.) );
    reactor_locations.push_back(TVector3(23.968670, 52.231329, 0.) );
    */
    

    std::vector<TVector3> reactor_positions;
    for (size_t i = 0; i < reactor_locations.size(); i++){
      reactor_positions.push_back(LLAtoECEF(reactor_locations[i][0], reactor_locations[i][1], reactor_locations[i][2]));
    }
    std::vector<std::string> reactor_types;
    for (size_t i = 0; i < reactor_names.size(); i++)
      reactor_types.push_back("PWR");

    if (reactor_names.size() != reactor_positions.size() || reactor_names.size() != thermal_powers.size()){
      std::cout<<"reactor name, powers and pos not equal!"<<std::endl;
      exit(0);
    }
    

    ///////////////////////////////////
    //// detector positions, names ////
    std::vector<std::string> detector_names;
    detector_names.push_back("1"); // Detector 1
    const size_t n_detectors = detector_names.size();

    std::vector<TVector3> detector_locations;
    //detector_locations.push_back(TVector3(38.236731, 127.090735, 0.) ); //SK NK border
    detector_locations.push_back(TVector3(40.522792, 124.905926, 0.) );  //NK China border

    //detector_locations.push_back(TVector3(29.795038, 48.255271, 0.) ); //Kuwait - Iran
    
    std::vector<TVector3> detector_positions;
    for (size_t i = 0; i < detector_locations.size(); i++)
      detector_positions.push_back(LLAtoECEF(detector_locations[i][0], detector_locations[i][1], detector_locations[i][2]));

    
    if (detector_names.size() != detector_positions.size()){
      std::cout<<"detector name and pos not equal!"<<std::endl;
      exit(0);
    }
      
    std::vector<Double_t> distances; // { [distances to 1] ,  [distances to 2] }
    std::cout<<"\n reactor distances to detector "<<detector_names[0]<<": "<<std::endl;
    for (size_t j = 0; j < n_reactors; j++){
      distances.push_back((detector_positions[0] - reactor_positions[j]).Mag());
      std::cout<<detector_names[0]<<"<->"<<reactor_names[j]<<" = "<<distances[j]<<"km"<<std::endl;
    }


    std::vector<std::vector<TH1D> > reactor_power_uncert_vs_d_vec; //{ [dist1: {reac1, reac2, ..}]  [dist2: {reac1, reac2, ..}]}
    size_t n_ds = 0;
    //for (Double_t reactor_separation_d = min_inter_d; reactor_separation_d <= max_inter_d; reactor_separation_d += del_inter_d){
    //d_vals.push_back(reactor_separation_d);
    n_ds += 1;
    std::vector<TH1D> reactor_power_uncert_at_d_vec;  //n_reactors
    
    for (size_t i = 0; i < n_reactors; i++){
      std::stringstream reactor_power_uncert_strm;
      reactor_power_uncert_strm<<reactor_names[i]<<"_power_uncert";//_d_"<<reactor_separation_d;
      TH1D h_reactor_total_power_uncert_d(reactor_power_uncert_strm.str().c_str(), reactor_power_uncert_strm.str().c_str(), 100, -1., 1.);
	
      reactor_power_uncert_at_d_vec.push_back(h_reactor_total_power_uncert_d);
    }    
    
    
    // read in constraint information
    std::vector<Double_t> constraint_means;
    std::vector<Double_t> constraint_mean_errs;
    std::vector<Double_t> constraint_sigmas;
    std::vector<Double_t> constraint_sigma_errs;

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
	printf("detector_name:%s, reactor_name:%s, fit_mean: %.3f\n", detector_names[0].c_str(),reactor_names[j].c_str(), reactor_distance_to_evs(distances[j],P_over_D_coeff,thermal_powers[j]));
      }
    }

    //const ULong64_t n_reactors = reactor_names.size();

    AxisCollection axes;
    axes.AddAxis(BinAxis("ev_prompt_fit", e_min, e_max, n_bins));

    ObsSet data_rep(0);
    BinnedED data_set_pdf("data_set_pdf", axes);
    data_set_pdf.SetObservables(data_rep);
    
    ////save objects to file
    printf("Save objects to file...\n");
    TFile *file_out = new TFile(out_filename_plots.c_str(), "RECREATE");
    bool fit_validity = 0;
    ULong64_t fit_try_max = 20;
    ULong64_t print_plots = 0;

    printf("running: d_21:%.9f(%.9f-%.9f) s_12:%.7f(%.7f-%.7f)\n", d21, param_d21_plot_min, param_d21_plot_max, s12, param_s12_plot_min, param_s12_plot_max);
    if (d21 >= param_d21_plot_min && d21 <= param_d21_plot_max && s12 >= param_s12_plot_min && s12<=param_s12_plot_max){
      printf("writing plots to: %s\n", out_filename_plots.c_str());
      print_plots++;
    }

    printf("print plots: %d\n\n ", print_plots);
    //printf("Fit number: %llu of %llu\n", i+1, n_parameter_sets);

    //double lh_value = 999999;

    double fit_geo_uth_norm = -999999;
    
    fit_validity = 0;

    std::stringstream hs_strm;
    hs_strm<<"Simulated Data: "<<reactor_to_monitor<<" (green) reactor + rest";
      
    THStack * hs = new THStack("hs", hs_strm.str().c_str());

    std::vector<Double_t> reactor_power_uncert(n_reactors, 0.);  
    for (ULong64_t n_fit=1; n_fit<=number_of_fits; n_fit++) {
      std::cout<<"\n\n"<<"  fit: "<<n_fit<<"/"<<number_of_fits<<std::endl;

      for (ULong64_t fit_try=1; fit_try<=fit_try_max; fit_try++) {
        reactor_power_uncert = LHFit_fit(spectrum_phwr_unosc_filepath,
                                         spectrum_pwr_unosc_filepath,
                                         spectrum_geo_uraniumthorium_unosc_filepath,
                                         spectrum_bkg_alphan_unosc_filepath,
                                         detector_names,
                                         reactor_names, reactor_types,
                                         distances,
                                         known_reactor_constraints, constraint_sigmas,
                                         file_out,
                                         d21, s12, s13,
                                         fit_validity, e_min, e_max, n_bins,
                                         flux_data, mc_scale_factor,
                                         fit_geo_uth_norm,
                                         param_d21_plot_min, param_d21_plot_max,
                                         param_s12_plot_min, param_s12_plot_max,
                                         P_over_D_coeff, thermal_powers,
                                         data_path, hs, n_fit, reactor_to_monitor,
                                         n_years);
	
        if (fit_validity==0)
          printf("Fit invalid... retrying (attempt no: %llu)\n", fit_try);
        else{
          printf("Fit valid. (attempt no: %llu)\n", fit_try);
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
    /*printf("writing to: %s\n", out_filename_csv.c_str());
      ofstream outcsvfile;
      outcsvfile.open(out_filename_csv.c_str(), std::ios_base::app);
      outcsvfile<<s12<<","<<d21<<","<<x_bin<<","<<y_bin<<","<<"\n";
      outcsvfile.close();
    */
    reactor_power_uncert_vs_d_vec.push_back(reactor_power_uncert_at_d_vec);
    
    TFile *file_res_plots_out = new TFile(out_filename_res_plots.c_str(), "RECREATE");
    if (job_id == 0){
      hs->Write();
    }
    for (size_t j = 0; j < n_ds; j++){
      for (size_t i = 0; i < n_reactors; i++){
        reactor_power_uncert_vs_d_vec[j][i].Write();
      }
    }
    file_res_plots_out->Close();    
    
    printf("End--------------------------------------\n");
    return 0; // completed successfully
  }
}
