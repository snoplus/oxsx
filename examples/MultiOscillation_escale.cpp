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
#include <GaussianERes.h>
#include <ParameterDict.h>
#include <ContainerTools.hpp>
#include <NuOsc.h>
#include <SurvProb.h>
#include "AntinuUtils.cpp"
#include "../util/oscillate_util_eres.cpp"

Double_t LHFit_fit(BinnedED &data_set_pdf, const std::string &spectrum_phwr_unosc_filepath,
  const std::string &spectrum_pwr_unosc_filepath,
  const std::string &spectrum_geo_uraniumthorium_unosc_filepath,
  const std::string &spectrum_bkg_alphan_unosc_filepath,
  std::vector<std::string> &reactor_names, std::vector<std::string> &reactor_types,
  std::vector<Double_t> &distances,
  std::vector<Double_t> &constraint_means, std::vector<Double_t> &constraint_sigmas,
  TFile *file_out,
  Double_t param_d21, Double_t param_s12, Double_t param_s13,
  bool &fit_validity,
  const double e_min, const double e_max, const size_t n_bins,
  const double flux_data, const double mc_scale_factor,
  double &fit_geo_uth_norm,
  const double param_d21_plot_min, const double param_d21_plot_max, const double param_s12_plot_min, 
  const double param_s12_plot_max, 
  const bool apply_energy_resolution_convolution,
  const Double_t del_energy,
  const Double_t e_resolution_estimate,
  const Double_t e_resolution_estimate_sigma,
  const bool apply_energy_scaling,
  const Double_t e_scaling_estimate,
  const Double_t e_scaling_estimate_sigma){

  printf("Begin fit--------------------------------------\n");
  printf("LHFit_fit:: del_21:%.9f, sin2_12:%.7f, sin2_13:%.7f\n", param_d21, param_s12, param_s13);

  char name[1000];
  TRandom3 *random_generator = new TRandom3();
  const ULong64_t n_pdf = reactor_names.size();

  // set up binning
  ObsSet data_rep(0);
  AxisCollection axes;
  axes.AddAxis(BinAxis("ev_prompt_fit", e_min, e_max, n_bins));

  // create LH function
  BinnedNLLH lh_function;
  lh_function.SetBufferAsOverflow(true);
  int buffL = 1;
  int buffR = 6;
  lh_function.SetBuffer(0, buffL, buffR);
  lh_function.SetDataDist(data_set_pdf); // initialise withe the data set

  // setup max and min ranges
  ParameterDict minima;
  ParameterDict maxima;
  ParameterDict initial_val;
  ParameterDict initial_err;

  TH1D *reactor_hist = new TH1D[n_pdf];

  BinnedED reactor_osc_pdf_fitosc_sum("reactor_osc_pdf_fitosc_sum",axes);
  reactor_osc_pdf_fitosc_sum.SetObservables(data_rep);

  BinnedED **reactor_unosc_pdf = new BinnedED*[n_pdf];
  BinnedED **reactor_osc_pdf = new BinnedED*[n_pdf];

  Double_t constraint_osc_mean_total = 0.;
  Double_t data_set_pdf_integral = data_set_pdf.Integral();

  bool geos_included = false;

  for (ULong64_t i = 0; i < n_pdf; i++){
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

      // load unoscillated reactor file (to oscillate, and to plot)
      //ROOTNtuple reactor_unosc_ntp(spectrum_unosc_filepath.c_str(), "nt"); // this would be made easier if this worked for specific branches!!
      TFile *f_in = new TFile(name);
      file_out->cd(); // switch to output file (for ntuple to use)
      TTree *reactor_unosc_ntp = (TTree*)f_in->Get("nt");
      TNtuple *reactor_osc_ntp = new TNtuple("nt", "Oscillated Prompt Energy", "ev_fit_energy_p1");

      // oscillate tree
      // if apply_energy_resolution_convolution set ev_fit_energy_p1 = mc_neutrino_energy + del_energy (1.022MeV) (e-e+ annhilation photons)
      // pdfs to be smeared by energy convolution
      if (apply_oscillation){
          ntOscillate_pruned(reactor_unosc_ntp, reactor_osc_ntp, param_d21, param_s12, param_s13, distances[i], apply_energy_resolution_convolution, del_energy);
	  std::cout<<" oscillate ";
      }else if (is_further_reactors){
	  ntOscillate_pruned_geo(reactor_unosc_ntp, reactor_osc_ntp, param_s12, apply_energy_resolution_convolution, del_energy);
	  std::cout<<" oscillate_geo";
      }else{
	  ntNoOscillate_pruned(reactor_unosc_ntp, reactor_osc_ntp, apply_energy_resolution_convolution, del_energy);
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
          Double_t normalisation_unosc = reactor_unosc_pdf[i]->Integral();
	  Double_t normalisation_reactor = reactor_osc_pdf[i]->Integral();
	  Double_t osc_loss = normalisation_reactor/normalisation_unosc;

	  Double_t constraint_osc_mean = constraint_means[i]*osc_loss*mc_scale_factor;
	  //Double_t constraint_osc_sigma = (constraint_sigmas[i]/constraint_means[i])*constraint_osc_mean;
	  Double_t constraint_osc_sigma = 0.15*constraint_osc_mean;

	  reactor_osc_pdf[i]->Normalise(); //remove number of events from mc
	  reactor_unosc_pdf[i]->Scale(1./flux_data); // osc pdf gets fitted, the unosc doesn't, scale it simply for plotting..

	  // Setting optimisation limits
	  sprintf(name, "%s_norm", reactor_names[i].c_str());
	  Double_t min = constraint_osc_mean-2.*constraint_osc_sigma; // let min and max float within 2 sigma
	  Double_t max = constraint_osc_mean+2.*constraint_osc_sigma;
	  if (min < 0) min = 0;
	  minima[name] = min;
	  maxima[name] = max;
	  printf("  added reactor %d/%d: %s, osc_survival: %.3f, norm_constraint: %.3f (min:%.3f max:%.3f) err: %.3f data_int:%.0f\n", i+1, n_pdf, reactor_names[i].c_str(), osc_loss,constraint_osc_mean, min, max, constraint_osc_sigma, data_set_pdf_integral);
	  Double_t random = random_generator->Uniform(0.5,1.5);
	  initial_val[name] = constraint_osc_mean*random;
	  initial_err[name] = constraint_osc_sigma;

	  lh_function.AddDist(*reactor_osc_pdf[i]);
	  lh_function.SetConstraint(name, constraint_osc_mean, constraint_osc_sigma);
      }
      else{
          // Setting optimisation limits
	  sprintf(name, "%s_norm", reactor_names[i].c_str());
	  Double_t min = 0.; // let min and max float within 2 sigma
	  Double_t max = 500.;
	  if (min < 0) min = 0;
	  minima[name] = min;
	  maxima[name] = max;
	  printf("  added reactor %d/%d: %s, norm: %.3f (min:%.3f max:%.3f) sigma: %.3f data_int:%.0f\n", i+1, n_pdf, reactor_names[i].c_str(), 0 , min, max, 0, data_set_pdf_integral);
	  Double_t random = random_generator->Uniform(0.5,1.5);
	  initial_val[name] = min + (max-min)*random;
	  initial_err[name] = min + (max-min)*random;

	  lh_function.AddDist(*reactor_osc_pdf[i]);
	  //lh_function.SetConstraint(name, constraint_osc_mean, constraint_osc_sigma); //no geo constraint
      }
  }

  ObsSet  obsSet(0);
  if (apply_energy_resolution_convolution) {
      Convolution* conv = new Convolution("conv");
      GaussianERes* gaus_energy_resolution = new GaussianERes(e_resolution_estimate,"gaus"); 
      gaus_energy_resolution->RenameParameter("eres_0","eres_");
      conv->SetEResolution(gaus_energy_resolution);
      conv->SetAxes(axes);
      conv->SetTransformationObs(obsSet);
      conv->SetDistributionObs(obsSet);
      conv->Construct();
  
      minima["eres_"] = 0;
      maxima["eres_"] = 0.1; 
      initial_val["eres_"] = e_resolution_estimate;
      initial_err["eres_"] = e_resolution_estimate_sigma;
      
      lh_function.AddSystematic(conv);
      //lh_function.SetConstraint("eres_", e_resolution_estimate, e_resolution_estimate_sigma);
  }

  if (apply_energy_scaling) {
      Scale* scale = new Scale("scale");
      scale->RenameParameter("scaleFactor","escale_");
      scale->SetScaleFactor(1.);
      scale->SetAxes(axes);
      scale->SetTransformationObs(obsSet);
      scale->SetDistributionObs(obsSet);
      scale->Construct();

      minima["escale_"] = e_scaling_estimate - e_scaling_estimate_sigma;
      maxima["escale_"] = e_scaling_estimate + e_scaling_estimate_sigma; 
      initial_val["escale_"] = e_scaling_estimate;
      initial_err["escale_"] = e_scaling_estimate_sigma;
      
      lh_function.AddSystematic(scale);
      lh_function.SetConstraint("escale_", e_scaling_estimate, e_scaling_estimate_sigma);
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
  
  Convolution* fitconv = new Convolution("fitconv");
  if (apply_energy_resolution_convolution) {
    GaussianERes* fitgaus = new GaussianERes(best_fit.at("eres_"),"fitgaus"); 
    fitconv->SetEResolution(fitgaus);
    fitconv->SetAxes(axes);
    fitconv->SetTransformationObs(obsSet);
    fitconv->SetDistributionObs(obsSet);
    fitconv->Construct();
  }
  
  Scale* fitscale = new Scale("fitscale");
  if (apply_energy_scaling) {
    fitscale->SetScaleFactor(best_fit.at("escale_"));
    //fitscale->SetScaleFactor(1.);
    fitscale->SetAxes(axes);
    fitscale->SetTransformationObs(obsSet);
    fitscale->SetDistributionObs(obsSet);
    fitscale->Construct();
  }

  for (ULong64_t j = 0; j < n_pdf; j++){
      sprintf(name, "%s_norm", reactor_names[j].c_str());
      reactor_osc_pdf[j]->Normalise();
      BinnedED pdf;
      if (apply_energy_resolution_convolution) {
	pdf = fitconv->operator()(*reactor_osc_pdf[j]);
	pdf.Normalise();
	//pdf.Scale(best_fit.at(name));
      }
      if (apply_energy_scaling) {
	if (!apply_energy_resolution_convolution)
	  pdf = fitscale->operator()(*reactor_osc_pdf[j]);
	else
	  pdf = fitscale->operator()(pdf);
	pdf.Normalise();
	
      }
      if (apply_energy_scaling || apply_energy_resolution_convolution){
	pdf.Scale(best_fit.at(name));
	reactor_osc_pdf_fitosc_sum.Add(pdf);
      }
      if (!apply_energy_scaling && !apply_energy_resolution_convolution){
	reactor_osc_pdf[j]->Scale(best_fit.at(name));
	reactor_osc_pdf_fitosc_sum.Add(*reactor_osc_pdf[j]);
      }
  }

  if (geos_included)
      fit_geo_uth_norm = best_fit.at("uraniumthorium_norm");
  else
      fit_geo_uth_norm = 0.;
  
  Double_t lh_val = 99999; // positive non-sensical value to return if fit is not valid
  if (fit_validity == true)
      lh_val = (-1.)*lh_function.Evaluate(); //lh_function.Evaluate();

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
      for (ULong64_t j = 0; j < n_pdf; j++){
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
      for (ULong64_t j = 0; j < n_pdf; j++){
          reactor_unosc_hist = DistTools::ToTH1D(*reactor_unosc_pdf[j]);
	  // reactor_unosc_hist.Sumw2();
	  sprintf(name, "reactor_unosc_pdf_%s_d21%.9f_s12%.7f_s13%.7f", reactor_names[j].c_str(), param_d21, param_s12, param_s13);
	  reactor_unosc_hist.SetName(name);
	  reactor_unosc_hist.GetYaxis()->SetTitle("Counts");
	  reactor_unosc_hist.GetXaxis()->SetTitle("Energy (MeV)");
	  reactor_unosc_hist.Write();
      }
  }

  printf("fit valid: %d, lh_value:%.9f\n", fit_validity, lh_val);
  printf("End fit--------------------------------------\n");
  return lh_val;
}

int main(int argc, char *argv[]) {

  if (argc != 25){
      std::cout<<"Error: 24 arguments expected."<<std::endl;
      return 1; // return>0 indicates error code
  }
  else{
    const std::string &data_path = argv[1];
    const std::string &info_file = argv[2];
    const std::string &spectrum_phwr_unosc_filepath = argv[3];
    const std::string &spectrum_pwr_unosc_filepath = argv[4];
    const std::string &spectrum_geo_uraniumthorium_unosc_filepath = argv[5];
    const std::string &spectrum_bkg_alphan_unosc_filepath = argv[6];
    const std::string &constraints_info_file = argv[7];
    const double s12 = atof(argv[8]);
    const double d21 = atof(argv[9]);
    const size_t x_bin = atoi(argv[10]);
    const size_t y_bin = atoi(argv[11]);
    const double s13 = atof(argv[12]);
    const double flux_data = atof(argv[13]);
    const double mc_scale_factor = atof(argv[14]);
    const std::string &out_filename_plots = argv[15];
    const std::string &out_filename_csv = argv[16];
    const double e_min = atof(argv[17]);
    const double e_max = atof(argv[18]);
    const size_t n_bins = atoi(argv[19]);
    const double param_d21_plot_min = atof(argv[20]);
    const double param_d21_plot_max = atof(argv[21]);
    const double param_s12_plot_min = atof(argv[22]);
    const double param_s12_plot_max = atof(argv[23]);
    const bool constrained_data = atoi(argv[24]);
    printf("Begin--------------------------------------\n");

    const bool apply_energy_resolution_convolution = false;
    const bool apply_energy_scaling = true;
    double e_resolution_estimate = 0.045;
    double e_resolution_estimate_sigma = 0.005;
    const Double_t del_energy = -0.784;  // energy added to antinu MC KE truth to convert to Eprompt truth
    double e_scaling_estimate = 1.;
    double e_scaling_estimate_sigma = 0.02;

    // read in reactor information
    std::vector<std::string> reactor_names;
    std::vector<Double_t> distances;
    std::vector<std::string> reactor_types;
    std::vector<ULong64_t> n_cores;
    std::vector<Double_t> powers;
    std::vector<Double_t> power_errs;
    readInfoFile(info_file, reactor_names, distances, reactor_types, n_cores, powers, power_errs);

    // read in constraint information
    std::vector<Double_t> constraint_means;
    std::vector<Double_t> constraint_mean_errs;
    std::vector<Double_t> constraint_sigmas;
    std::vector<Double_t> constraint_sigma_errs;

    // read constraint info for each reactor in the info file (one at time to ensure they match correctly)
    for (size_t i=0; i<(size_t)reactor_names.size(); i++){
        double fit_mean, fit_mean_err, fit_sigma, fit_sigma_err;
	readConstraintsInfoFile(constraints_info_file, reactor_names[i].c_str(), fit_mean, fit_mean_err, fit_sigma, fit_sigma_err);
	constraint_means.push_back(fit_mean);
	constraint_mean_errs.push_back(fit_mean_err);
	constraint_sigmas.push_back(fit_sigma);
	constraint_sigma_errs.push_back(fit_sigma_err);
    }

    for (size_t i=0; i<(size_t)reactor_names.size(); i++)
      printf("i:%llu, reactor_name:%s, fit_mean: %.3f, fit_sigma: %.3f\n", i, reactor_names[i].c_str(), constraint_means[i], constraint_sigmas[i]);

    const ULong64_t n_pdf = reactor_names.size();

    AxisCollection axes;
    axes.AddAxis(BinAxis("ev_prompt_fit", e_min, e_max, n_bins));
    BinnedED data_set_pdf("data_set_pdf", axes);

    // initialise data
    if (constrained_data == false)
        LHFit_initialise(data_set_pdf, data_path, flux_data, e_min, e_max, n_bins);
    else
        LHFit_load_fake_data(data_set_pdf, data_path, flux_data, e_min, e_max, n_bins);

    // save objects to file
    printf("Save objects to file...\n");
    TFile *file_out = new TFile(out_filename_plots.c_str(), "RECREATE");
    bool fit_validity = 0;
    ULong64_t fit_try_max = 20;
    ULong64_t print_plots = 0;

    printf("running: d_21:%.9f(%.9f-%.9f) s_12:%.7f(%.7f-%.7f)\n", d21, param_d21_plot_min, param_d21_plot_max, s12, param_s12_plot_min, param_s12_plot_max);
    if (d21 >= param_d21_plot_min && d21 <= param_d21_plot_max && s12 >= param_s12_plot_min && s12<=param_s12_plot_max){
        printf("writing plots to: %s\n", out_filename_plots.c_str());
	print_plots+=1;
    }

    printf("print plots: %d\n\n ", print_plots);
    //printf("Fit number: %llu of %llu\n", i+1, n_parameter_sets);

    double lh_value = 999999;

    double fit_geo_uth_norm = -999999;

    fit_validity = 0;
    for (ULong64_t fit_try=1; fit_try<=fit_try_max; fit_try++) {
      lh_value = LHFit_fit(data_set_pdf, spectrum_phwr_unosc_filepath,
              spectrum_pwr_unosc_filepath,
              spectrum_geo_uraniumthorium_unosc_filepath,
              spectrum_bkg_alphan_unosc_filepath,
              reactor_names, reactor_types,
              distances,
              constraint_means, constraint_sigmas,
              file_out,
              d21, s12, s13,
              fit_validity, e_min, e_max, n_bins,
              flux_data, mc_scale_factor,
	      fit_geo_uth_norm,
              param_d21_plot_min, param_d21_plot_max,
	      param_s12_plot_min, param_s12_plot_max,
	      apply_energy_resolution_convolution,
	      del_energy, e_resolution_estimate,
	      e_resolution_estimate_sigma,
	      apply_energy_scaling,
	      e_scaling_estimate,
	      e_scaling_estimate_sigma);

      if (fit_validity==0)
          printf("Fit invalid... retrying (attempt no: %llu)\n", fit_try);
      else{
          printf("Fit valid. (attempt no: %llu)\n", fit_try);
	  fit_try = fit_try_max+1;
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
    printf("writing to: %s\n", out_filename_csv.c_str());
    ofstream outcsvfile;
    outcsvfile.open(out_filename_csv.c_str(), std::ios_base::app);
    outcsvfile<<s12<<","<<d21<<","<<x_bin<<","<<y_bin<<","<<lh_value<<","<<fit_validity<<","<<fit_geo_uth_norm<<"\n";
    outcsvfile.close();

    printf("End--------------------------------------\n");
    return 0; // completed successfully
  }
}
