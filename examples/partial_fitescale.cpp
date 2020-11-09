// A fit in energy for signal and a background
#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include <Rand.h>
#include <fstream>
#include <iostream>

#include <TFile.h>
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

Double_t LHFit_fit(BinnedED *data_pdf, BinnedED *mc_pdf,
                   bool &fit_validity,
                   const Double_t nhit_scaling_estimate,
                   const Double_t nhit_scaling_estimate_sigma,
                   AxisCollection axes){

  printf("Begin fit--------------------------------------\n");

  TRandom3 *random_generator = new TRandom3();

  // create LH function
  BinnedNLLH lh_function;
  lh_function.SetBufferAsOverflow(true);
  int buffL = 2;
  int buffR = 2;
  lh_function.SetBuffer(0, buffL, buffR);
  lh_function.SetDataDist(*data_pdf); // initialise withe the data set

  // setup max and min ranges
  ParameterDict minima;
  ParameterDict maxima;
  ParameterDict initial_val;
  ParameterDict initial_err;

  //Double_t constraint_osc_mean_total = 0.;
  //Double_t data_set_pdf_integral = data_set_pdf.Integral();

  Double_t minimum = 0;
  Double_t max = 5;;
  if (minimum < 0) minimum = 0;
  minima["mc_norm"] = minimum;
  maxima["mc_norm"] = max;
  
  Double_t random = random_generator->Uniform(0.5,1.5);
  initial_val["mc_norm"] = random;
  initial_err["mc_norm"] = 0.05;

  lh_function.AddDist(*mc_pdf);
  lh_function.SetConstraint("mc_norm", 1., 0.05);
  
  ObsSet  obsSet(0);
  /*if (apply_energy_resolution_convolution) {
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
  }*/

  //if (apply_energy_scaling) {
  Scale* scale = new Scale("scale");
  scale->RenameParameter("scaleFactor","nhitscale_");
  scale->SetScaleFactor(1.);
  scale->SetAxes(axes);
  scale->SetTransformationObs(obsSet);
  scale->SetDistributionObs(obsSet);
  scale->Construct();

  minima["nhitscale_"] = 0.95;//nhit_scaling_estimate - nhit_scaling_estimate_sigma;
  maxima["nhitscale_"] = 1.3;//nhit_scaling_estimate + nhit_scaling_estimate_sigma; 
  initial_val["nhitscale_"] = nhit_scaling_estimate;
  initial_err["nhitscale_"] = nhit_scaling_estimate_sigma;
      
  lh_function.AddSystematic(scale);
  lh_function.SetConstraint("nhitscale_", nhit_scaling_estimate, nhit_scaling_estimate_sigma);
  //}
  
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
  
  /*Convolution* fitconv = new Convolution("fitconv");
  if (apply_energy_resolution_convolution) {
    GaussianERes* fitgaus = new GaussianERes(best_fit.at("eres_"),"fitgaus"); 
    fitconv->SetEResolution(fitgaus);
    fitconv->SetAxes(axes);
    fitconv->SetTransformationObs(obsSet);
    fitconv->SetDistributionObs(obsSet);
    fitconv->Construct();
    }*/
  
  Scale* fitscale = new Scale("fitscale");
  //if (apply_energy_scaling) {
  fitscale->SetScaleFactor(best_fit.at("nhitscale_"));
  fitscale->SetAxes(axes);
  fitscale->SetTransformationObs(obsSet);
  fitscale->SetDistributionObs(obsSet);
  fitscale->Construct();
    //}
  
  //for (ULong64_t j = 0; j < n_pdf; j++){
  //sprintf("mc_norm", "%s_norm", reactor_names[j].c_str());
  //reactor_osc_pdf[j]->Normalise();
  BinnedED pdf;
  pdf = *mc_pdf;
  pdf.Normalise();
  pdf = fitscale->operator()(pdf);
  pdf.Scale(best_fit.at("mc_norm"));
  /*if (apply_energy_resolution_convolution) {
	pdf = fitconv->operator()(*reactor_osc_pdf[j]);
	pdf.Normalise();
	//pdf.Scale(best_fit.at(name));
      }
      if (apply_energy_scaling) {
	if (!apply_energy_resolution_convolution)
	  pdf = fitscale->operator()(*reactor_osc_pdf[j]);
	else
	  
	pdf.Normalise();
	
    }*/
  
  //if (apply_energy_scaling || apply_energy_resolution_convolution){
  //pdf.Scale(best_fit.at("mc_norm"));
  //}
  //if (!apply_energy_scaling && !apply_energy_resolution_convolution){
  
  Double_t lh_val = 99999; // positive non-sensical value to return if fit is not valid
  if (fit_validity == true)
    lh_val = (-1.)*lh_function.Evaluate(); //lh_function.Evaluate();

  // write plots to file (only 'good' plots - those with the best fit values)
  //if (param_d21>=param_d21_plot_min && param_d21<=param_d21_plot_max && param_s12>=param_s12_plot_min && param_s12<=param_s12_plot_max){
  TFile * f_out = new TFile("test_nhitscale.root","RECREATE");
  f_out->cd();
  TCanvas* c_data_fit = new TCanvas("c_data_fit");
  c_data_fit->cd();

  
  // data set
  TH1D data_set_hist = DistTools::ToTH1D(*data_pdf);
  // data_hist.Sumw2();
  data_set_hist.SetName("data_set_hist");
  data_set_hist.GetYaxis()->SetTitle("Counts");
  data_set_hist.GetXaxis()->SetTitle("nhits");
  data_set_hist.SetLineColor(kBlue);
  data_set_hist.Write();
  data_set_hist.Draw("same");

  TH1D mc_set_hist = DistTools::ToTH1D(pdf);
  mc_set_hist.SetName("mc_set_hist");
  mc_set_hist.GetYaxis()->SetTitle("Counts");
  mc_set_hist.GetXaxis()->SetTitle("nhits");
  mc_set_hist.SetLineColor(kRed);
  mc_set_hist.Write();
  mc_set_hist.Draw("same");

  f_out->cd();
  c_data_fit->Write();
  TH1D mc_hist = DistTools::ToTH1D(*mc_pdf);
  mc_hist.Write();

  f_out->Close();
    
  printf("fit valid: %d, lh_value:%.9f\n", fit_validity, lh_val);
  printf("End fit--------------------------------------\n");
  return lh_val;
}

int main(int argc, char *argv[]) {
//int main() {
  /*if (argc != 25){
    std::cout<<"Error: 24 arguments expected."<<std::endl;
      return 1; // return>0 indicates error code
      }
  else{*/
    std::string data_path = "/data/snoplus/PrunedBiPos/output/combined_out_z_257669_259063_norm.root";
    std::string mc_path = "/data/snoplus/PrunedBiPos/output/combined_mc_out_PartialScintBipo214_ScintRun.root";

    printf("Begin--------------------------------------\n");
    
    TFile * f = new TFile(mc_path.c_str());
    //f->ls();
    TH1D* h_mc_nhit_Bi = (TH1D*)f->Get("h_nhit_Bi_r5000");

    const size_t nhit_min = 500;//h_mc_nhit_Bi->GetXaxis()->GetXmin();
    const size_t nhit_max = 900;//h_mc_nhit_Bi->GetXaxis()->GetXmax();
    const size_t n_bins = (nhit_max-nhit_min)/8.;//h_mc_nhit_Bi->GetXaxis()->GetNbins();

    AxisCollection axes;
    axes.AddAxis(BinAxis("nhits", nhit_min, nhit_max, n_bins));

    BinnedED *nhit_mc_pdf = new BinnedED("mc", axes);
    nhit_mc_pdf->SetObservables(0);
    
    for(size_t j = 0; j < h_mc_nhit_Bi->GetXaxis()->GetNbins(); j++){
      double centre = h_mc_nhit_Bi->GetXaxis()->GetBinCenter(j+1);
      double content = h_mc_nhit_Bi->GetBinContent(j+1);
      nhit_mc_pdf->Fill(centre,content);
    }

    //nhit_mc_pdf->Normalise();
    f->Close();
    
    
    TFile * f_data = new TFile(data_path.c_str());
    //f->ls();
    TH1D* h_data_nhit_Bi = (TH1D*)f_data->Get("h_nhit_Bi_r5000_norm");

    BinnedED *nhit_data_pdf = new BinnedED("data", axes);
    nhit_data_pdf->SetObservables(0);
      
    for(size_t j = 0; j < h_data_nhit_Bi->GetXaxis()->GetNbins() ; j++){
      double centre = h_data_nhit_Bi->GetXaxis()->GetBinCenter(j+1);
      double content = h_data_nhit_Bi->GetBinContent(j+1);
      nhit_data_pdf->Fill(centre,content);
    }
    
    //std::cout<<nhit_data_pdf->Integral()<<"  "<<nhit_mc_pdf->Integral()<<std::endl;
    //const bool apply_energy_resolution_convolution = false;
    //const bool apply_energy_scaling = true;
    //double e_resolution_estimate = 0.045;
    //double e_resolution_estimate_sigma = 0.005;
    //const Double_t del_energy = -0.784;  // energy added to antinu MC KE truth to convert to Eprompt truth
    double e_scaling_estimate = 1.0;
    double e_scaling_estimate_sigma = 0.05;


    // save objects to file
    //printf("Save objects to file...\n");
    //TFile *file_out = new TFile(out_filename_plots.c_str(), "RECREATE");
    bool fit_validity = 0;
    ULong64_t fit_try_max = 20;
    ULong64_t print_plots = 0;

    //printf("running: d_21:%.9f(%.9f-%.9f) s_12:%.7f(%.7f-%.7f)\n", d21, param_d21_plot_min, param_d21_plot_max, s12, param_s12_plot_min, param_s12_plot_max);
    //printf("writing plots to: %s\n", out_filename_plots.c_str());
	//print_plots+=1;
   
    //printf("print plots: %d\n\n ", print_plots);
    
    double lh_value = 999999;

    fit_validity = 0;
    /*for (ULong64_t fit_try=1; fit_try<=fit_try_max; fit_try++) {
      lh_value = LHFit_fit(nhit_data_pdf, nhit_mc_pdf,
                           fit_validity, e_scaling_estimate,
                           e_scaling_estimate_sigma,axes);
  
      if (fit_validity==0)
          printf("Fit invalid... retrying (attempt no: %llu)\n", fit_try);
      else{
          printf("Fit valid. (attempt no: %llu)\n", fit_try);
	  fit_try = fit_try_max+1;
      }
      }*/
    ///////////////////

    ///////////////////

    // close unoscillated reactor file
    f_data->Close();

    printf("End--------------------------------------\n");
    return 0; // completed successfully
}
//}
