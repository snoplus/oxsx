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
#include <ParameterDict.h>
#include <ContainerTools.hpp>
#include <NuOsc.h>
#include <SurvProb.h>
#include "AntinuUtils.cpp"
#include "../util/oscillate_util.cpp"

Double_t LHFit_fit(BinnedED &data_set_pdf, BinnedED **spectra_ke_pdf, BinnedED **spectra_ev_pdf, std::vector<std::string> &reactor_names, std::vector<Double_t> &distances, std::vector<std::string> &reactor_types, std::vector<Double_t> &fit_means, std::vector<Double_t> &fit_mean_errs, std::vector<Double_t> &fit_sigmas, std::vector<Double_t> &fit_sigma_errs, TFile *file_out, Double_t param_d21, Double_t param_s12, Double_t param_s13, bool &fit_validity){

    printf("Begin fit--------------------------------------\n");
    printf("LHFit_fit:: del_21:%.9f, sin2_12:%.7f, sin2_13:%.7f\n", param_d21, param_s12, param_s13);

    char name[1000];
    char spectrum_unosc_filepath[1000];
    char spectrum_ke_filepath[1000];
    char spectrum_prompt_filepath[1000];
    const ULong64_t n_pdf = reactor_names.size()-1;
    ObsSet data_rep(0);
    // set up binning
    AxisCollection axes;
    Double_t e_min = 0.5; //2; //*6 for kamland paper #1, *2 for paper #2
    Double_t e_max = 8; //8;
    Int_t n_bins = (8-0.5)*10; //(8-2)*10; //13 for paper #1, 17 for paper #2
    axes.AddAxis(BinAxis("mc_neutrino_energy", e_min, e_max, n_bins));

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

    TH1D *reactor_hist = new TH1D[n_pdf];

    BinnedED reactor_pdf_fitosc_sum("reactor_pdf_fitosc_sum",axes);
    reactor_pdf_fitosc_sum.SetObservables(data_rep);

    // setup survival probability
    printf("Setup survival probability...\n");
    //SurvProb *surv_prob[n_pdf];
    BinnedED **reactor_osc_pdf = new BinnedED*[n_pdf];
    BinnedED **reactor_pdf = new BinnedED*[n_pdf];
    for (ULong64_t i = 0; i < n_pdf; i++){
        // for each reactor, load spectrum pdf for reactor type
        reactor_pdf[i] = new BinnedED(reactor_names[i], axes);
        reactor_pdf[i]->SetObservables(0);

        sprintf(spectrum_unosc_filepath, "/data/snoplusmc/lidgard/OXSX_rat6169_penergy_un/flux1/combinedpwr_flux1_year100_cleanround1.root");
        sprintf(spectrum_ke_filepath, "/home/lidgard/antinu_analysis/sensitivity_plot/processed/data/osc_ntp_ke_%s_%d_%d_%d.root", reactor_names[i].c_str(), (int)(param_d21*1e9), (int)(param_s12*1e7), (int)(param_s13*1e7));
        sprintf(spectrum_prompt_filepath, "/home/lidgard/antinu_analysis/sensitivity_plot/processed/data/osc_ntp_prompt_%s_%d_%d_%d.root", reactor_names[i].c_str(), (int)(param_d21*1e9), (int)(param_s12*1e7), (int)(param_s13*1e7));

        if ((reactor_types[i]=="PWR")||(reactor_types[i]=="BWR")){
            printf("adding pwr reactor: ");
            if ((test_file_exists(spectrum_ke_filepath)==false)||(test_file_exists(spectrum_prompt_filepath)==false))
                write_file_pruned(spectrum_unosc_filepath, spectrum_ke_filepath, spectrum_prompt_filepath, param_d21, param_s12, param_s13, distances[i]);
        }

        else if (reactor_types[i]=="PHWR"){
            printf("adding phwr reactor: ");

            if ((test_file_exists(spectrum_ke_filepath)==false)||(test_file_exists(spectrum_prompt_filepath)==false))
                write_file_pruned(spectrum_unosc_filepath, spectrum_ke_filepath, spectrum_prompt_filepath, param_d21, param_s12, param_s13, distances[i]); // currently using PWR!!!!
        }
        else{
            printf("Throw: Reactor doesn't match any loaded type...\n");
            exit(0); // throw std::exception(); //continue;
        }

	    // load newly created oscillated ntuple
        ULong64_t flux_mc = 100;
        printf("adding unosc reactor: ");
        sprintf(name, "/data/snoplusmc/lidgard/OXSX_kamland_rat6169_penergy_un/flux1/combinedpwr_flux1_year100_cleanround1_prompt_oxsx.root");
        ROOTNtuple reactor_osc_ntp(name, "nt");
        reactor_osc_pdf[i] = new BinnedED(reactor_names[i], axes);
        reactor_osc_pdf[i]->SetObservables(0);
        for(size_t j = 0; j < reactor_osc_ntp.GetNEntries(); j++)
            reactor_osc_pdf[i]->Fill(reactor_osc_ntp.GetEntry(j));
        reactor_osc_pdf[i]->Scale(1./flux_mc);
        
        ROOTNtuple reactor_ntp(spectrum_prompt_filepath, "nt");
        for(size_t j = 0; j < reactor_ntp.GetNEntries(); j++)
            reactor_pdf[i]->Fill(reactor_ntp.GetEntry(j));
        reactor_pdf[i]->Scale(1./flux_mc);

        // reactor pdf is the mc spectrum scaled to the kamland unosc spectrum. total events.
        // here we have reactor individual contributions. so scale by contribution to total. find this with our mc spectra.
        // this is the average number of events to the total in mc (before scaling by kamland).

        Double_t normalisation_unosc = reactor_osc_pdf[i]->Integral();
        Double_t normalisation_reactor = reactor_pdf[i]->Integral();
        Double_t normalisation_loss = normalisation_reactor/normalisation_unosc;
	    Double_t contribution_factor = fit_means[i]/fit_means[n_pdf];
        Double_t contribution_factor_err = pow(pow(fit_sigmas[i]/fit_means[i],2)+pow(fit_sigmas[n_pdf]/fit_means[n_pdf],2),0.5)*contribution_factor;

        reactor_pdf[i]->Normalise(); //remove number of events from mc
	    reactor_pdf[i]->Scale(scale_factor_unosc); // scale to kamland
        reactor_pdf[i]->Scale(contribution_factor); // scale to individual reactor
        reactor_pdf[i]->Scale(normalisation_loss); // scale to individual reactor

        printf("%s:\tscale:%.4f\tscale_combined:%.4f\tc_factor: %.4f\tc_factor_err: %.4f\tnorm_reactor: %.0f\tnorm_unosc: %.0f\tn_factor: %.4f\tnorm: %.1f\n", reactor_names[i].c_str(), fit_means[i], fit_means[n_pdf], contribution_factor, contribution_factor_err, normalisation_reactor, normalisation_unosc, normalisation_loss, reactor_pdf[i]->Integral());

        // Setting optimisation limits
        sprintf(name, "%s_norm", reactor_names[i].c_str());
        Double_t min = scale_factor_unosc*normalisation_loss*(contribution_factor-2.*contribution_factor_err); // let min and max float within 2 sigma (but constrained later)
        Double_t max = scale_factor_unosc*normalisation_loss*(contribution_factor+2.*contribution_factor_err);
        if (min < 0) min = 0;
        minima[name] = min;
        maxima[name] = max;
        printf("  added reactor %d/%d: %s, norm: %.3f (min:%.3f max:%.3f) err: %.3f\n", i+1, n_pdf, reactor_names[i].c_str(), scale_factor_unosc*normalisation_loss*contribution_factor, min, max, scale_factor_unosc*normalisation_loss*contribution_factor_err);
        initial_val[name] = scale_factor_unosc*normalisation_loss*contribution_factor;
        initial_err[name] = scale_factor_unosc*normalisation_loss*contribution_factor_err;

        lh_function.AddDist(*reactor_pdf[i]);
        lh_function.SetConstraint(name, scale_factor_unosc*normalisation_loss*contribution_factor, scale_factor_unosc*normalisation_loss*contribution_factor_err);

        reactor_pdf_fitosc_sum.Add(*reactor_pdf[i]);

        // create histograms to save
        if (param_d21>=7.48e-05 && param_d21<=7.89e-05 && param_s12>=0.333 && param_s12<=0.375){
            reactor_hist[i] = DistTools::ToTH1D(*reactor_pdf[i]);
            sprintf(name, "%s_hist_d21%.5f_s12%.5f_s13%.5f", reactor_names[i].c_str(), param_d21, param_s12, param_s13);
            reactor_hist[i].SetName(name);
            reactor_hist[i].GetXaxis()->SetTitle("Energy (MeV)");
            reactor_hist[i].GetYaxis()->SetTitle("Counts");
            reactor_hist[i].SetLineColor(kRed);
            file_out->cd();
            reactor_hist[i].Write();
        }

        //delete oscillated ntp root files
        remove(spectrum_prompt_filepath);
        remove(spectrum_ke_filepath);
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

    Double_t lh_val = 99999; // positive non-sensical value to return if fit is not valid
    if (fit_validity == true)
        lh_val =(-1)*lh_function.Evaluate();


    // write plots to file (only 'good' plots - those with the best fit values)
    if (param_d21>=7.48e-05 && param_d21<=7.89e-05 && param_s12>=0.333 && param_s12<=0.375){
        file_out->cd();
        // and their sum
        TH1D reactor_hist_fitosc_sum = DistTools::ToTH1D(reactor_pdf_fitosc_sum);
        sprintf(name, "reactor_hist_fitosc_sum");
        reactor_hist_fitosc_sum.SetName(name);
        reactor_hist_fitosc_sum.GetXaxis()->SetTitle("Energy (MeV)");
        reactor_hist_fitosc_sum.GetYaxis()->SetTitle("Counts");
        reactor_hist_fitosc_sum.SetLineColor(kRed);
        reactor_hist_fitosc_sum.Write();

        // data set
        TH1D data_set_hist = DistTools::ToTH1D(data_set_pdf);
        // data_hist.Sumw2();
        sprintf(name, "data_set_hist");
        data_set_hist.SetName(name);
        data_set_hist.GetYaxis()->SetTitle("Counts");
        data_set_hist.GetXaxis()->SetTitle("Energy (MeV)");
        data_set_hist.Write();

        // pdfs of spectra
        TH1D pwr_spectrum_hist = DistTools::ToTH1D(*reactor_osc_pdf[0]);
        // pwr_spectrum_hist.Sumw2();
        sprintf(name, "pwr_spectrum_hist");
        pwr_spectrum_hist.SetName(name);
        pwr_spectrum_hist.GetYaxis()->SetTitle("Counts");
        pwr_spectrum_hist.GetXaxis()->SetTitle("Energy (MeV)");
        pwr_spectrum_hist.Write();

        TH1D phwr_spectrum_hist = DistTools::ToTH1D(*spectra_ev_pdf[1]);
        // phwr_spectrum_hist.Sumw2();
        sprintf(name, "phwr_spectrum_hist");
        phwr_spectrum_hist.SetName(name);
        phwr_spectrum_hist.GetYaxis()->SetTitle("Counts");
        phwr_spectrum_hist.GetXaxis()->SetTitle("Energy (MeV)");
        phwr_spectrum_hist.Write();
    }

    printf("fit valid: %d, lh_value:%.9f\n", fit_validity, lh_val);
    printf("End fit--------------------------------------\n");
    return lh_val;
}

int main(int argc, char *argv[]) {

    if (argc != 9){
        std::cout<<"Error: 8 arguments expected."<<std::endl;
        return 1; // return>0 indicates error code
    }
    else{
        const std::string &in_path = argv[1];
        const std::string &data_path = argv[2];
        const std::string &info_file = argv[3];
        const std::string &constraints_info_file = argv[4];
        const std::string &parameter_file = argv[5];
        const size_t flux_data = atoi(argv[6]);
        const std::string &out_filename_plots = argv[7];
        const std::string &out_filename_csv = argv[8];
        printf("Begin--------------------------------------\n");

        // read in reactor information
        std::vector<std::string> reactor_names;
        std::vector<Double_t> distances;
        std::vector<std::string> reactor_types;
        std::vector<ULong64_t> n_cores;
        std::vector<Double_t> powers;
        std::vector<Double_t> power_errs;
        readInfoFile(info_file, reactor_names, distances, reactor_types, n_cores, powers, power_errs);
        reactor_names.push_back("combinedpwr");

        // read in constraint information
        std::vector<Double_t> fit_means;
        std::vector<Double_t> fit_mean_errs;
        std::vector<Double_t> fit_sigmas;
        std::vector<Double_t> fit_sigma_errs;

        // read constraint info for each reactor in the info file (one at time to ensure they match correctly)
        for (size_t i=0; i<(size_t)reactor_names.size(); i++){
            double fit_mean, fit_mean_err, fit_sigma, fit_sigma_err;
            readConstraintsInfoFile(constraints_info_file, reactor_names[i].c_str(), fit_mean, fit_mean_err, fit_sigma, fit_sigma_err);
            fit_means.push_back(fit_mean);
            fit_mean_errs.push_back(fit_mean_err);
            fit_sigmas.push_back(fit_sigma);
            fit_sigma_errs.push_back(fit_sigma_err);
        }

        for (size_t i=0; i<(size_t)reactor_names.size(); i++)
            printf("i:%llu, reactor_name:%s, fit_mean: %.3f, fit_mean_err: %.3f, fit_sigma: %.3f, fit_sigma_err: %.3f\n", i, reactor_names[i].c_str(), fit_means[i], fit_mean_errs[i], fit_sigmas[i], fit_sigma_errs[i]);

        // read in parameter information
        std::vector<Double_t> d_21s;
        std::vector<Double_t> s_12s;
        std::vector<Double_t> s_13s;
        readParameterFile(parameter_file, d_21s, s_12s, s_13s);

        const ULong64_t n_pdf = reactor_names.size()-1;
        BinnedED **spectra_ke_pdf = new BinnedED*[n_pdf]; // PWR=0, PHWR=1
        BinnedED **spectra_ev_pdf = new BinnedED*[n_pdf]; // PWR=0, PHWR=1
        const ULong64_t n_parameter_sets = d_21s.size();
        Double_t lh_values[n_parameter_sets];

        BinnedED data_set_pdf = LHFit_initialise(spectra_ke_pdf, spectra_ev_pdf, in_path, data_path, flux_data);

        // save objects to file
        printf("Save objects to file...\n");
        TFile *file_out = 0;

        bool fit_validity = 0;
        for (ULong64_t i=0; i<n_parameter_sets; i++) {

            if (d_21s[i]>=7.48e-05 && d_21s[i]<=7.89e-05 && s_12s[i]>=0.333 && s_12s[i]<=0.375)
                if (file_out==0) {
                    printf("writing plots to: %s\n", out_filename_plots.c_str());
                    file_out = new TFile(out_filename_plots.c_str(), "RECREATE");
                }

            printf("Fit number: %llu of %llu\n", i+1, n_parameter_sets);
            lh_values[i] = LHFit_fit(data_set_pdf, spectra_ke_pdf, spectra_ev_pdf, reactor_names, distances, reactor_types, fit_means, fit_mean_errs, fit_sigmas, fit_sigma_errs, file_out, d_21s[i], s_12s[i], s_13s[i], fit_validity);
        }

        if (file_out!=0) file_out->Close();

        //Write fit coefficients to txt file
        printf("writing to: %s\n", out_filename_csv.c_str());
        FILE *fOut = fopen(out_filename_csv.c_str(), "w");
        fprintf(fOut,"d21,s12,s13,lh_value,fitValidity\n");
        for (ULong64_t i=0; i<n_parameter_sets; i++)
            fprintf(fOut,"%.9f,%.7f,%.7f,%.9f,%llu\n", d_21s[i], s_12s[i], s_13s[i], lh_values[i], fit_validity);
        fclose(fOut);

        printf("End--------------------------------------\n");
        return 0; // completed successfully
    }
}
