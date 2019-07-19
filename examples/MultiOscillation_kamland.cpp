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

Double_t LHFit_fit(BinnedED &data_set_pdf, const std::string &spectrum_unosc_filepath,
    std::vector<std::string> &reactor_names, std::vector<Double_t> &distances, std::vector<std::string> &reactor_types,
    std::vector<Double_t> &constraint_means, std::vector<Double_t> &constraint_sigmas,
    TFile *file_out,
    Double_t param_d21, Double_t param_s12, Double_t param_s13,
    bool &fit_validity,
    const double e_min, const double e_max, const size_t n_bins, const double flux_data){

    printf("Begin fit--------------------------------------\n");
    printf("LHFit_fit:: del_21:%.9f, sin2_12:%.7f, sin2_13:%.7f\n", param_d21, param_s12, param_s13);

    char name[1000];
    char spectrum_osc_prompt_filepath[1000];
    const ULong64_t n_pdf = reactor_names.size();

    // set up binning
    ObsSet data_rep(0);
    AxisCollection axes;
    axes.AddAxis(BinAxis("ev_prompt_fit", e_min, e_max, n_bins));

    // create LH function
    BinnedNLLH lh_function;
    lh_function.SetBufferAsOverflow(true);
    int buff = 0;
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

    BinnedED **reactor_unosc_pdf = new BinnedED*[n_pdf];
    BinnedED **reactor_pdf = new BinnedED*[n_pdf];
    
    Double_t constraint_osc_mean_total = 0.;
    std::vector<Double_t> constraint_osc_means;
    Double_t data_set_pdf_integral = data_set_pdf.Integral();

    for (ULong64_t i = 0; i < n_pdf; i++){
        // for each reactor, load spectrum pdf for reactor type
        reactor_pdf[i] = new BinnedED(reactor_names[i], axes);
        reactor_pdf[i]->SetObservables(0);

        sprintf(spectrum_osc_prompt_filepath, "/home/lidgard/antinu_analysis/sensitivity_plot/processed/data/osc_ntp_prompt_%s_%d_%d_%d.root", reactor_names[i].c_str(), (int)(param_d21*1e9), (int)(param_s12*1e7), (int)(param_s13*1e7)); // change this to use the output path of the csv file

        if ((reactor_types[i]=="PWR")||(reactor_types[i]=="BWR")){
            printf("adding pwr reactor: ");
            if (test_file_exists(spectrum_osc_prompt_filepath)==false)
                write_file_pruned(spectrum_unosc_filepath.c_str(), spectrum_osc_prompt_filepath, param_d21, param_s12, param_s13, distances[i]);
        }

        else if (reactor_types[i]=="PHWR"){
            printf("adding phwr reactor: ");
            if (test_file_exists(spectrum_osc_prompt_filepath)==false)
                write_file_pruned(spectrum_unosc_filepath.c_str(), spectrum_osc_prompt_filepath, param_d21, param_s12, param_s13, distances[i]); // currently using PWR, change filename for reactor type!!!!
        }
        else{
            printf("Throw: Reactor doesn't match any loaded type...\n");
            exit(0); // throw std::exception(); //continue;
        }

	    // load newly created oscillated ntuple
        printf("adding unosc reactor: ");
        
        //ROOTNtuple reactor_unosc_ntp(spectrum_unosc_filepath.c_str(), "nt"); // make this work for specific branches
        TFile *f_in = new TFile(spectrum_unosc_filepath.c_str());
        TTree *reactor_unosc_ntp = (TTree*)f_in->Get("nt");
        Double_t ev_energy_p1;
        reactor_unosc_ntp->SetBranchAddress("ev_fit_energy_p1", &ev_energy_p1);
        reactor_unosc_pdf[i] = new BinnedED(reactor_names[i], axes);
        reactor_unosc_pdf[i]->SetObservables(0);
        for(size_t j = 0; j < reactor_unosc_ntp->GetEntries(); j++){
            reactor_unosc_ntp->GetEntry(j);
            reactor_unosc_pdf[i]->Fill(ev_energy_p1);
        }
        f_in->Close();
            //reactor_unosc_pdf[i]->Fill(reactor_unosc_ntp.GetEntry(j));
        //reactor_unosc_pdf[i]->Scale(1./flux_data);

        ROOTNtuple reactor_ntp(spectrum_osc_prompt_filepath, "nt");
        for(size_t j = 0; j < reactor_ntp.GetNEntries(); j++)
            reactor_pdf[i]->Fill(reactor_ntp.GetEntry(j));
        //reactor_pdf[i]->Scale(1./flux_data);
        
        // work out total oscillated integral of constraints
        Double_t normalisation_unosc = reactor_unosc_pdf[i]->Integral();
        Double_t normalisation_reactor = reactor_pdf[i]->Integral();
        Double_t osc_loss = normalisation_reactor/normalisation_unosc;
        
        constraint_osc_means.push_back(constraint_means[i]*osc_loss);
        constraint_osc_mean_total += constraint_means[i]*osc_loss;
    }
    
    Double_t constraint_osc_normalised_total = 0.; 
    for (ULong64_t i = 0; i < n_pdf; i++){
        
        // reactor pdf is the mc spectrum scaled to the kamland unosc spectrum. total events.
        // here we have reactor individual contributions. so scale by contribution to total. find this with our mc spectra.
        // this is the average number of events to the total in mc (before scaling by kamland).

        reactor_pdf[i]->Normalise(); //remove number of events from mc
        Double_t constraint_osc_normalised = constraint_osc_means[i]/constraint_osc_mean_total*data_set_pdf_integral;
        constraint_osc_normalised_total += constraint_osc_normalised;

        Double_t constraint_osc_normalised_sigma = constraint_osc_normalised*constraint_sigmas[i]/constraint_means[i]; //pow(pow(constraint_sigmas[i]/constraint_means[i],2)+pow(constraint_sigmas[n_pdf]/constraint_means[n_pdf],2),0.5)*contribution_factor;

        printf("reactor_names:%s\tconstraint_osc_normalised:%.4f\tconstraint_osc_normalised_sigma:%.4f\tdata_set_pdf_integral:%.4f\tconstraint_osc_normalised_total:%.4f\n", 
            reactor_names[i].c_str(), constraint_osc_normalised, constraint_osc_normalised_sigma, data_set_pdf_integral, constraint_osc_normalised_total);

        // Setting optimisation limits
        sprintf(name, "%s_norm", reactor_names[i].c_str());
        Double_t min = constraint_osc_normalised-2.*constraint_osc_normalised_sigma; // let min and max float within 2 sigma (but constrained later)
        Double_t max = constraint_osc_normalised+2.*constraint_osc_normalised_sigma;
        if (min < 0) min = 0;
        //if (max < 0) max = 1000;
        minima[name] = min;
        maxima[name] = max;
        printf("  added reactor %d/%d: %s, norm: %.3f (min:%.3f max:%.3f) err: %.3f\n", i+1, n_pdf, reactor_names[i].c_str(), constraint_osc_normalised, min, max, constraint_osc_normalised_sigma);
        initial_val[name] = constraint_osc_normalised;
        initial_err[name] = 0.1*constraint_osc_normalised;

        lh_function.AddDist(*reactor_pdf[i]);
        lh_function.SetConstraint(name, constraint_osc_normalised, constraint_osc_normalised_sigma);

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
        remove(spectrum_osc_prompt_filepath);
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
        TH1D pwr_spectrum_hist = DistTools::ToTH1D(*reactor_unosc_pdf[0]);
        // pwr_spectrum_hist.Sumw2();
        sprintf(name, "pwr_spectrum_hist");
        pwr_spectrum_hist.SetName(name);
        pwr_spectrum_hist.GetYaxis()->SetTitle("Counts");
        pwr_spectrum_hist.GetXaxis()->SetTitle("Energy (MeV)");
        pwr_spectrum_hist.Write();

        // TH1D phwr_spectrum_hist = DistTools::ToTH1D(*spectra_ev_pdf[1]);
        // // phwr_spectrum_hist.Sumw2();
        // sprintf(name, "phwr_spectrum_hist");
        // phwr_spectrum_hist.SetName(name);
        // phwr_spectrum_hist.GetYaxis()->SetTitle("Counts");
        // phwr_spectrum_hist.GetXaxis()->SetTitle("Energy (MeV)");
        // phwr_spectrum_hist.Write();
    }

    printf("fit valid: %d, lh_value:%.9f\n", fit_validity, lh_val);
    printf("End fit--------------------------------------\n");
    return lh_val;
}

int main(int argc, char *argv[]) {

    if (argc != 12){
        std::cout<<"Error: 11 arguments expected."<<std::endl;
        return 1; // return>0 indicates error code
    }
    else{
        const std::string &in_path = argv[1];
        const std::string &info_file = argv[2];
        const std::string &spectrum_unosc_filepath = argv[3];
        const std::string &constraints_info_file = argv[4];
        const std::string &parameter_file = argv[5];
        const double flux_data = atoi(argv[6]);
        const std::string &out_filename_plots = argv[7];
        const std::string &out_filename_csv = argv[8];
        const double e_min = atof(argv[9]);
        const double e_max = atof(argv[10]);
        const size_t n_bins = atoi(argv[11]);
        printf("Begin--------------------------------------\n");

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

        // read in parameter information
        std::vector<Double_t> d_21s;
        std::vector<Double_t> s_12s;
        std::vector<Double_t> s_13s;
        readParameterFile(parameter_file, d_21s, s_12s, s_13s);

        const ULong64_t n_pdf = reactor_names.size();
        const ULong64_t n_parameter_sets = d_21s.size();
        Double_t lh_values[n_parameter_sets];

        AxisCollection axes;
        axes.AddAxis(BinAxis("ev_prompt_fit", e_min, e_max, n_bins));
        BinnedED data_set_pdf("data_set_pdf", axes);

        LHFit_initialise_kamland(data_set_pdf, e_min, e_max, n_bins);

        ////save objects to file
        printf("Save objects to file...\n");
        TFile *file_out = 0;
        bool fit_validity = 0;
        size_t fit_try = 0;
        size_t fit_try_max = 10;
        //while ((fit_try <= fit_try_max)&&(fit_validity==0)) {
            fit_try++;
            for (ULong64_t i=0; i<n_parameter_sets; i++) {

                if (d_21s[i]>=7.48e-05 && d_21s[i]<=7.89e-05 && s_12s[i]>=0.333 && s_12s[i]<=0.375)
                    if (file_out==0) {
                        printf("writing plots to: %s\n", out_filename_plots.c_str());
                        file_out = new TFile(out_filename_plots.c_str(), "RECREATE");
                    }

                printf("Fit number: %llu of %llu\n", i+1, n_parameter_sets);
                lh_values[i] = LHFit_fit(data_set_pdf, spectrum_unosc_filepath,
                                reactor_names, distances, reactor_types,
                                constraint_means, constraint_sigmas,
                                file_out,
                                d_21s[i], s_12s[i], s_13s[i],
                                fit_validity, e_min, e_max, n_bins,
                                flux_data);
            }
        //}

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
