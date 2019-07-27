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

Double_t LHFit_fit(BinnedED &data_set_pdf, const std::string &spectrum_phwr_unosc_filepath,
    const std::string &spectrum_pwr_unosc_filepath, const std::string &spectrum_further_unosc_filepath,
    std::vector<std::string> &reactor_names, std::vector<std::string> &reactor_types,
    std::vector<Double_t> &constraint_means, std::vector<Double_t> &constraint_sigmas,
    TFile *file_out,
    Double_t param_d21, Double_t param_s12, Double_t param_s13,
    bool &fit_validity,
    const double e_min, const double e_max, const size_t n_bins,
    const double flux_data, const double mc_scale_factor,
    const double param_d21_plot_min, const double param_d21_plot_max, const double param_s12_plot_min, const double param_s12_plot_max){

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
    int buff = 0;
    lh_function.SetBuffer(0, buff, buff);
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

    for (ULong64_t i = 0; i < n_pdf; i++){
        // for each reactor, load spectrum pdf for reactor type
        sprintf(name, "%s_unosc", reactor_names[i].c_str());
        reactor_unosc_pdf[i] = new BinnedED(name, axes);
        reactor_unosc_pdf[i]->SetObservables(0);
        reactor_osc_pdf[i] = new BinnedED(reactor_names[i], axes);
        reactor_osc_pdf[i]->SetObservables(0);

        if ((reactor_types[i]=="PWR")||(reactor_types[i]=="BWR"))
            sprintf(name, "%s", spectrum_pwr_unosc_filepath.c_str());
        else if (reactor_types[i]=="PHWR")
            sprintf(name, "%s", spectrum_phwr_unosc_filepath.c_str());
        else if (reactor_types[i]=="further_reactors") // if reactor type is the set of further reactors
            sprintf(name, "%s", spectrum_further_unosc_filepath.c_str());
        else{
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
        if (reactor_types[i]=="further_reactors")
            ntOscillate_pruned(reactor_unosc_ntp, reactor_osc_ntp, param_d21, param_s12, param_s13, 10000);
        else
            ntOscillate_pruned(reactor_unosc_ntp, reactor_osc_ntp, param_d21, param_s12, param_s13);

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

        // work out total oscillated integral of constraints
        Double_t normalisation_unosc = reactor_unosc_pdf[i]->Integral();
        Double_t normalisation_reactor = reactor_osc_pdf[i]->Integral();
        Double_t osc_loss = normalisation_reactor/normalisation_unosc;

        Double_t constraint_osc_mean = constraint_means[i]*osc_loss*mc_scale_factor;
        Double_t constraint_osc_sigma = (constraint_sigmas[i]/constraint_means[i])*constraint_osc_mean; //pow(pow(constraint_sigmas[i]/constraint_means[i],2)+pow(constraint_sigmas[n_pdf]/constraint_means[n_pdf],2),0.5)*contribution_factor;

        reactor_osc_pdf[i]->Normalise(); //remove number of events from mc
        reactor_unosc_pdf[i]->Scale(1./flux_data); // osc pdf gets fitted, the unosc doesn't, scale it simply for plotting..

        // Setting optimisation limits
        sprintf(name, "%s_norm", reactor_names[i].c_str());
        Double_t min = constraint_osc_mean-2.*constraint_osc_sigma; // let min and max float within 2 sigma
        Double_t max = constraint_osc_mean+2.*constraint_osc_sigma;
        if (min < 0) min = 0;
        minima[name] = min;
        maxima[name] = max;
        printf("  added reactor %d/%d: %s, norm: %.3f (min:%.3f max:%.3f) err: %.3f data_int:%.0f\n", i+1, n_pdf, reactor_names[i].c_str(), constraint_osc_mean, min, max, constraint_osc_sigma, data_set_pdf_integral);
        Double_t random = random_generator->Uniform(0.5,1.5);
        initial_val[name] = constraint_osc_mean*random;
        initial_err[name] = constraint_osc_sigma;

        lh_function.AddDist(*reactor_osc_pdf[i]);
        lh_function.SetConstraint(name, constraint_osc_mean, constraint_osc_sigma);
    }

    // fit
    printf("Built LH function, fitting...\n");
    Minuit min;
    min.SetMethod("Migrad");
    min.SetMaxCalls(1000000);
    min.SetTolerance(3);
    min.SetMinima(minima);
    min.SetMaxima(maxima);
    min.SetInitialValues(initial_val);
    min.SetInitialErrors(initial_err);

    FitResult fit_result = min.Optimise(&lh_function);
    fit_result.SetPrintPrecision(9);
    ParameterDict best_fit = fit_result.GetBestFit();
    fit_result.Print();
    fit_validity = fit_result.GetValid();

    for (ULong64_t j = 0; j < n_pdf; j++){
        sprintf(name, "%s_norm", reactor_names[j].c_str());
        reactor_osc_pdf[j]->Normalise();
        reactor_osc_pdf[j]->Scale(best_fit.at(name));
        reactor_osc_pdf_fitosc_sum.Add(*reactor_osc_pdf[j]);
    }

    Double_t lh_val = 99999; // positive non-sensical value to return if fit is not valid
    if (fit_validity == true)
        lh_val =(-1)*lh_function.Evaluate();

    // write plots to file (only 'good' plots - those with the best fit values)
    if (param_d21>=param_d21_plot_min && param_d21<=param_d21_plot_max && param_s12>=param_s12_plot_min && param_s12<=param_s12_plot_max){
        file_out->cd();
        // and their sum
        TH1D reactor_hist_fitosc_sum = DistTools::ToTH1D(reactor_osc_pdf_fitosc_sum);
        sprintf(name, "reactor_hist_fitosc_sum_d21%.9f_s12%.7f_s13%.7f", param_d21, param_s12, param_s13);
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

    if (argc != 20){
        std::cout<<"Error: 19 arguments expected."<<std::endl;
        return 1; // return>0 indicates error code
    }
    else{
        const std::string &in_path = argv[1];
        const std::string &data_path = argv[2];
        const std::string &info_file = argv[3];
        const std::string &spectrum_phwr_unosc_filepath = argv[4];
        const std::string &spectrum_pwr_unosc_filepath = argv[5];
        const std::string &spectrum_further_unosc_filepath = argv[6];
        const std::string &constraints_info_file = argv[7];
        const std::string &parameter_file = argv[8];
        const double flux_data = atof(argv[9]);
        const double mc_scale_factor = atof(argv[10]);
        const std::string &out_filename_plots = argv[11];
        const std::string &out_filename_csv = argv[12];
        const double e_min = atof(argv[13]);
        const double e_max = atof(argv[14]);
        const size_t n_bins = atoi(argv[15]);
        const double param_d21_plot_min = atof(argv[16]);
        const double param_d21_plot_max = atof(argv[17]);
        const double param_s12_plot_min = atof(argv[18]);
        const double param_s12_plot_max = atof(argv[19]);
        printf("Begin--------------------------------------\n");

        // read in reactor information
        std::vector<std::string> reactor_names;
        std::vector<Double_t> distances;
        std::vector<std::string> reactor_types;
        std::vector<ULong64_t> n_cores;
        std::vector<Double_t> powers;
        std::vector<Double_t> power_errs;
        readInfoFile(info_file, reactor_names, distances, reactor_types, n_cores, powers, power_errs);
        //reactor_names.push_back("reactors_greater1000km"); //comment these 2 lines out to remove the further reactors
        //reactor_types.push_back("further_reactors");

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

        // initialise data
        LHFit_initialise(data_set_pdf, data_path, flux_data, e_min, e_max, n_bins);

        ////save objects to file
        printf("Save objects to file...\n");
        TFile *file_out = new TFile(out_filename_plots.c_str(), "RECREATE");
        bool fit_validity = 0;
        ULong64_t fit_try_max = 20;
        ULong64_t print_plots = 0;

        for (ULong64_t i=0; i<n_parameter_sets; i++) {

            printf("running: d_21:%.9f(%.9f-%.9f) s_12:%.7f(%.7f-%.7f)\n", d_21s[i], param_d21_plot_min, param_d21_plot_max, s_12s[i], param_s12_plot_min, param_s12_plot_max);
            if (d_21s[i]>=param_d21_plot_min && d_21s[i]<=param_d21_plot_max && s_12s[i]>=param_s12_plot_min && s_12s[i]<=param_s12_plot_max){
                printf("writing plots to: %s\n", out_filename_plots.c_str());
                print_plots++;
            }

            printf("Fit number: %llu of %llu\n", i+1, n_parameter_sets);
            fit_validity = 0;
            for (ULong64_t fit_try=1; fit_try<=fit_try_max; fit_try++) {
                lh_values[i] = LHFit_fit(data_set_pdf, spectrum_phwr_unosc_filepath, spectrum_pwr_unosc_filepath, spectrum_further_unosc_filepath,
                                reactor_names, reactor_types,
                                constraint_means, constraint_sigmas,
                                file_out,
                                d_21s[i], s_12s[i], s_13s[i],
                                fit_validity, e_min, e_max, n_bins,
                                flux_data, mc_scale_factor,
                                param_d21_plot_min, param_d21_plot_max, param_s12_plot_min, param_s12_plot_max);

                if (fit_validity==0)
                    printf("Fit invalid... retrying (attempt no: %llu)\n", fit_try);
                else{
                    printf("Fit valid. (attempt no: %llu)\n", fit_try);
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
