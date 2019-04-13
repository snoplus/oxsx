// A simple fit in energy for signal and a background
#include <string>
#include <vector>
#include <math.h>
#include <Rand.h>
#include <fstream>

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <ROOTNtuple.h>

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
#include <TH1D.h>
#include <TRandom3.h>
#include "AntinuUtils.cpp"

void LHFit_fit(BinnedED &data_set_pdf, BinnedED **spectra_pdf, Double_t *reactor_scale, Double_t *reactor_scale_err, std::vector<std::string> &reactor_names, std::vector<Double_t> &distances, std::vector<std::string> &reactor_types, const std::string out_filename, bool fit_validity){

    printf("Begin fit--------------------------------------\n");
    printf("LHFit_fit...\n");

    char name[1000];
    const ULong64_t n_pdf = reactor_names.size();
    ObsSet data_rep(0);
    // set up binning
    AxisCollection axes;
    Double_t e_min = 0.425*2; //2; //*6 for kamland paper #1, *2 for paper #2
    Double_t e_max = 0.425*19; //8;
    Int_t n_bins = 17; //(8-2)*10; //13 for paper #1, 17 for paper #2
    axes.AddAxis(BinAxis("mc_neutrino_energy", e_min, e_max, n_bins));


    // create LH function
    BinnedNLLH lh_function;
    lh_function.SetBufferAsOverflow(true);
    int buff = 2;
    lh_function.SetBuffer(0, buff, buff);
    lh_function.SetDataDist(data_set_pdf); // initialise withe the data set

    // setup max and min ranges
    ParameterDict minima;
    ParameterDict maxima;
    ParameterDict initial_val;
    ParameterDict initial_err;

    double param_d21 = 6.9e-5;//kamland#1=6.9e-5;//us=7.58e-5;
    double param_s12 = 0.5359;
    double param_s13 = 0.02303;
    minima["d21"] = param_d21*0.01;
    maxima["d21"] = param_d21*100;
    initial_val["d21"] = param_d21;
    initial_err["d21"] = 0.1*param_d21;
    minima["s12"] = 0.1;
    maxima["s12"] = 0.5;
    initial_val["s12"] = param_s12;
    initial_err["s12"] = 1*param_s12;

    // setup survival probability
    printf("Setup survival probability...\n");
    SurvProb *surv_prob[n_pdf];
    NuOsc *reactor_systematic[n_pdf];
    BinnedED **reactor_pdf = new BinnedED*[n_pdf];
    for (ULong64_t i = 0; i < n_pdf; i++){
        // for each reactor, load spectrum pdf for reactor type
        reactor_pdf[i] = new BinnedED(reactor_names[i], axes);
        reactor_pdf[i]->SetObservables(0);
        if ((reactor_types[i]=="PWR")||(reactor_types[i]=="BWR"))
            reactor_pdf[i]->Add(*spectra_pdf[0],1); //use PWR pdf
        else if (reactor_types[i]=="PHWR")
                reactor_pdf[i]->Add(*spectra_pdf[1],1); //use PHWR pdf
            else{
               printf("Throw: Reactor doesn't match any loaded type...\n");
               continue;
            }
        //reactor_pdf[i]->Normalise();
        //reactor_pdf[i]->Scale(reactor_scale[i]); // normalise to integral for each reactor

        // setup survival probability
        sprintf(name, "%s_survival", reactor_names[i].c_str());
        surv_prob[i] = new SurvProb(param_d21, param_s12, distances[i], name);
        surv_prob[i]->Setsinsqrtheta13s(param_s13); // manual, fixed setting of theta13
        surv_prob[i]->RenameParameter("delmsqr21_0", "d21"); // rename all parameters to the same for the fit
        surv_prob[i]->RenameParameter("sinsqrtheta12_0", "s12");
        sprintf(name, "%s_systematic", reactor_names[i].c_str());
        reactor_systematic[i] = new NuOsc(name);
        reactor_systematic[i]->SetFunction(surv_prob[i]);
        reactor_systematic[i]->SetAxes(axes);
        reactor_systematic[i]->SetTransformationObs(data_rep);
        reactor_systematic[i]->SetDistributionObs(data_rep);

        // Setting optimisation limits
        //std::cout << " scale" << reactor_scale[i] << " err" << reactor_scale_err[i] << " min" << reactor_scale[i]-reactor_scale_err[i]*reactor_scale[i] << " max" << reactor_scale[i]+reactor_scale_err[i]*reactor_scale[i] << std::endl;
        sprintf(name, "%s_norm", reactor_names[i].c_str());
        Double_t min = reactor_scale[i]-0.196*reactor_scale_err[i]*reactor_scale[i];
        Double_t max = reactor_scale[i]+0.196*reactor_scale_err[i]*reactor_scale[i];
        if (min < 0) min = 0;
        minima[name] = min;
        maxima[name] = max;
        printf("  added reactor %d/%d: %s, norm: %.3f (min:%.3f max:%.3f)\n", i+1, n_pdf, reactor_names[i].c_str(), reactor_scale[i], min, max);
        initial_val[name] = reactor_scale[i];
        initial_err[name] = reactor_scale_err[i]*reactor_scale[i];

        sprintf(name,"group%d",i);
        lh_function.AddSystematic(reactor_systematic[i],name);
        lh_function.AddDist(*reactor_pdf[i],std::vector<std::string>(1,name), true);

        sprintf(name, "%s_norm", reactor_names[i].c_str());
        lh_function.SetConstraint(name, reactor_scale[i], reactor_scale_err[i]*reactor_scale[i]);
    }

    //lh_function.SetConstraint("d21", 6.9e-5, 1e-7); //7.58e-5;
    //lh_function.SetConstraint("s12", 0.45, 0.3); //0.359;

    printf("Built LH function, fitting...\n");
    // fit
    Minuit min;
    min.SetMethod("Migrad");
    min.SetMaxCalls(10000);
    min.SetTolerance(0.01);
    min.SetMinima(minima);
    min.SetMaxima(maxima);
    min.SetInitialValues(initial_val);
    min.SetInitialErrors(initial_err);

    FitResult fit_result = min.Optimise(&lh_function);
    fit_result.SetPrintPrecision(6);
    ParameterDict best_fit = fit_result.GetBestFit();
    fit_result.Print();
    fit_validity = fit_result.GetValid();
    printf("best_fit[d21]:%.4e best_fit[s12]:%.4f\n", best_fit["d21"], best_fit["s12"]);

    // determine the oscillated normalisation (the fit returns normalisation of unoscillated spectra)
    std::vector<double> normalisation_osc = lh_function.GetWorkingNormalisations();
    for (ULong64_t i = 0; i < n_pdf; i++){
        sprintf(name, "%s_norm", reactor_names[i].c_str());
        printf("%s:\tscale:%.4f\tnorm_unosc: %.3f\tnorm_osc: %.4f\n", name, reactor_scale[i], normalisation_osc.at(i), reactor_scale[i]*normalisation_osc.at(i));
    }

    // apply fitted oscillations to reactor pdf's
    SurvProb *surv_prob_fit[n_pdf];
    BinnedED **reactor_pdf_fitosc = new BinnedED*[n_pdf];
    BinnedED reactor_pdf_fitosc_sum("reactor_pdf_fitosc_sum",axes);
    reactor_pdf_fitosc_sum.SetObservables(data_rep);
    
    for (ULong64_t i = 0; i < n_pdf; i++){
        sprintf(name, "%s_pdf_fitosc", reactor_names[i].c_str());
        reactor_pdf_fitosc[i] = new BinnedED(name, axes);
        reactor_pdf_fitosc[i]->SetObservables(data_rep);

        NuOsc OscResult("OscResult");
        sprintf(name, "%s_survival_fit", reactor_names[i].c_str());
        surv_prob_fit[i] = new SurvProb(best_fit.at("d21"), best_fit.at("s12"), distances[i], name);
        surv_prob_fit[i]->Setsinsqrtheta13s(param_s13); // manual, fixed setting of theta13
        OscResult.SetFunction(surv_prob_fit[i]);

        OscResult.SetAxes(axes);
        OscResult.SetTransformationObs(data_rep);
        OscResult.SetDistributionObs(data_rep);
        OscResult.Construct();

        reactor_pdf_fitosc[i]->Add(OscResult(*reactor_pdf[i]),1);

        sprintf(name,"%s_norm",reactor_names[i].c_str());
        reactor_pdf_fitosc[i]->Scale(best_fit.at(name)*normalisation_osc.at(i));

        reactor_pdf_fitosc_sum.Add(*reactor_pdf_fitosc[i],1);
    }

    // save objects to file
    printf("Save objects to file...\n");
    TFile *file_out = new TFile(out_filename.c_str(), "RECREATE");

    // oscillated reactor pdf's
    TH1D *reactor_hist = new TH1D[n_pdf];
    TH1D *reactor_hist_fitosc = new TH1D[n_pdf];
    for (ULong64_t i = 0; i< n_pdf; i++){
        reactor_hist[i] = DistTools::ToTH1D(*reactor_pdf[i]);
        sprintf(name, "%s_hist", reactor_names[i].c_str());
        reactor_hist[i].SetName(name);
        reactor_hist[i].Write();

        reactor_hist_fitosc[i] = DistTools::ToTH1D(*reactor_pdf_fitosc[i]);
        sprintf(name, "%s_hist_fitosc", reactor_names[i].c_str());
        reactor_hist_fitosc[i].SetName(name);
        reactor_hist_fitosc[i].GetXaxis()->SetTitle("Energy (MeV)");
        reactor_hist_fitosc[i].GetYaxis()->SetTitle("Counts");
        reactor_hist_fitosc[i].SetLineColor(kRed);
        reactor_hist_fitosc[i].Write();
    }
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
    TH1D pwr_spectrum_hist = DistTools::ToTH1D(*spectra_pdf[0]);
    // pwr_spectrum_hist.Sumw2();
    sprintf(name, "pwr_spectrum_hist");
    pwr_spectrum_hist.SetName(name);
    pwr_spectrum_hist.GetYaxis()->SetTitle("Counts");
    pwr_spectrum_hist.GetXaxis()->SetTitle("Energy (MeV)");
    pwr_spectrum_hist.Write();

    TH1D phwr_spectrum_hist = DistTools::ToTH1D(*spectra_pdf[1]);
    // phwr_spectrum_hist.Sumw2();
    sprintf(name, "phwr_spectrum_hist");
    phwr_spectrum_hist.SetName(name);
    phwr_spectrum_hist.GetYaxis()->SetTitle("Counts");
    phwr_spectrum_hist.GetXaxis()->SetTitle("Energy (MeV)");
    phwr_spectrum_hist.Write();

    file_out->Close();

    //write output to png image
    TCanvas *c = new TCanvas();
    data_set_hist.Draw();
    reactor_hist_fitosc_sum.Draw("same");
    sprintf(name, "%s.png", out_filename.c_str());
    c->SaveAs(name);

    printf("End fit--------------------------------------\n");
}

int main(int argc, char *argv[]){

    if (argc != 6){
        printf("Error: 5 arguments expected.\n");
        return 1; // return>0 indicates error code
    }
    else{
        const std::string &in_path = argv[1];
        const std::string &data_path = argv[2];
        const std::string &info_file = argv[3];
        const size_t flux_data = atoi(argv[4]);
        const std::string &out_file = argv[5];

        printf("--------------------------------------\n");

        // read in reactor information
        std::vector<std::string> reactor_names;
        std::vector<Double_t> distances;
        std::vector<std::string> reactor_types;
        std::vector<ULong64_t> n_cores;
        std::vector<Double_t> powers;
        std::vector<Double_t> power_errs;
        readInfoFile(info_file, reactor_names, distances, reactor_types, n_cores, powers, power_errs);

        const ULong64_t n_pdf = reactor_names.size();
        BinnedED **spectra_pdf = new BinnedED*[n_pdf]; // PWR=0, PHWR=1
        Double_t *reactor_scale = new Double_t[n_pdf];
        Double_t *reactor_scale_err = new Double_t[n_pdf];
        bool fit_validity = false;

        for (size_t i=0; i<(size_t)reactor_names.size(); i++)
            reactor_scale_err[i] = power_errs[i]/powers[i];

        BinnedED data_set_pdf = LHFit_initialise(spectra_pdf, reactor_scale, in_path, data_path, reactor_names, flux_data);

        LHFit_fit(data_set_pdf, spectra_pdf, reactor_scale, reactor_scale_err, reactor_names, distances, reactor_types, out_file, fit_validity);

        printf("--------------------------------------\n");
        return 0;
    }
}
