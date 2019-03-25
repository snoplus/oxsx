// A fit in energy for signal and a background
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

void readInfoFile(const std::string &runInfoFileName, std::vector<std::string> &reactor_names, std::vector<Double_t> &distances, std::vector<std::string> &reactor_types, std::vector<ULong64_t> &n_cores, std::vector<Double_t> &powers, std::vector<Double_t> &power_errs ) {
    // Read couchDB run-info text file
    std::ifstream in;
    in.open(runInfoFileName.c_str());
    //std::cout << "opening file: " << runInfoFileName.c_str() << std::endl;

    std::fill(reactor_names.begin(), reactor_names.end(), "");
    std::fill(distances.begin(), distances.end(), 0.);
    std::fill(reactor_types.begin(), reactor_types.end(), "");
    std::fill(n_cores.begin(), n_cores.end(), 0);
    std::fill(powers.begin(), powers.end(), 0.);
    std::fill(power_errs.begin(), power_errs.end(), 0.);

    std::string reactor_name,distance,reactor_type,n_core,power,power_err;
    ULong64_t line_no = 0;

    // read until end of file.
    while(in.good()){
        std::getline(in,reactor_name,',');
        std::getline(in,distance,',');
        std::getline(in,reactor_type,',');
        std::getline(in,n_core,',');
        std::getline(in,power,',');
        std::getline(in,power_err,'\n');

        if (line_no>0){ //skip csv header
            if (strcmp(reactor_name.c_str(),"")!=0) {
                reactor_names.push_back(reactor_name);
                distances.push_back(atof(distance.c_str()));
                reactor_types.push_back(reactor_type.c_str());
                n_cores.push_back(atoi(n_core.c_str()));
                powers.push_back(atof(power.c_str()));
                power_errs.push_back(atof(power_err.c_str()));

                //std::cout << "v: reactor_name: " << reactor_names[line_no-1] << ", distance: " << distances[line_no-1] << ", reactor_type: " << reactor_types[line_no-1] << ", n_core: " << n_cores[line_no-1] << ", power: " << powers[line_no-1] << ", power_err: " << power_errs[line_no-1] << std::endl; //debug check ('-1' for header)
            }
        }
        line_no++;
    }
    in.close();

    // print out read info
    //for (size_t i=0; i<(size_t)reactor_names.size(); i++)
    //    printf("i:%llu,reactor_names[i]:%s, distance: %.5f, type: %s, n_cores: %llu, power: %.5f, power_err: %.5f \n", i, reactor_names[i].c_str(), distances[i], reactor_types[i].c_str(), n_cores[i], powers[i], power_errs[i]);
}

void readParameterFile(const std::string &runParameterFileName, std::vector<Double_t> &d_21s, std::vector<Double_t> &s_12s, std::vector<Double_t> &s_13s) {
    // Read text file containing oscillation parameters
    std::ifstream in;
    in.open(runParameterFileName.c_str());
    //std::cout << "opening file: " << runParameterFileName.c_str() << std::endl;

    std::fill(d_21s.begin(), d_21s.end(), 0.);
    std::fill(s_12s.begin(), s_12s.end(), 0.);
    std::fill(s_13s.begin(), s_13s.end(), 0.);

    std::string d_21, s_12, s_13;
    ULong64_t line_no = 0;

    // read until end of file.
    // format of file: 'd_21,s_12,s_13\n'
    while(in.good()){
        std::getline(in,d_21,',');
        std::getline(in,s_12,',');
        std::getline(in,s_13,'\n');

        if (line_no>0){ //skip csv header
            if (strcmp(d_21.c_str(),"")!=0) {
                d_21s.push_back(atof(d_21.c_str()));
                s_12s.push_back(atof(s_12.c_str()));
                s_13s.push_back(atof(s_13.c_str()));

                //std::cout << "v: d_21: " << d_21s[line_no-1] << ", s_12: " << s_12s[line_no-1] << ", s_13: " << s_13s[line_no-1] << std::endl; //debug check ('-1' for header)
            }
        }
        line_no++;
    }
    in.close();

    // // print out read parameters
    // for (size_t i=0; i<(size_t)d_21s.size(); i++)
        // printf("i:%llu, d_21:%.5f, s_12:%.5f, s_13:%.5f\n", i, d_21s[i], s_12s[i], s_13s[i]);
}

BinnedED LHFit_initialise(BinnedED **spectra_pdf, Double_t *reactor_scale, const std::string &in_path, const std::string &data_path, std::vector<std::string> &reactor_names, ULong64_t flux_data){
    
    printf("Begin init--------------------------------------\n");
    printf("LHFit_initialise...\n");

    char name[1000];
    ULong64_t flux_mc = 100;

    // setup spectra filenames
    const std::string pwr_file = std::string("combinedpwr");//std::string("ohi");//std::string("combinedpwr");
    const std::string phwr_file = std::string("combinedphwr");//std::string("ohi");//std::string("combinedphwr");
    printf("Reading individual reactor spectra from: %s\n", in_path.c_str());
    sprintf(name, "%sflux%llu/[reactor_name]_flux%llu_day360_cleanround1_ke_oxsx.root", in_path.c_str(), flux_data, flux_data);
    printf("Using unoscillated spectrum for reactors:\n\t%s\n",name);
    printf("Using unoscillated spectrum for:\n");
    sprintf(name, "%sflux%llu/%s_flux%llu_day360_cleanround1_ke_oxsx.root", in_path.c_str(), flux_mc, pwr_file.c_str(), flux_mc);
    printf("\ttype PWR & BWR: %s\n", name);
    const std::string pwr_unosc_path = std::string(name);
    sprintf(name, "%sflux%llu/%s_flux%llu_day360_cleanround1_ke_oxsx.root", in_path.c_str(), flux_mc, phwr_file.c_str(), flux_mc);
    printf("\ttype PHWR: %s\n", name);
    const std::string phwr_unosc_path = std::string(name);    
    
    // setup (oscillated) data filename
    sprintf(name, "%s", data_path.c_str());
    printf("Loading data spectrum: %s\n\n", name);

    // setup ntuple
    ObsSet data_rep(0);

    // set up binning
    AxisCollection axes;
    Double_t e_min = 2;//0.425*2; //2; //*6 for kamland paper #1, *2 for paper #2
    Double_t e_max = 8;//0.425*19; //8;
    Int_t n_bins = 60;//17; //(8-2)*10; //13 for paper #1, 17 for paper #2
    axes.AddAxis(BinAxis("mc_neutrino_energy", e_min, e_max, n_bins));

    // load (oscillated) data ntuple
    BinnedED data_set_pdf("dataSetPdf", axes);
    data_set_pdf.SetObservables(data_rep);
    ROOTNtuple data_ntp(data_path, "nt");
    for(ULong64_t i = 0; i < data_ntp.GetNEntries(); i++)
        data_set_pdf.Fill(data_ntp.GetEntry(i));
    data_set_pdf.Scale(1./flux_data);
    Double_t data_set_integral = data_set_pdf.Integral(); //record the integral before normalising
    //data_set_pdf.Normalise();
    printf("Loading data: %s (osc)integral:%.3f\n", data_path.c_str(), data_set_pdf.Integral());

    // number of reactors
    const ULong64_t n_pdf = reactor_names.size();

    // load unoscillated reactor ntuples (for scale factors)
    Double_t reactor_integrals[n_pdf];
    for(ULong64_t i = 0; i < n_pdf; i++){
        sprintf(name, "%sflux%llu/%s_flux%llu_day360_cleanround1_ke_oxsx.root", in_path.c_str(), flux_mc, reactor_names[i].c_str(), flux_mc);
        ROOTNtuple reactor_ntp(name, "nt");
        sprintf(name, "%s_pdf", reactor_names[i].c_str());
        BinnedED *reactor_pdf = new BinnedED(name, axes);
        reactor_pdf->SetObservables(0);
        for(ULong64_t j = 0; j < reactor_ntp.GetNEntries(); j++)
            reactor_pdf->Fill(reactor_ntp.GetEntry(j));
        reactor_pdf->Scale(1./flux_data);
        reactor_integrals[i] = reactor_pdf->Integral(); //record the integral
        reactor_scale[i] = reactor_integrals[i];
    }
    
    for(ULong64_t i = 0; i < n_pdf; i++){
    	printf("Loading reactor: %s\t(unosc)integral:%.3f\n", reactor_names[i].c_str(), reactor_scale[i]);
    }

    // load unoscillated spectra
    // pwr spectrum
    ROOTNtuple pwr_spectrum_ntp(pwr_unosc_path, "nt");
    sprintf(name, "pwr_spectrum_pdf");
    spectra_pdf[0] = new BinnedED(name, axes);
    spectra_pdf[0]->SetObservables(0);
    for(ULong64_t i = 0; i < pwr_spectrum_ntp.GetNEntries(); i++)
        spectra_pdf[0]->Fill(pwr_spectrum_ntp.GetEntry(i));
    spectra_pdf[0]->Scale(1./flux_mc);
    Double_t pwr_spectrum_integral = spectra_pdf[0]->Integral(); //record the integral before normalising
    spectra_pdf[0]->Normalise();

    // phwr spectrum
    ROOTNtuple phwr_spectrum_ntp(phwr_unosc_path, "nt");
    sprintf(name, "phwr_spectrum_pdf");
    spectra_pdf[1] = new BinnedED(name, axes);
    spectra_pdf[1]->SetObservables(0);
    for(ULong64_t i = 0; i < phwr_spectrum_ntp.GetNEntries(); i++)
        spectra_pdf[1]->Fill(phwr_spectrum_ntp.GetEntry(i));
    spectra_pdf[1]->Scale(1./flux_mc);
    Double_t phwr_spectrum_integral = spectra_pdf[1]->Integral(); //record the integral before normalising
    spectra_pdf[1]->Normalise();

    printf("End init--------------------------------------\n");
    return data_set_pdf;
}

Double_t LHFit_fit(BinnedED &data_set_pdf, BinnedED **spectra_pdf, BinnedNLLH &lh_function, SurvProb **surv_prob, NuOsc **reactor_systematic, BinnedED **reactor_pdf, Double_t *reactor_scale, Double_t *reactor_scale_err, std::vector<std::string> &reactor_names, std::vector<Double_t> &distances, std::vector<std::string> &reactor_types, Double_t param_d21, Double_t param_s12, Double_t param_s13, bool &fit_validity, ParameterDict &best_fit){

    printf("Begin fit--------------------------------------\n");
    printf("LHFit_fit:: del_21:%.7f, sin2_12:%.5f, sin2_13:%.5f\n", param_d21, param_s12, param_s13);

    char name[1000];
    const ULong64_t n_pdf = reactor_names.size();
    ObsSet data_rep(0);
    // set up binning
    AxisCollection axes;
    Double_t e_min = 2;//0.425*2; //2; //*6 for kamland paper #1, *2 for paper #2
    Double_t e_max = 8;//0.425*19; //8;
    Int_t n_bins = 60;//17; //(8-2)*10; //13 for paper #1, 17 for paper #2
    axes.AddAxis(BinAxis("mc_neutrino_energy", e_min, e_max, n_bins));

    // create LH function
    //BinnedNLLH lh_function;
    lh_function.SetBufferAsOverflow(true);
    int buff = 2;
    lh_function.SetBuffer(0, buff, buff);
    lh_function.SetDataDist(data_set_pdf); // initialise withe the data set
    
    // setup max and min ranges
    ParameterDict minima;
    ParameterDict maxima;
    ParameterDict initial_val;
    ParameterDict initial_err;
    
    //double param_d21 = 6.9e-5;//kamland#1=6.9e-5;//us=7.58e-5;
    //double param_s12 = 0.5359;
    //double param_s13 = 0.02303;
    minima["d21"] = param_d21*0.01;
    maxima["d21"] = param_d21*100;
    initial_val["d21"] = param_d21;
    initial_err["d21"] = 1*param_d21;
    minima["s12"] = 0.2;
    maxima["s12"] = 0.4;
    initial_val["s12"] = param_s12;
    initial_err["s12"] = 1*param_s12;

    // setup survival probability
    printf("Setup survival probability...\n");
    //SurvProb *surv_prob[n_pdf];
    //NuOsc *reactor_systematic[n_pdf];
    //BinnedED **reactor_pdf = new BinnedED*[n_pdf];
std::vector<double> normalisation_osc;
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
        reactor_pdf[i]->Normalise();
        reactor_pdf[i]->Scale(reactor_scale[i]); // normalise to integral for each reactor

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
        Double_t min = 0;//reactor_scale[i]-1.96*reactor_scale_err[i]*reactor_scale[i];
        Double_t max = 10;//reactor_scale[i]+1.96*reactor_scale_err[i]*reactor_scale[i];
        if (min < 0) min = 0;
        minima[name] = min;
        maxima[name] = max;
        printf("\tadded reactor %d/%d: %s\tnorm: %.3f+/-%.3f\t(min:%.3f max:%.3f)\n", i+1, n_pdf, reactor_names[i].c_str(), reactor_scale[i], reactor_scale_err[i]*reactor_scale[i], min, max);
        initial_val[name] = reactor_scale[i];
        initial_err[name] = reactor_scale_err[i]*reactor_scale[i];

        sprintf(name,"group%d",i);
        lh_function.AddSystematic(reactor_systematic[i], name);
        lh_function.AddDist(*reactor_pdf[i], std::vector<std::string>(1, name),true);  

        normalisation_osc = lh_function.GetWorkingNormalisations();
        printf("size:%d\tnpdf:%d\n", name, normalisation_osc.size(), n_pdf);

        sprintf(name, "%s_norm", reactor_names[i].c_str());
        lh_function.SetConstraint(name, reactor_scale[i], reactor_scale_err[i]*reactor_scale[i]);
    }

    //lh_function.SetConstraint("d21", 7e-5, 9e-5); // no constraints for likelihood map
    //lh_function.SetConstraint("s12", 0.5, 0.3);

    // fit
    printf("Built LH function, fitting...\n");
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
    best_fit = fit_result.GetBestFit();
    fit_result.Print();
    fit_validity = fit_result.GetValid();
    printf("bestFit[d21]:%.4e bestFit[s12]:%.4f\n", best_fit["d21"], best_fit["s12"]);

    Double_t lh_val = 99999; // positive non-sensical value to return if fit is not valid
    if (fit_validity == true)
        lh_val =(-1)*lh_function.Evaluate();

    printf("fit valid: %d, lh_value:%.5f\n", fit_validity, lh_val);
    printf("End fit--------------------------------------\n");
    return lh_val;
}

void LHFit_produce_histograms(BinnedED &data_set_pdf, BinnedED **spectra_pdf, BinnedNLLH &lh_function, SurvProb **surv_prob, BinnedED **reactor_pdf, Double_t *reactor_scale, std::vector<std::string> &reactor_names, std::vector<Double_t> &distances, ParameterDict &best_fit, const std::string out_filename){

    char name[1000];
    const ULong64_t n_pdf = reactor_names.size();
    ObsSet data_rep(0);
    // set up binning
    AxisCollection axes;
    Double_t e_min = 2;//0.425*2; //2; //*6 for kamland paper #1, *2 for paper #2
    Double_t e_max = 8;//0.425*19; //8;
    Int_t n_bins = 60;//17; //(8-2)*10; //13 for paper #1, 17 for paper #2
    axes.AddAxis(BinAxis("mc_neutrino_energy", e_min, e_max, n_bins));

    double param_s13 = 0.02303;

    // determine the oscillated normalisation (the fit returns normalisation of unoscillated spectra)
    printf("bestFit[d21]:%.4e bestFit[s12]:%.4f\n", best_fit["d21"], best_fit["s12"]);
    //double normalisation_constants[n_pdf];
    //double normalisation_osc[n_pdf];
    std::vector<double> normalisation_osc = lh_function.GetWorkingNormalisations();
    printf("size:%d\tnpdf:%d\n", name, normalisation_osc.size(), n_pdf);
    for (ULong64_t i = 0; i < n_pdf; i++){
        sprintf(name, "%s_norm", reactor_names[i].c_str());
        //normalisation_constants[i] = lh_function.GetWorkingNormalisations().at(i) / reactor_scale[i];
        //normalisation_osc = lh_function.GetWorkingNormalisations();//normalisation_constants[i] * reactor_scale[i];
        printf("%s:\tscale:%.4f\tnorm_osc: %.3f\tratio: %.4f\n", name, reactor_scale[i], normalisation_osc.at(i), reactor_scale[i]);
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

        //sprintf(name,"%s_norm",reactor_names[i].c_str()); //don't use the value directly from the fit, use the oscillated normalisation
        reactor_pdf_fitosc[i]->Scale(best_fit.at(name));
        //reactor_pdf_fitosc[i]->Scale(reactor_scale[i]);

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
