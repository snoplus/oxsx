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

void readInfoFile(const std::string &runInfoFileName, std::vector<std::string> &reactor_names, std::vector<Double_t> &distances, std::vector<std::string> &reactor_types, std::vector<ULong64_t> &n_cores, std::vector<Double_t> &powers ) {
    // Read couchDB run-info text file
    std::ifstream in;
    in.open(runInfoFileName.c_str());
    //std::cout << "opening file: " << runInfoFileName.c_str() << std::endl;

    std::fill(reactor_names.begin(), reactor_names.end(), "");
    std::fill(distances.begin(), distances.end(), 0.);
    std::fill(reactor_types.begin(), reactor_types.end(), "");
    std::fill(n_cores.begin(), n_cores.end(), 0);
    std::fill(powers.begin(), powers.end(), 0.);

    std::string reactor_name,distance,reactor_type,n_core,power;
    ULong64_t line_no = 0;

    // read until end of file.
    while(in.good()){
        std::getline(in,reactor_name,',');
        std::getline(in,distance,',');
        std::getline(in,reactor_type,',');
        std::getline(in,n_core,',');
        std::getline(in,power,'\n');

        if (line_no>0){ //skip csv header
            if (strcmp(reactor_name.c_str(),"")!=0) {
                reactor_names.push_back(reactor_name);
                distances.push_back(atof(distance.c_str()));
                reactor_types.push_back(reactor_type.c_str());
                n_cores.push_back(atoi(n_core.c_str()));
                powers.push_back(atof(power.c_str()));

                //std::cout << "v: reactor_name: " << reactor_names[line_no-1] << ", distance: " << distances[line_no-1] << ", reactor_type: " << reactor_types[line_no-1] << ", n_core: " << n_cores[line_no-1] << ", power: " << powers[line_no-1] << std::endl; //debug check ('-1' for header)
            }
        }
        line_no++;
    }
    in.close();
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
}

BinnedED LHFit_initialise(BinnedED **spectra_pdf, Double_t *reactor_scale, const std::string &in_path, const std::string &data_path, std::vector<std::string> &reactor_names, std::vector<Double_t> &distances, std::vector<std::string> &reactor_types){
    
    printf("Begin init--------------------------------------\n");
    printf("LHFit_initialise...\n");

    char name[1000];
    ULong64_t flux = 100;

    // setup spectra filenames
    const std::string pwr_file = std::string("combinedpwr");
    const std::string phwr_file = std::string("combinedphwr");
    printf("Reading individual reactor spectra from: %s\n", in_path.c_str());
    sprintf(name, "%sflux%llu/[reactor_name]_flux%llu_day360_cleanround1_oxsx.root", in_path.c_str(), flux, flux);
    printf("Using unoscillated spectrum for reactors:\n\t%s\n",name);
    printf("Using unoscillated spectrum for:\n");
    sprintf(name, "%sflux%llu/%s_flux%llu_day360_cleanround1_oxsx.root", in_path.c_str(), flux, pwr_file.c_str(), flux);
    printf("\ttype PWR & BWR: %s\n", name);
    const std::string pwr_unosc_path = std::string(name);
    sprintf(name, "%sflux%llu/%s_flux%llu_day360_cleanround1_oxsx.root", in_path.c_str(), flux, phwr_file.c_str(), flux);
    printf("\ttype PHWR: %s\n", name);
    const std::string phwr_unosc_path = std::string(name);

    // setup (oscillated) data filename
    sprintf(name, "%s", data_path.c_str());
    printf("Loading data spectrum: %s\n", name);

    // setup ntuple
    ObsSet data_rep(0);

    // set up binning
    AxisCollection axes;
    Double_t e_min = 2;
    Double_t e_max = 8;
    Int_t n_bins = (8-2)*10;
    axes.AddAxis(BinAxis("mc_neutrino_energy", e_min, e_max, n_bins));

    // load (oscillated) data ntuple
    BinnedED data_set_pdf("dataSetPdf", axes);
    data_set_pdf.SetObservables(data_rep);
    ROOTNtuple data_ntp(data_path, "nt");
    for(ULong64_t i = 0; i < data_ntp.GetNEntries(); i++)
        data_set_pdf.Fill(data_ntp.GetEntry(i));
    Double_t data_set_integral = data_set_pdf.Integral(); //record the integral before normalising
    data_set_pdf.Normalise();

    // number of reactors
    const ULong64_t n_pdf = reactor_names.size();

    // load oscillated reactor ntuples (for scale factors)
    Double_t reactor_integrals[n_pdf];
    for(ULong64_t i = 0; i < n_pdf; i++){
        sprintf(name, "%sflux%llu/osc2/%s_flux%llu_day360_cleanround1_osc2_oxsx.root", in_path.c_str(), flux, reactor_names[i].c_str(), flux);
        ROOTNtuple reactor_ntp(name, "nt");
        sprintf(name, "%s_pdf", reactor_names[i].c_str());
        BinnedED *reactor_pdf = new BinnedED(name, axes);
        reactor_pdf->SetObservables(0);
        for(ULong64_t j = 0; j < reactor_ntp.GetNEntries(); j++)
            reactor_pdf->Fill(reactor_ntp.GetEntry(j));
        reactor_integrals[i] = reactor_pdf->Integral(); //record the integral
        reactor_scale[i] = data_set_integral/reactor_integrals[i];
	printf("Loading reactor: %s integral:%.0f/data integral:%.0f => Ratio:%4f (scale:%4f)\n", reactor_names[i].c_str(), reactor_integrals[i], data_set_integral, 1./reactor_scale[i], reactor_scale[i]);
    }

    // load unoscillated spectra
    // pwr spectrum
    ROOTNtuple pwr_spectrum_ntp(pwr_unosc_path, "nt");
    sprintf(name, "pwr_spectrum_pdf");
    spectra_pdf[0] = new BinnedED(name, axes);
    spectra_pdf[0]->SetObservables(0);
    for(ULong64_t i = 0; i < pwr_spectrum_ntp.GetNEntries(); i++)
        spectra_pdf[0]->Fill(pwr_spectrum_ntp.GetEntry(i));
    Double_t pwr_spectrum_integral = spectra_pdf[0]->Integral(); //record the integral before normalising
    spectra_pdf[0]->Normalise();

    // phwr spectrum
    ROOTNtuple phwr_spectrum_ntp(phwr_unosc_path, "nt");
    sprintf(name, "phwr_spectrum_pdf");
    spectra_pdf[1] = new BinnedED(name, axes);
    spectra_pdf[1]->SetObservables(0);
    for(ULong64_t i = 0; i < phwr_spectrum_ntp.GetNEntries(); i++)
        spectra_pdf[1]->Fill(phwr_spectrum_ntp.GetEntry(i));
    Double_t phwr_spectrum_integral = spectra_pdf[1]->Integral(); //record the integral before normalising
    spectra_pdf[1]->Normalise();

    printf("End init--------------------------------------\n");
    return data_set_pdf;
}

Double_t LHFit_fit(BinnedED &data_set_pdf, BinnedED **spectra_pdf, Double_t *reactor_scale, std::vector<std::string> &reactor_names, std::vector<Double_t> &distances, std::vector<std::string> &reactor_types, Double_t param_d21, Double_t param_s12, Double_t param_s13){

    printf("Begin fit--------------------------------------\n");
    printf("LHFit_fit:: del_21:%.5f, sin2_12:%.5f, sin2_13:%.5f\n", param_d21, param_s12, param_s13);

    char name[1000];
    const ULong64_t n_pdf = reactor_names.size();
    ObsSet data_rep(0);
    // set up binning
    AxisCollection axes;
    Double_t e_min = 2;
    Double_t e_max = 8;
    Int_t n_bins = (8-2)*10;
    axes.AddAxis(BinAxis("mc_neutrino_energy", e_min, e_max, n_bins));

    // create LH function
    BinnedNLLH lh_function;
    lh_function.SetBufferAsOverflow(true);
    int buff = 5;
    lh_function.SetBuffer(0, buff, buff);
    lh_function.SetDataDist(data_set_pdf); // initialise withe the data set
    
    // setup max and min ranges
    ParameterDict minima;
    ParameterDict maxima;
    ParameterDict initial_val;
    ParameterDict initial_err;

    // setup survival probability
    printf("Setup survival probability...\n");
    SurvProb *surv_prob[n_pdf];
    //NuOsc *reactor_systematic[n_pdf];
    BinnedED **reactor_pdf = new BinnedED*[n_pdf];
    for (ULong64_t i = 0; i < n_pdf; i++){
        // for each reactor, load spectrum pdf for reactor type
        reactor_pdf[i] = new BinnedED(reactor_names[i], axes);
        reactor_pdf[i]->SetObservables(0);

        // setup survival probability
        sprintf(name, "%s_survival", reactor_names[i].c_str());
        surv_prob[i] = new SurvProb(param_d21, param_s12, distances[i]);
        surv_prob[i]->Setsinsqrtheta13s(param_s13); // manual, fixed setting of theta13
        surv_prob[i]->RenameParameter("delmsqr21_0", "d21"); // rename all parameters to the same for the fit
        surv_prob[i]->RenameParameter("sinsqrtheta12_0", "s12");
        sprintf(name, "%s_systematic", reactor_names[i].c_str());
        NuOsc reactor_systematic(name);
        reactor_systematic.SetFunction(surv_prob[i]);
        reactor_systematic.SetAxes(axes);
        reactor_systematic.SetTransformationObs(data_rep);
        reactor_systematic.SetDistributionObs(data_rep);
        reactor_systematic.Construct();

        if ((reactor_types[i]=="PWR")||(reactor_types[i]=="BWR"))
            reactor_pdf[i]->Add(reactor_systematic(*spectra_pdf[0]),1); //use PWR pdf
        else if (reactor_types[i]=="PHWR")
                reactor_pdf[i]->Add(reactor_systematic(*spectra_pdf[1]),1); //use PHWR pdf
            else{
               printf("Throw: Reactor doesn't match any loaded type...\n");
               continue;
            }
        reactor_pdf[i]->Normalise();
        reactor_pdf[i]->Scale(1./reactor_scale[i]); // normalise to integral for each reactor

        // Setting optimisation limits
        sprintf(name, "%s_norm", reactor_names[i].c_str());
        Double_t min = 0.99/reactor_scale[i];
        Double_t max = 1.01/reactor_scale[i];
        if (min < 0.0001) min = 0.0001;
        if (max > 0.5) max = 0.5; // no reactor is more than 50% of the signal
        //if (max > 1.0) max = 1.0;
        minima[name] = min;
        maxima[name] = max;
        printf("  added reactor %d/%d: %s, min:%.3f max:%.3f\n", i+1, n_pdf, reactor_names[i].c_str(), min, max);
        initial_val[name] = 1./reactor_scale[i];
        initial_err[name] = 0.001*initial_val[name];

        lh_function.AddDist(*reactor_pdf[i]);

        //sprintf(name, "%s_norm", reactor_names[i].c_str());
        //lh_function.SetConstraint(name,7591,380); // use scale information here too?
    }

    //lh_function.SetConstraint("d21", 7e-5, 9e-5); // constraints?
    //lh_function.SetConstraint("s12", 0.5, 0.3);

    // fit
    printf("Built LH function, fitting...\n");
    Minuit min;
    min.SetMethod("Migrad");
    min.SetMaxCalls(10000000);
    min.SetMinima(minima);
    min.SetMaxima(maxima);
    min.SetInitialValues(initial_val);
    min.SetInitialErrors(initial_err);

    FitResult fit_result = min.Optimise(&lh_function);
    ParameterDict best_fit = fit_result.GetBestFit();
    fit_result.Print();

    Double_t lh_val =(-1)*lh_function.Evaluate();

    printf("lh_value:%.5f\n", lh_val);
    printf("End fit--------------------------------------\n");
    return lh_val;
}

int main(int argc, char *argv[]) {

    if (argc != 6){
        std::cout<<"Error: 5 arguments expected."<<std::endl;
        return 1; // return>0 indicates error code
    }
    else{
        const std::string &in_path = argv[1];
        const std::string &data_path = argv[2];
        const std::string &infoFile = argv[3];
        const std::string &parameterFile = argv[4];
        const std::string &outFile = argv[5];
        printf("Begin--------------------------------------\n");

        // read in reactor information
        std::vector<std::string> reactor_names;
        std::vector<Double_t> distances;
        std::vector<std::string> reactor_types;
        std::vector<ULong64_t> n_cores;
        std::vector<Double_t> powers;
        readInfoFile(infoFile, reactor_names, distances, reactor_types, n_cores, powers);

        // read in parameter information
        std::vector<Double_t> d_21s;
        std::vector<Double_t> s_12s;
        std::vector<Double_t> s_13s;
        readParameterFile(parameterFile, d_21s, s_12s, s_13s);

        // // print out read info
        // for (size_t i=0; i<(size_t)reactor_names.size(); i++)
            // printf("i:%llu,reactor_names[i]:%s, distance: %.5f, type: %s, n_cores: %llu, power: %.5f \n", i, reactor_names[i].c_str(), distances[i], reactor_types[i].c_str(), n_cores[i], powers[i]);

        // // print out read parameters
        // for (size_t i=0; i<(size_t)d_21s.size(); i++)
            // printf("i:%llu, d_21:%.5f, s_12:%.5f, s_13:%.5f\n", i, d_21s[i], s_12s[i], s_13s[i]);

        const ULong64_t n_pdf = reactor_names.size();
        BinnedED **spectra_pdf = new BinnedED*[n_pdf]; // PWR=0, PHWR=1
        Double_t *reactor_scale = new Double_t[n_pdf];
        const ULong64_t n_parameter_sets = d_21s.size();
        Double_t lh_values[n_parameter_sets];
        
        BinnedED data_set_pdf = LHFit_initialise(spectra_pdf, reactor_scale, in_path, data_path, reactor_names, distances, reactor_types);

        for (ULong64_t i=0; i<n_parameter_sets; i++) {
            lh_values[i] = LHFit_fit(data_set_pdf, spectra_pdf, reactor_scale, reactor_names, distances, reactor_types, d_21s[i], s_12s[i], s_13s[i]);
        }

        // write output to file
        ULong64_t fitValidity = 1; //make this do something useful.. (i.e. when minuit throws out an error this should be 0.)

        //Write fit coefficients to txt file
        printf("writing to: %s\n", outFile.c_str());
        FILE *fOut = fopen(outFile.c_str(), "w");
        fprintf(fOut,"d21,s12,s13,lh_value,fitValidity\n");
        for (ULong64_t i=0; i<n_parameter_sets; i++)
            fprintf(fOut,"%.9f,%.7f,%.7f,%.5f,%llu\n", d_21s[i], s_12s[i], s_13s[i], lh_values[i], fitValidity);
        fclose(fOut);

        printf("End--------------------------------------\n");
        return 0; // completed successfully
    }
}
