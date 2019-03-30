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
    Double_t e_min = 2;//0.425*6;//0.425*2; //2; //*6 for kamland paper #1, *2 for paper #2
    Double_t e_max = 8;//0.425*19;//0.425*19; //8;
    Int_t n_bins = 60;//13;//17; //(8-2)*10; //13 for paper #1, 17 for paper #2
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
