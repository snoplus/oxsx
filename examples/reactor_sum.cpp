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
//#include "AntinuUtilsKamland.cpp"
#include "AntinuUtils.cpp"

int main(int argc, char *argv[]){

    if (argc != 6){
        printf("Error: 5 arguments expected.\n");
        return 1; // return>0 indicates error code
    }
    else{
	// get command line arguments
        const std::string &in_path = argv[1];
        const std::string &data_path = argv[2];
        const std::string &info_file = argv[3];
        const size_t flux_data = atoi(argv[4]);
        const std::string &out_filename = argv[5];

        printf("--------------------------------------\n");

        // read in reactor information
        std::vector<std::string> reactor_names;
        std::vector<Double_t> distances;
        std::vector<std::string> reactor_types;
        std::vector<ULong64_t> n_cores;
        std::vector<Double_t> powers;
        std::vector<Double_t> power_errs;
        readInfoFile(info_file, reactor_names, distances, reactor_types, n_cores, powers, power_errs);

	    // initialise variables for pdf's
        const ULong64_t n_pdf = reactor_names.size();
        BinnedED **spectra_pdf = new BinnedED*[n_pdf]; // PWR=0, PHWR=1
        Double_t *reactor_scale = new Double_t[n_pdf];
        Double_t *reactor_scale_err = new Double_t[n_pdf];

        for (size_t i=0; i<(size_t)reactor_names.size(); i++) // use uncertainty on power for uncertainty onflux for now
            reactor_scale_err[i] = power_errs[i]/powers[i];

	    // initialise variables for fitting
	    double param_d21 = 7.58e-5;//kamland#1=6.9e-5;//us=7.58e-5;
	    double param_s12 = 0.359;
	    double param_s13 = 0.02303;
	    BinnedNLLH lh_function;
	    SurvProb **surv_prob = new SurvProb*[n_pdf];
	    NuOsc **reactor_systematic = new NuOsc*[n_pdf];
	    BinnedED **reactor_pdf = new BinnedED*[n_pdf];
	    bool fit_validity = false;
	    ParameterDict best_fit;

        // intialise pdf's
        BinnedED data_set_pdf = LHFit_initialise(spectra_pdf, reactor_scale, in_path, data_path, reactor_names, flux_data);

        // fit pdf's
	    LHFit_fit(data_set_pdf, spectra_pdf, lh_function, surv_prob, reactor_systematic, reactor_pdf, reactor_scale, reactor_scale_err, reactor_names, distances, reactor_types, param_d21, param_s12,  param_s13, fit_validity, best_fit);

        // produce histograms
        LHFit_produce_histograms(data_set_pdf, spectra_pdf, lh_function, surv_prob, reactor_pdf, reactor_scale, reactor_names, distances, best_fit, out_filename);

        printf("--------------------------------------\n");
        return 0;
    }
}
