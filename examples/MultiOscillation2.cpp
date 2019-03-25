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
#include "AntinuUtils.cpp"

int main(int argc, char *argv[]) {

    if (argc != 7){
        std::cout<<"Error: 6 arguments expected."<<std::endl;
        return 1; // return>0 indicates error code
    }
    else{
        // get command line arguments
        const std::string &in_path = argv[1];
        const std::string &data_path = argv[2];
        const std::string &info_file = argv[3];
        const std::string &parameter_file = argv[4];
        const size_t flux_data = atoi(argv[5]);
        const std::string &out_file = argv[6];
        printf("Begin--------------------------------------\n");

        // read in reactor information
        std::vector<std::string> reactor_names;
        std::vector<Double_t> distances;
        std::vector<std::string> reactor_types;
        std::vector<ULong64_t> n_cores;
        std::vector<Double_t> powers;
        std::vector<Double_t> power_errs;
        readInfoFile(info_file, reactor_names, distances, reactor_types, n_cores, powers, power_errs);

        // read in parameter information
        std::vector<Double_t> d_21s;
        std::vector<Double_t> s_12s;
        std::vector<Double_t> s_13s;
        readParameterFile(parameter_file, d_21s, s_12s, s_13s);

        // initialise variables for pdf's
        const ULong64_t n_pdf = reactor_names.size();
        BinnedED **spectra_pdf = new BinnedED*[n_pdf]; // PWR=0, PHWR=1
        Double_t *reactor_scale = new Double_t[n_pdf];
        Double_t *reactor_scale_err = new Double_t[n_pdf];

        for (size_t i=0; i<(size_t)reactor_names.size(); i++) // use uncertainty on power for uncertainty onflux for now
            reactor_scale_err[i] = power_errs[i]/powers[i];

        // intialise pdf's
        BinnedED data_set_pdf = LHFit_initialise(spectra_pdf, reactor_scale, in_path, data_path, reactor_names, flux_data);

        // fit pdf's
        bool fit_validity;
        const ULong64_t n_parameter_sets = d_21s.size();
        double lh_values[n_parameter_sets];
        for (ULong64_t i=0; i<n_parameter_sets; i++) {

            // initialise variables for fitting
            SurvProb **surv_prob = new SurvProb*[n_pdf];
            NuOsc **reactor_systematic = new NuOsc*[n_pdf];
            BinnedED **reactor_pdf = new BinnedED*[n_pdf];
            fit_validity = false; //reset
            ParameterDict best_fit;
            BinnedNLLH lh_function;

            // fit
            LHFit_fit(data_set_pdf, spectra_pdf, lh_function, surv_prob, reactor_systematic, reactor_pdf, reactor_scale, reactor_scale_err, reactor_names, distances, reactor_types, d_21s[i], s_12s[i], s_13s[i], fit_validity, best_fit);

            // evaluate LH (this changes parameters - must get parameters for oscillation norm before doing this)
            if (fit_validity == true)
                lh_values[i] =(-1)*lh_function.Evaluate();
            else
                lh_values[i] = 9999; // positive non-sensical value to return if fit is not valid
            printf("fit_validity: %d\td_21^2: %.9f\ts_12^2: %.5f\ts_13^2: %.5f\tlh_value: %.4f\n", fit_validity, d_21s[i], s_12s[i], s_13s[i], lh_values[i]);
        }


        //Write fit coefficients to txt file
        printf("writing to: %s\n", out_file.c_str());
        FILE *fOut = fopen(out_file.c_str(), "w");
        fprintf(fOut,"d21,s12,s13,lh_value,fitValidity\n");
        for (ULong64_t i=0; i<n_parameter_sets; i++)
            fprintf(fOut,"%.9f,%.7f,%.7f,%.5f,%llu\n", d_21s[i], s_12s[i], s_13s[i], lh_values[i], fit_validity);
        fclose(fOut);

        printf("End--------------------------------------\n");
        return 0; // completed successfully
    }
}

