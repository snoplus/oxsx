#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <string>
#include <TTree.h>
#include <TVector3.h>
#include <iostream>
#include <TObject.h>
#include <math.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TNtuple.h>
#include "antinu_cleaning_util.cpp"

Int_t main(Int_t argc, char *argv[]) {

    if (argc != 15) {
        std::cout<<"Error: 14 arguments expected. Got: "<<argc-1<<std::endl;
        return 1; // return>0 indicates error code
    }
    else {
        const std::string &filename_input = argv[1];
        const std::string &filename_output = argv[2];
        double energy_ep_min = atof(argv[3]);
        double energy_ep_max = atof(argv[4]);
        double energy_n_min = atof(argv[5]);
        double energy_n_max = atof(argv[6]);
        double deltaTmin = atof(argv[7]);
        double deltaTmax = atof(argv[8]);
        double promptRmax = atof(argv[9]);
        double lateRmax = atof(argv[10]);
        double deltaRmax = atof(argv[11]);
        const size_t pass_min = atoi(argv[12]);
        const size_t pass_max = atoi(argv[13]);
        const size_t cut_no = atoi(argv[14]);

        char filename_input_i[1000];
        char filename_output_i[1000];

        for (size_t pass_i = pass_min; pass_i <= pass_max; pass_i++){
            printf("cleaning: %s_pass%d.root\n", filename_input.c_str(), pass_i);
            sprintf(filename_input_i, "%s_pass%d.ntuple.root", filename_input.c_str(), pass_i);
            sprintf(filename_output_i, "%s_pass%d_cleanround%d.ntuple.root", filename_output.c_str(), pass_i, cut_no);

            //write csv output file to show process has begun (values filled upon completion)
            char *name = new char[1000];
            sprintf(name, "%s.csv",filename_output.c_str());
            FILE *fOut = fopen(name,"w");
            fprintf(fOut,"filename_input,filename_output,energy_ep_min,energy_ep_max,energy_n_min,energy_n_max,deltaTmin,deltaTmax,promptRmax,lateRmax,deltaRmax,initial_entries,final_entries,finished\n");
            fprintf(fOut,"%s,%s,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i\n", filename_input_i, filename_output_i, -9000, -9000, -9000, -9000, -9000, -9000, -9000 -9000, -9000, -9000, -9000, -9000,0);
            fclose(fOut);

            process_cuts(filename_input_i, filename_output_i, energy_ep_min, energy_ep_max, energy_n_min, energy_n_max, deltaTmin, deltaTmax, promptRmax, lateRmax, deltaRmax);
        }

        return 0; // completed successfully
    }
}
