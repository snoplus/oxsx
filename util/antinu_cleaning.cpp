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

    if (argc != 14) {
        std::cout<<"Error: 13 arguments expected. Got: "<<argc-1<<std::endl;
        return 1; // return>0 indicates error code
    }
    else {
        const std::string &filename_input = argv[1];
        const std::string &filename_output = argv[2];
        double energy_ep_min = atof(argv[3]);
        double energy_ep_max = atof(argv[4]);
        double energy_n_min = atof(argv[5]);
        double energy_n_max = atof(argv[6]);
        double Rmax = atof(argv[7]);
        double Zmin = atof(argv[8]);
        double Zmax = atof(argv[9]);
        double deltaTmin = atof(argv[10]);
        double deltaTmax = atof(argv[11]);
        double deltaRmax = atof(argv[12]);
        int alphan_low_high = atoi(argv[13]);

        //write csv output file to show process has begun (values filled upon completion)
        char *name = new char[1000];
        sprintf(name, "%s.csv",filename_output.c_str());
        std::cout<<name<<std::endl;
        FILE *fOut = fopen(name,"w");
        fprintf(fOut,"filename_input,filename_output,energy_ep_min,energy_ep_max,energy_n_min,energy_n_max,Rmax,Zmin,Zmax,deltaTmin,deltaTmax,deltaRmax,initial_entries,final_entries,finished\n");
        fprintf(fOut,"%s,%s,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i\n", filename_input.c_str(), filename_output.c_str(), -9000, -9000, -9000, -9000, -9000, -9000, -9000 -9000, -9000, -9000, -9000, -9000,0);
        fclose(fOut);
        process_cuts(filename_input, filename_output, Rmax, Zmax, Zmin, energy_ep_min, energy_ep_max, energy_n_min, energy_n_max, deltaTmin, deltaTmax, deltaRmax, alphan_low_high);

        return 0; // completed successfully
    }
}
