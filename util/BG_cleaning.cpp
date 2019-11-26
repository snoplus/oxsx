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
#include "BG_cleaning_util.cpp"

Int_t main(Int_t argc, char *argv[]) {

    if (argc != 6) {
        std::cout<<"Error: 5 arguments expected. Got: "<<argc-1<<std::endl;
        return 1; // return>0 indicates error code
    }
    else {
        const std::string &filename_input = argv[1];
        const std::string &filename_output = argv[2];
        double energy_min = atof(argv[3]);
        double energy_max = atof(argv[4]);
        double Rmax = atof(argv[5]);

        //write csv output file to show process has begun (values filled upon completion)
        char *name = new char[1000];
        sprintf(name, "%s.csv",filename_output.c_str());
	std::cout<<name<<std::endl;
	FILE *fOut = fopen(name,"w");
	fprintf(fOut,"filename_input,filename_output,energy_min,energy_max,Rmax,initial_entries,final_entries,finished\n");
	fprintf(fOut,"%s,%s,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i\n", filename_input.c_str(), filename_output.c_str(), -9000, -9000, -9000, -9000, 0);
	fclose(fOut);
	process_cuts(filename_input, filename_output, energy_min, energy_max, Rmax);

        return 0; // completed successfully
    }
}
