#include <TFile.h>
#include <TMath.h>
#include <string>
#include <TNtuple.h>
#include <iostream>
#include <TObject.h>
#include <TRandom3.h>
#include "oscillate_util.cpp"

Int_t main(int argc, char* argv[]){
  if (argc != 4)
    std::cout<<"3 arguments expected!"<<std::endl;
  else{
    const char* root_in = argv[1];
    const char* root_out = argv[2];

    double sin_sqr_theta_12 = atof(argv[3]);
    write_file_geo(root_in, root_out, sin_sqr_theta_12);
  }
}
