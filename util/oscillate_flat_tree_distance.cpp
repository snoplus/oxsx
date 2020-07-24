#include <TFile.h>
#include <TMath.h>
#include <string>
#include <TNtuple.h>
#include <iostream>
#include <TObject.h>
#include <TRandom3.h>
#include "oscillate_util.cpp"

Int_t main(int argc, char* argv[]){
  if (argc != 7)
    std::cout<<"6 arguments expected!"<<std::endl;
  else{
    const char* root_in = argv[1];
    const char* root_out = argv[2];

    double del_m_sqr_21 = atof(argv[3]);
    double sin_sqr_theta_12 = atof(argv[4]);
    double sin_sqr_theta_13 = atof(argv[5]);
    const double distance = atof(argv[6]);
    write_file(root_in, root_out, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13, distance);
  }
}
