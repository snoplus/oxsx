#include <TFile.h>
#include <TMath.h>
#include <string>
#include <TNtuple.h>
#include <iostream>
#include <TObject.h>
#include <TRandom3.h>
#include "oscillate_util.cpp"

Int_t main(int argc, char* argv[]){
  if (argc != 10)
    std::cout<<"9 arguments expected!"<<std::endl;
  else{
    const std::string root_in = argv[1];
    const std::string root_out = argv[2];
    const size_t pass_min = atoi(argv[6]);
    const size_t pass_max = atoi(argv[7]);
    const std::string cut_label = argv[8];
    const std::string osc_label = argv[9];

    char root_in_i[1000];
    char root_out_i[1000];

    for (size_t pass_i = pass_min; pass_i <= pass_max; pass_i++){
        printf("oscillating: %s_pass%d_%s.root\n", root_in.c_str(), pass_i, cut_label.c_str());
        sprintf(root_in_i, "%s_pass%d_%s.root", root_in.c_str(), pass_i, cut_label.c_str());
        sprintf(root_out_i, "%s_pass%d_%s_%s.root", root_out.c_str(), pass_i, cut_label.c_str(), osc_label.c_str());

        double del_m_sqr_21 = atof(argv[3]);
        double sin_sqr_theta_12 = atof(argv[4]);
        double sin_sqr_theta_13 = atof(argv[5]);
        write_file(root_in_i, root_out_i, del_m_sqr_21, sin_sqr_theta_12, sin_sqr_theta_13);
    }
  }
}
