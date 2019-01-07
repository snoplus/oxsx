#include <TFile.h>
#include <RAT/DB.hh>
#include <TH1D.h>
#include <TMath.h>
#include <string>
#include <TNtuple.h>
#include <iostream>
#include <RAT/DU/Utility.hh>
#include <TObject.h>

void round_1_cuts(const std::string filename_input, const std::string filename_output, double nhitMinPrompt, double nhitMinLate, double deltaT = 1000000, double deltaP = 6000, const bool use_mc = false) {
    // checks for event pairs with time difference < 1000micros and both events within 6000mm
    TFile *file_input = new TFile(filename_input.c_str());
    TTree *tree_input = (TTree*)file_input->Get("tt");
    int n_entries = tree_input->GetEntries();
    Int_t days, nextDays;
    Int_t secs, nextSec;
    Int_t nsecs, nextNSec;
    Int_t nhit, nextnhit;
    double timeD;
    float position_r;
    bool flag1 = true; // assume that the first event is invalidated by an event we didn't catch
    bool flag2, isCoincidence, time_window_pass, position_r_pass;
 
    std::vector<float> mc_positron_position_r;
    //std::vector<float> mc_neutron_position_r;
    std::vector<bool> ev_fit_validity;
    std::vector<int> ev_nhit;
    std::vector<float> ev_fit_position_r;
    std::vector<unsigned int> ev_time_days;
    std::vector<unsigned int> ev_time_seconds;
    std::vector<float> ev_time_nanoseconds;
    
    if (use_mc){        
        tree_input->SetBranchAddress("mc_positron_position_r",&mc_positron_position_r);
        //tree_input->SetBranchAddress("mc_neutron_position_r",&mc_neutron_position_r);
    }
    else{
        tree_input->SetBranchAddress("ev_fit_position_r",&ev_fit_position_r);
    }
    
    tree_input->SetBranchAddress("ev_fit_validity",&ev_fit_validity);
    tree_input->SetBranchAddress("ev_nhit",&ev_nhit);
    tree_input->SetBranchAddress("ev_time_days",&ev_time_days);
    tree_input->SetBranchAddress("ev_time_seconds",&ev_time_seconds);
    tree_input->SetBranchAddress("ev_time_nanoseconds",&ev_time_nanoseconds);

    TFile *file_output = new TFile(filename_output.c_str(),"RECREATE");
    TTree *tree_output = tree_input->CloneTree(0);

    for (size_t i=0; i < (n_entries-1); i++){
        tree_input->GetEntry(i);
        nhit = 0;
        nextnhit = 0;
        days = 0; //reset
        nextDays = 0;
        secs = 0;
        nextSec = 0;
        nsecs = 0;
        nextNSec = 0;
        timeD = 0.;
        position_r = 0.;
        isCoincidence = false;
        flag2 = false;
        time_window_pass = false;
        position_r_pass = false;

        // get data from vectors
        days = ev_time_days[i];
        nextDays = ev_time_days[i+1];
        secs = ev_time_seconds[i];
        nextSec = ev_time_seconds[i+1];
        nsecs = ev_time_nanoseconds[i];
        nextNSec = ev_time_nanoseconds[i+1];
        
        // specific for mc data
        if (use_mc){
            position_r = mc_positron_position_r[i];
        }
        else{ //specific for data
            position_r = ev_fit_position_r[i];
        }
        
        if(std::abs(nextDays - days) > 0) //begin time tests
            isCoincidence = false;
        else { // days are equal
            if(std::abs(nextSec - secs) > 0)
                isCoincidence = false;
            else { // seconds are equal
                if(std::abs(nextNSec - nsecs) < deltaT){ // nanosec are to within deltaT
                    isCoincidence = true;
                    timeD = std::abs(nextNSec - nsecs);
                }
                else
                    isCoincidence = false; // so close..
            }
        }
        if (isCoincidence){
            if (nextNSec>=nsecs)
                flag2 = (isCoincidence && (nhit >= nhitMinPrompt) && (nextnhit >= nhitMinLate));
            else
                flag2 = (isCoincidence && (nhit >= nhitMinLate) && (nextnhit >= nhitMinPrompt));
        }

        // conditions for event to pass cuts
        if ((position_r < deltaP)&&ev_fit_validity[i]){
            position_r_pass = true;
            time_window_pass = (flag1 || flag2);
        }

        if(time_window_pass && position_r_pass && ev_fit_validity[i])
            tree_output->Fill();

        // shuffle the flag along, if this event is tagged by it's coincidence with the next (flag2)
        // the next one is flagged by its coincidence with the one before (flag1)
        flag1 = flag2;
    }

    tree_output->Print();
    tree_output->AutoSave();
    file_output->Close();
}

void round_2_cuts(const std::string filename_input, const std::string filename_output, double nhitMinPrompt, double nhitMinLate, double deltaT, double deltaP, const bool use_mc) {
// checks for event pairs with time difference < 275micros, both events within 5000mm and position difference between events (in pair) < 1250mm
    TFile *file_input = new TFile(filename_input.c_str());
    TTree *tree_input = (TTree*)file_input->Get("tt");
    int n_entries = tree_input->GetEntries();
 
    if (use_mc){
        std::vector<float> mc_positron_position_r;
        std::vector<float> mc_positron_position_x;
        std::vector<float> mc_positron_position_y;
        std::vector<float> mc_positron_position_z;
        std::vector<float> mc_neutron_position_r;
        std::vector<float> mc_neutron_position_x;
        std::vector<float> mc_neutron_position_y;
        std::vector<float> mc_neutron_position_z;
        std::vector<unsigned int> mc_positron_time_days; // this requires MC tracking information?
        std::vector<unsigned int> mc_positron_time_seconds;
        std::vector<float> mc_positron_time_nanoseconds;
        std::vector<unsigned int> mc_neutron_time_days;
        std::vector<unsigned int> mc_neutron_time_seconds;
        std::vector<float> mc_neutron_time_nanoseconds;
        tree_input->SetBranchAddress("mc_positron_position_r",&mc_positron_position_r);
        tree_input->SetBranchAddress("mc_positron_position_x",&mc_positron_position_x);
        tree_input->SetBranchAddress("mc_positron_position_y",&mc_positron_position_y);
        tree_input->SetBranchAddress("mc_positron_position_z",&mc_positron_position_z);
        tree_input->SetBranchAddress("mc_neutron_position_r",&mc_neutron_position_r);
        tree_input->SetBranchAddress("mc_neutron_position_x",&mc_neutron_position_x);
        tree_input->SetBranchAddress("mc_neutron_position_y",&mc_neutron_position_y);
        tree_input->SetBranchAddress("mc_neutron_position_z",&mc_neutron_position_z);
        tree_input->SetBranchAddress("ev_time_days",&mc_positron_time_days);
        tree_input->SetBranchAddress("ev_time_seconds",&mc_positron_time_seconds);
        tree_input->SetBranchAddress("ev_time_nanoseconds",&mc_positron_time_nanoseconds);
        tree_input->SetBranchAddress("ev_time_days",&mc_neutron_time_days);
        tree_input->SetBranchAddress("ev_time_seconds",&mc_neutron_time_seconds);
        tree_input->SetBranchAddress("ev_time_nanoseconds",&mc_neutron_time_nanoseconds);
    }
    else{
        std::vector<bool> ev_fit_validity;
        std::vector<float> ev_fit_position_r;
        std::vector<float> ev_fit_position_x;
        std::vector<float> ev_fit_position_y;
        std::vector<float> ev_fit_position_z;
        std::vector<unsigned int> ev_time_days;
        std::vector<unsigned int> ev_time_seconds;
        std::vector<float> ev_time_nanoseconds;
        tree_input->SetBranchAddress("ev_fit_validity",&ev_fit_validity);
        tree_input->SetBranchAddress("ev_fit_position_r",&ev_fit_position_r);
        tree_input->SetBranchAddress("ev_fit_position_x",&ev_fit_position_x);
        tree_input->SetBranchAddress("ev_fit_position_y",&ev_fit_position_y);
        tree_input->SetBranchAddress("ev_fit_position_z",&ev_fit_position_z);
        tree_input->SetBranchAddress("ev_time_days",&ev_time_days);
        tree_input->SetBranchAddress("ev_time_seconds",&ev_time_seconds);
        tree_input->SetBranchAddress("ev_time_nanoseconds",&ev_time_nanoseconds);
    }

    TFile *file_output = new TFile(filename_output.c_str(),"RECREATE");
    TTree *tree_output = tree_input->CloneTree(0);

    for (int i = 0; i < n_entries; i++){
        tree_input->GetEntry(i);
        if (use_mc){
            tree_output->Fill();
        }
        else{
            tree_output->Fill();
        }
    }

    tree_output->Print();
    tree_output->AutoSave();
    file_output->Close();
}

int main(int argc, char *argv[]) {

    if (argc != 9){
        std::cout<<"Error: 2 arguments expected."<<std::endl;
        return 1; // return>0 indicates error code 
    }
    else {
        const std::string &filename_input = argv[1];
        const std::string &filename_output = argv[2];
        double nhitMinPrompt = atof(argv[3]);
        double nhitMinLate = atof(argv[4]);
        double deltaT = atof(argv[5]);
        double deltaP = atof(argv[6]);
        const bool use_mc = argv[7];
        const int round_selection = atoi(argv[8]);

        if (round_selection==1)
            round_1_cuts(filename_input, filename_output, nhitMinPrompt, nhitMinLate, deltaT, use_mc);
        if (round_selection==2)
            round_2_cuts(filename_input, filename_output, nhitMinPrompt, nhitMinLate, deltaT, deltaP, use_mc);
        if (round_selection==3){
            round_1_cuts(filename_input, filename_output, nhitMinPrompt, nhitMinLate, deltaT, use_mc);
            round_2_cuts(filename_input, filename_output, nhitMinPrompt, nhitMinLate, deltaT, deltaP, use_mc);
        }
        return 0; // completed successfully
    }
}
