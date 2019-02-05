#include <TFile.h>
#include <RAT/DB.hh>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <string>
#include <TNtuple.h>
#include <iostream>
#include <RAT/DU/Utility.hh>
#include <TObject.h>
#include <math.h>

void process_cuts(const std::string filename_input, const std::string filename_output, double energy_ep_min, double energy_ep_max, double energy_n_min, double energy_n_max, double deltaT = 1000000, double deltaP = 6000, double deltaD = 6000, const bool use_mc = false){

    TFile *file_input = new TFile(filename_input.c_str());
    TTree *tree_input = (TTree*)file_input->Get("tt");

    // load neutron pdf
    TFile *file_input_n = new TFile("/data/snoplus/lidgard/OXSX/flux100/myNeutron_flux100.root");
    TH2D *h2_after_cut_fit_neutron = (TH2D*)file_input_n->Get("h2_after_cut_efit_neutron");
    TFile *file_output = new TFile(filename_output.c_str(), "RECREATE");
    TTree *tree_output = new TTree("tt", "AntinuTree");

    const unsigned int n_slices = 100;
    char *name = new char[1000];
    TH1D *h_after_cut_fit_neutron[n_slices];
    //TF1 *f_after_cut_fit_neutron[n_slices];
    for (unsigned int i=0; i<n_slices; i++){
        sprintf(name, "h_after_cut_fit_neutron%d",i);
        h_after_cut_fit_neutron[i] = h2_after_cut_fit_neutron->ProjectionY(name,i*1-3,i*1+3);
        sprintf(name, "h_after_cut_fit_neutron%d ",i);
        h_after_cut_fit_neutron[i]->SetTitle(name);
        //if (h_after_cut_fit_neutron[i]->GetEntries()>0){
        //    sprintf(name, "f_after_cut_fit_neutron%d",i);
            //f_after_cut_fit_neutron[i] = new TF1(name,"gaus",0,3);
            //h_after_cut_fit_neutron[i]->Fit(f_after_cut_fit_neutron[i], "RQN");
        //}
    }

    size_t n_passed = 0;
    Int_t day_p1, day_p2;
    Int_t sec_p1, sec_p2;
    Int_t nsec_p1, nsec_p2;
    //Int_t nhit_p1, nhit_p2;
    bool fit_validity;
    double time_ns_diff;

    //float position_r, position_x, position_y, position_z;
    float position_ep_x;//, position_ep_y, position_ep_z, position_ep_r;
    float position_n_x;//, position_n_y, position_n_z, position_n_r;
    float position_p1_r, position_p1_x, position_p1_y, position_p1_z;
    float position_p2_r, position_p2_x, position_p2_y, position_p2_z;
    float energy_p1, energy_p1_corrected, energy_p2, energy_ep, energy_n, energy_nu, ev_particle_distance;
    bool all_pass, energy_pass, coincidence_pass, time_window_pass, position_r_pass, particle_distance_pass, any_ev_passed, ev_passed;

    //TH1D h_before_cut_nhit_p1("h_before_cut_nhit_p1", "h_before_cut_nhit_p1", 1000, 0, 5000);
    //TH1D h_before_cut_nhit_p2("h_before_cut_nhit_p2", "h_before_cut_nhit_p2", 1000, 0, 5000);
    //TH2D h2_before_cut_nhit("h2_before_cut_nhit", "h2_before_cut_nhit", 1000, 0, 5000, 1000, 0, 5000);
    TH2D h2_before_cut_emc("h2_before_cut_emc", "h2_before_cut_emc", 10000, 0, 10, 1000, 0, 1);
    TH1D h_before_cut_emc("h_before_cut_emc", "h_before_cut_emc", 1000, 0, 10);
    TH1D h_before_cut_emc_nu("h_before_cut_emc_nu", "h_before_cut_emc_nu", 1000, 0, 10);
    TH2D h2_before_cut_efit("h2_before_cut_efit", "h2_before_cut_efit", 1000, 0, 10, 1000, 0, 10);
    TH2D h2_before_cut_efit_corrected("h2_before_cut_efit_corrected", "h2_before_cut_efit_corrected", 1000, 0, 10, 1000, 0, 10);
    TH1D h_before_cut_efit_p1_corrected("h_before_cut_efit_p1_corrected", "h_before_cut_efit_p1_corrected", 1000, 0, 10);
    TH2D h2_before_cut_efit_neutron("h2_before_cut_efit_neutron", "h2_before_cut_efit_neutron", 100, 0, 10, 300, 0, 3);
    TH1D h_before_cut_efit_p2("h_before_cut_efit_p2", "h_before_cut_efit_p2", 1000, 0, 10);
    TH1D h_before_cut_time_diff("h_before_cut_time_diff", "h_before_cut_time_diff", 1000, 0, 5000000);
    TH2D h2_before_cut_energy_resolution("h2_before_cut_energy_resolution", "h2_before_cut_energy_resolution", 200, -2, 2, 201, -2, 2);
    TH2D h2_before_cut_position_resolution("h2_before_cut_position_resolution", "h2_before_cut_position_resolution", 200, -2, 2, 201, -2, 2);
    TH1D h_before_cut_position_displacement("h_before_cut_position_displacement", "h_before_cut_position_displacement", 1000, 0, 10000);
    TH2D h2_before_cut_time_diff_displacement("h2_before_cut_time_diff_displacement", "h2_before_cut_time_diff_displacement", 1000, 0, 5000000, 1000, 0, 10000);

    //TH1D h_after_cut_nhit_p1("h_after_cut_nhit_p1", "h_after_cut_nhit_p1", 1000, 0, 5000);
    //TH1D h_after_cut_nhit_p2("h_after_cut_nhit_p2", "h_after_cut_nhit_p2", 1000, 0, 5000);
    //TH2D h2_after_cut_nhit("h2_after_cut_nhit", "h2_after_cut_nhit", 1000, 0, 5000, 1000, 0, 5000);
    TH2D h2_after_cut_emc("h2_after_cut_emc", "h2_after_cut_emc", 10000, 0, 10, 1000, 0, 1);
    TH1D h_after_cut_emc("h_after_cut_emc", "h_after_cut_emc", 1000, 0, 10);
    TH1D h_after_cut_emc_nu("h_after_cut_emc_nu", "h_after_cut_emc_nu", 1000, 0, 10);
    TH2D h2_after_cut_efit("h2_after_cut_efit", "h2_after_cut_efit", 1000, 0, 10, 1000, 0, 10);
    TH2D h2_after_cut_efit_corrected("h2_after_cut_efit_corrected", "h2_after_cut_efit_corrected", 1000, 0, 10, 1000, 0, 10);
    TH1D h_after_cut_efit_p1_corrected("h_after_cut_efit_p1_corrected", "h_after_cut_efit_p1_corrected", 1000, 0, 10);
    TH2D h2_after_cut_efit_neutron("h2_after_cut_efit_neutron", "h2_after_cut_efit_neutron", 100, 0, 10, 300, 0, 3);
    TH1D h_after_cut_efit_p2("h_after_cut_efit_p2", "h_after_cut_efit_p2", 1000, 0, 10);
    TH1D h_after_cut_time_diff("h_after_cut_time_diff", "h_after_cut_time_diff", 1000, 0, 5000000);
    TH2D h2_after_cut_energy_resolution("h2_after_cut_energy_resolution", "h2_after_cut_energy_resolution", 200, -2, 2, 201, -2, 2);
    TH2D h2_after_cut_position_resolution("h2_after_cut_position_resolution", "h2_after_cut_position_resolution", 200, -2, 2, 201, -2, 2);
    TH1D h_after_cut_position_displacement("h_after_cut_position_displacement", "h_after_cut_position_displacement", 1000, 0, 10000);
    TH2D h2_after_cut_time_diff_displacement("h2_after_cut_time_diff_displacement", "h2_after_cut_time_diff_displacement", 1000, 0, 5000000, 1000, 0, 10000);

    double neutron_capture_energy = 1.87;
    double e_rem = 0.784; //1.374 +0.47

    std::vector<float> output_mc_quench;
    std::vector<float> output_mc_neutrino_energy;
    std::vector<float> output_mc_positron_energy;
    std::vector<float> output_mc_neutron_energy;
    std::vector<float> output_mc_neutrino_position_r;
    std::vector<float> output_mc_neutrino_position_x;
    std::vector<float> output_mc_neutrino_position_y;
    std::vector<float> output_mc_neutrino_position_z;
    std::vector<float> output_mc_positron_position_r;
    std::vector<float> output_mc_positron_position_x;
    std::vector<float> output_mc_positron_position_y;
    std::vector<float> output_mc_positron_position_z;
    std::vector<float> output_mc_neutron_position_r;
    std::vector<float> output_mc_neutron_position_x;
    std::vector<float> output_mc_neutron_position_y;
    std::vector<float> output_mc_neutron_position_z;
    std::vector<bool> output_ev_fit_validity;
    std::vector<float> output_ev_fit_energy;
    std::vector<float> output_ev_fit_position_r;
    std::vector<float> output_ev_fit_position_x;
    std::vector<float> output_ev_fit_position_y;
    std::vector<float> output_ev_fit_position_z;
    std::vector<int> output_ev_nhit;
    std::vector<unsigned int> output_ev_time_days;
    std::vector<unsigned int> output_ev_time_seconds;
    std::vector<float> output_ev_time_nanoseconds;
    std::vector<float> output_reactor_info_latitude;
    std::vector<float> output_reactor_info_longitude;
    std::vector<float> output_reactor_info_altitude;
    std::vector<float> output_reactor_info_distance;
    tree_output->Branch("mc_quench",&output_mc_quench);
    tree_output->Branch("mc_neutrino_energy",&output_mc_neutrino_energy);
    tree_output->Branch("mc_positron_energy",&output_mc_positron_energy);
    tree_output->Branch("mc_neutron_energy",&output_mc_neutron_energy);
    tree_output->Branch("mc_neutrino_position_r",&output_mc_neutrino_position_r);
    tree_output->Branch("mc_neutrino_position_x",&output_mc_neutrino_position_x);
    tree_output->Branch("mc_neutrino_position_y",&output_mc_neutrino_position_y);
    tree_output->Branch("mc_neutrino_position_z",&output_mc_neutrino_position_z);
    tree_output->Branch("mc_positron_position_r",&output_mc_positron_position_r);
    tree_output->Branch("mc_positron_position_x",&output_mc_positron_position_x);
    tree_output->Branch("mc_positron_position_y",&output_mc_positron_position_y);
    tree_output->Branch("mc_positron_position_z",&output_mc_positron_position_z);
    tree_output->Branch("mc_neutron_position_r",&output_mc_neutron_position_r);
    tree_output->Branch("mc_neutron_position_x",&output_mc_neutron_position_x);
    tree_output->Branch("mc_neutron_position_y",&output_mc_neutron_position_y);
    tree_output->Branch("mc_neutron_position_z",&output_mc_neutron_position_z);
    tree_output->Branch("ev_fit_validity",&output_ev_fit_validity);
    tree_output->Branch("ev_fit_energy",&output_ev_fit_energy);
    tree_output->Branch("ev_fit_position_r",&output_ev_fit_position_r);
    tree_output->Branch("ev_fit_position_x",&output_ev_fit_position_x);
    tree_output->Branch("ev_fit_position_y",&output_ev_fit_position_y);
    tree_output->Branch("ev_fit_position_z",&output_ev_fit_position_z);
    tree_output->Branch("ev_nhit",&output_ev_nhit);
    tree_output->Branch("ev_time_days",&output_ev_time_days);
    tree_output->Branch("ev_time_seconds",&output_ev_time_seconds);
    tree_output->Branch("ev_time_nanoseconds",&output_ev_time_nanoseconds);
    tree_output->Branch("reactor_info_latitude",&output_reactor_info_latitude);
    tree_output->Branch("reactor_info_longitude",&output_reactor_info_longitude);
    tree_output->Branch("reactor_info_altitude",&output_reactor_info_altitude);
    tree_output->Branch("reactor_info_distance",&output_reactor_info_distance);

    //size_t ientry=0;
    std::vector<float> *mc_quench = 0;
    std::vector<float> *mc_neutrino_energy = 0;
    std::vector<float> *mc_positron_energy = 0;
    std::vector<float> *mc_neutron_energy = 0;
    std::vector<float> *mc_neutrino_position_r = 0;
    std::vector<float> *mc_neutrino_position_x = 0;
    std::vector<float> *mc_neutrino_position_y = 0;
    std::vector<float> *mc_neutrino_position_z = 0;
    std::vector<float> *mc_positron_position_r = 0;
    std::vector<float> *mc_positron_position_x = 0;
    std::vector<float> *mc_positron_position_y = 0;
    std::vector<float> *mc_positron_position_z = 0;
    std::vector<float> *mc_neutron_position_r = 0;
    std::vector<float> *mc_neutron_position_x = 0;
    std::vector<float> *mc_neutron_position_y = 0;
    std::vector<float> *mc_neutron_position_z = 0;
    std::vector<bool> *ev_fit_validity = 0;
    std::vector<float> *ev_fit_energy = 0;
    std::vector<float> *ev_fit_position_r = 0;
    std::vector<float> *ev_fit_position_x = 0;
    std::vector<float> *ev_fit_position_y = 0;
    std::vector<float> *ev_fit_position_z = 0;
    std::vector<int> *ev_nhit = 0;
    std::vector<unsigned int> *ev_time_days = 0;
    std::vector<unsigned int> *ev_time_seconds = 0;
    std::vector<float> *ev_time_nanoseconds = 0;
    std::vector<float> *reactor_info_latitude = 0;
    std::vector<float> *reactor_info_longitude = 0;
    std::vector<float> *reactor_info_altitude = 0;
    std::vector<float> *reactor_info_distance = 0;

    tree_input->SetBranchAddress("mc_quench",&mc_quench);
    tree_input->SetBranchAddress("mc_neutrino_energy",&mc_neutrino_energy);
    tree_input->SetBranchAddress("mc_positron_energy",&mc_positron_energy);
    tree_input->SetBranchAddress("mc_neutron_energy",&mc_neutron_energy);
    tree_input->SetBranchAddress("mc_neutrino_position_r",&mc_neutrino_position_r);
    tree_input->SetBranchAddress("mc_neutrino_position_x",&mc_neutrino_position_x);
    tree_input->SetBranchAddress("mc_neutrino_position_y",&mc_neutrino_position_y);
    tree_input->SetBranchAddress("mc_neutrino_position_z",&mc_neutrino_position_z);
    tree_input->SetBranchAddress("mc_positron_position_r",&mc_positron_position_r);
    tree_input->SetBranchAddress("mc_positron_position_x",&mc_positron_position_x);
    tree_input->SetBranchAddress("mc_positron_position_y",&mc_positron_position_y);
    tree_input->SetBranchAddress("mc_positron_position_z",&mc_positron_position_z);
    tree_input->SetBranchAddress("mc_neutron_position_r",&mc_neutron_position_r);
    tree_input->SetBranchAddress("mc_neutron_position_x",&mc_neutron_position_x);
    tree_input->SetBranchAddress("mc_neutron_position_y",&mc_neutron_position_y);
    tree_input->SetBranchAddress("mc_neutron_position_z",&mc_neutron_position_z);
    tree_input->SetBranchAddress("ev_fit_validity",&ev_fit_validity);
    tree_input->SetBranchAddress("ev_fit_energy",&ev_fit_energy);
    tree_input->SetBranchAddress("ev_fit_position_r",&ev_fit_position_r);
    tree_input->SetBranchAddress("ev_fit_position_x",&ev_fit_position_x);
    tree_input->SetBranchAddress("ev_fit_position_y",&ev_fit_position_y);
    tree_input->SetBranchAddress("ev_fit_position_z",&ev_fit_position_z);
    tree_input->SetBranchAddress("ev_nhit",&ev_nhit);
    tree_input->SetBranchAddress("ev_time_days",&ev_time_days);
    tree_input->SetBranchAddress("ev_time_seconds",&ev_time_seconds);
    tree_input->SetBranchAddress("ev_time_nanoseconds",&ev_time_nanoseconds);
    tree_input->SetBranchAddress("reactor_info_latitude",&reactor_info_latitude);
    tree_input->SetBranchAddress("reactor_info_longitude",&reactor_info_longitude);
    tree_input->SetBranchAddress("reactor_info_altitude",&reactor_info_altitude);
    tree_input->SetBranchAddress("reactor_info_distance",&reactor_info_distance);

    size_t n_entries = tree_input->GetEntries();
    int i_percent = (n_entries/100)*10;
    std::cout<<"initial entries: "<<n_entries<<std::endl;
    for (size_t i=0; i < (n_entries-1); i++){

        // print progress
        if(i_percent && !(i%i_percent))
            std::cout << i/i_percent*10 << "% done " << std::endl;

        tree_input->GetEntry(i);
        any_ev_passed = false;

        // reset vectors to zero length for each event
        output_mc_quench.clear();
        output_mc_neutrino_energy.clear();
        output_mc_positron_energy.clear();
        output_mc_neutron_energy.clear();
        output_mc_neutrino_position_r.clear();
        output_mc_neutrino_position_x.clear();
        output_mc_neutrino_position_y.clear();
        output_mc_neutrino_position_z.clear();
        output_mc_positron_position_r.clear();
        output_mc_positron_position_x.clear();
        output_mc_positron_position_y.clear();
        output_mc_positron_position_z.clear();
        output_mc_neutron_position_r.clear();
        output_mc_neutron_position_x.clear();
        output_mc_neutron_position_y.clear();
        output_mc_neutron_position_z.clear();
        output_ev_fit_validity.clear();
        output_ev_fit_energy.clear();
        output_ev_fit_position_r.clear();
        output_ev_fit_position_x.clear();
        output_ev_fit_position_y.clear();
        output_ev_fit_position_z.clear();
        output_ev_nhit.clear();
        output_ev_time_days.clear();
        output_ev_time_seconds.clear();
        output_ev_time_nanoseconds.clear();
        output_reactor_info_latitude.clear();
        output_reactor_info_longitude.clear();
        output_reactor_info_altitude.clear();
        output_reactor_info_distance.clear();

        // get mc information
        //position_ep_r = mc_positron_position_r->at(0);
        position_ep_x = mc_positron_position_x->at(0);
        //position_ep_y = mc_positron_position_y->at(0);
        //position_ep_z = mc_positron_position_z->at(0);
        //position_n_r = mc_neutron_position_r->at(0);
        position_n_x = mc_neutron_position_x->at(0);
        //position_n_y = mc_neutron_position_y->at(0);
        //position_n_z = mc_neutron_position_z->at(0);
        energy_ep = mc_positron_energy->at(0);
        energy_n = mc_neutron_energy->at(0);
        energy_nu = mc_neutrino_energy->at(0);

        // this should never print. i.e. in MC data there is only 1 parent.
        if (mc_positron_position_r->size() > 1)
            std::cout << i << "more than one MC parent: " << mc_positron_position_r->size() << std::endl;

        // first test all ev entries, then test between entries
        size_t n_ev = ev_fit_validity->size();
        for (size_t j=0; j < n_ev; j++){ // look through all triggers

            ev_passed = false; // if this trigger (j) is in coincidence with any other trigger (and passes cuts, then we will mark this as true)

            // get data from vectors for this trigger
            fit_validity = ev_fit_validity->at(j);
            day_p1 = ev_time_days->at(j);
            sec_p1 = ev_time_seconds->at(j);
            nsec_p1 = ev_time_nanoseconds->at(j);
            position_p1_r = ev_fit_position_r->at(j);
            position_p1_x = ev_fit_position_x->at(j);
            position_p1_y = ev_fit_position_y->at(j);
            position_p1_z = ev_fit_position_z->at(j);
            energy_p1 = ev_fit_energy->at(j);
            //nhit_p1 = ev_nhit->at(j);

            for (size_t k=0; k < n_ev; k++){ // compare with each other trigger, starting with the next trigger up to the final trigger

                if (k==j) continue; // test trigger for coincidence against every other trigger except itself

                // reset pass flags
                all_pass = false;
                coincidence_pass = false;
                time_window_pass = false;
                position_r_pass = false;
                energy_pass = false;
                particle_distance_pass = false;

                // get data from vectors for this trigger
                day_p2 = ev_time_days->at(k);
                sec_p2 = ev_time_seconds->at(k);
                nsec_p2 = ev_time_nanoseconds->at(k);
                position_p2_r = ev_fit_position_r->at(k);
                position_p2_x = ev_fit_position_x->at(k);
                position_p2_y = ev_fit_position_y->at(k);
                position_p2_z = ev_fit_position_z->at(k);
                energy_p2 = ev_fit_energy->at(k);
                //nhit_p2 = ev_nhit->at(k);

                // inter particle distance
                TVector3 position_p1_pos(position_p1_x, position_p1_y, position_p1_z);
                TVector3 position_p2_pos(position_p2_x, position_p2_y, position_p2_z);
                ev_particle_distance = (position_p2_pos - position_p1_pos).Mag();

                // calculate time difference
                if(std::abs(day_p2 - day_p1) > 0){ //begin time tests
                    coincidence_pass = false;
                    time_ns_diff = std::abs(day_p2 - day_p1)*24*60*60*1e9 + std::abs(sec_p2 - sec_p1)*1e9 + std::abs(nsec_p2 - nsec_p1);
                }
                else { // days are equal
                    if(std::abs(sec_p2 - sec_p1) > 0){
                        coincidence_pass = false;
                        time_ns_diff = std::abs(sec_p2 - sec_p1)*1e9 + std::abs(nsec_p2 - nsec_p1);
                    }
                    else { // seconds are equal
                        time_ns_diff = std::abs(nsec_p2 - nsec_p1);
                        if(std::abs(nsec_p2 - nsec_p1) < deltaT) // nanosec are to within deltaT
                            coincidence_pass = true;
                        else
                            coincidence_pass = false; // so close..
                    }
                }

                // conditions for event to pass cuts:
                // both particles to be fitted within radius
                if ((position_p1_r < deltaP)&&(position_p2_r < deltaP))
                    position_r_pass = true;

                // inter-particle distance cut
                if (ev_particle_distance < deltaD)
                    particle_distance_pass = true;

                // time window cut
                time_window_pass = coincidence_pass; //(coincidence_pass_previous || coincidence_pass);

                // energy cut
                if (nsec_p1 <= nsec_p2){ // p1 is first (prompt particle i.e. assumed positron)
                    if (((energy_p1+e_rem) > energy_ep_min) && ((energy_p1+e_rem) < energy_ep_max) && (energy_p2 > energy_n_min) && (energy_p2 < energy_n_max))
                        energy_pass = true;
                }
                else{
                    if (((energy_p2+e_rem) > energy_ep_min) && ((energy_p2+e_rem) < energy_ep_max) && (energy_p1 > energy_n_min) && (energy_p1 < energy_n_max))
                        energy_pass = true;
                }

                all_pass = fit_validity && time_window_pass && position_r_pass && particle_distance_pass && energy_pass;

                if (all_pass){
                    ev_passed = true; // if any trigger passes, throw a flag (flag resets once per 'j')
                    any_ev_passed = true; // if any ev passes, throw a flags (flag reset once per 'i')

                    // fill histograms with all triggers which pass cuts
                    if (nsec_p1 <= nsec_p2){ // take care of 'reversed' nanosecond times
                        // if there is nothing in the projection, then return zero for e_n
                        double e_n = 0;
                        unsigned int my_n = round(energy_p1*10);
                        if (my_n>=n_slices) {
                            //std::cout<< "warning: slice too high (nsec_p1 <= nsec_p2): " << my_n << " e_p1: " << energy_p1 << std::endl;
                            my_n = n_slices-1;
                        }
                        if (my_n<0) my_n = 0;
                        if (h_after_cut_fit_neutron[my_n]->GetEntries()>0)
                            e_n = h_after_cut_fit_neutron[my_n]->GetRandom();
                        //h_after_cut_nhit_p1.Fill(nhit_p1);
                        //h_after_cut_nhit_p2.Fill(nhit_p2);
                        //h2_after_cut_nhit.Fill(nhit_p1, nhit_p2);
                        h2_after_cut_efit.Fill(energy_p1, energy_p2);
                        energy_p1_corrected = energy_p1+e_rem+e_n;
                        h2_after_cut_efit_corrected.Fill(energy_p1_corrected, energy_p2);
                        h_after_cut_efit_p1_corrected.Fill(energy_p1_corrected);
                        h2_after_cut_efit_neutron.Fill(energy_p1, energy_nu - (energy_p1+e_rem) );
                        h_after_cut_efit_p2.Fill(energy_p2);
                        h_after_cut_time_diff.Fill(time_ns_diff);
                        h_after_cut_position_displacement.Fill(ev_particle_distance);
                        h2_after_cut_time_diff_displacement.Fill(time_ns_diff, ev_particle_distance);
                        h2_after_cut_energy_resolution.Fill((energy_p1_corrected-energy_nu)/energy_nu, (energy_p2-neutron_capture_energy)/neutron_capture_energy);
                        h2_after_cut_position_resolution.Fill((position_p1_x-position_ep_x)/position_ep_x, (position_p2_x-position_n_x)/position_n_x);
                    }
                    else{
                        double e_n = 0;
                        unsigned int my_n = round(energy_p2*10);
                        if (my_n>=n_slices) {
                            //std::cout<< "warning: slice too high (nsec_p1 > nsec_p2): " << my_n << " e_p1: " << energy_p2 << std::endl;
                            my_n = n_slices-1;
                        }
                        if (my_n<0) my_n = 0;
                        if (h_after_cut_fit_neutron[my_n]->GetEntries()>0)
                            e_n = h_after_cut_fit_neutron[my_n]->GetRandom();
                        //h_after_cut_nhit_p1.Fill(nhit_p2);
                        //h_after_cut_nhit_p2.Fill(nhit_p1);
                        //h2_after_cut_nhit.Fill(nhit_p2, nhit_p1);
                        h2_after_cut_efit.Fill(energy_p2, energy_p1);
                        energy_p1_corrected = energy_p2+e_rem+e_n;
                        h2_after_cut_efit_corrected.Fill(energy_p1_corrected, energy_p1);
                        h_after_cut_efit_p1_corrected.Fill(energy_p1_corrected);
                        h2_after_cut_efit_neutron.Fill(energy_p2, energy_nu - (energy_p2+e_rem) );
                        h_after_cut_efit_p2.Fill(energy_p1);
                        h_after_cut_time_diff.Fill(time_ns_diff);
                        h_after_cut_position_displacement.Fill(ev_particle_distance);
                        h2_after_cut_time_diff_displacement.Fill(time_ns_diff, ev_particle_distance);
                        h2_after_cut_energy_resolution.Fill((energy_p1_corrected-energy_nu)/energy_nu, (energy_p1-neutron_capture_energy)/neutron_capture_energy);
                        h2_after_cut_position_resolution.Fill((position_p2_x-position_ep_x)/position_ep_x, (position_p1_x-position_n_x)/position_n_x);
                    }
                }

                // fill histograms with every trigger
                if (nsec_p1 <= nsec_p2){ // take care of 'reversed' nanosecond times
                    double e_n = 0;
                    unsigned int my_n = round(energy_p1*10);
                    if (my_n>=n_slices) {
                        //std::cout<< "warning: slice too high (nsec_p1 <= nsec_p2): " << my_n << " e_p1: " << energy_p1 << std::endl;
                        my_n = n_slices-1;
                    }
                    if (my_n<0) my_n = 0;
                    if (h_after_cut_fit_neutron[my_n]->GetEntries()>0)
                        e_n = h_after_cut_fit_neutron[my_n]->GetRandom();
                    //h_before_cut_nhit_p1.Fill(nhit_p1);
                    //h_before_cut_nhit_p2.Fill(nhit_p2);
                    //h2_before_cut_nhit.Fill(nhit_p1, nhit_p2);
                    h2_before_cut_efit.Fill(energy_p1, energy_p2);
                    energy_p1_corrected = energy_p1+e_rem+e_n;
                    h2_before_cut_efit_corrected.Fill(energy_p1_corrected, energy_p2);
                    h_before_cut_efit_p1_corrected.Fill(energy_p1_corrected);
                    h2_before_cut_efit_neutron.Fill(energy_p1, energy_nu - (energy_p1+e_rem) );
                    h_before_cut_efit_p2.Fill(energy_p2);
                    h_before_cut_time_diff.Fill(time_ns_diff);
                    h_before_cut_position_displacement.Fill(ev_particle_distance);
                    h2_before_cut_time_diff_displacement.Fill(time_ns_diff, ev_particle_distance);
                    h2_before_cut_energy_resolution.Fill((energy_p1_corrected-energy_nu)/energy_nu, (energy_p2-neutron_capture_energy)/neutron_capture_energy);
                    h2_before_cut_position_resolution.Fill((position_p1_x-position_ep_x)/position_ep_x, (position_p2_x-position_n_x)/position_n_x);
                }
                else{
                    double e_n = 0;
                    unsigned int my_n = round(energy_p2*10);
                    if (my_n>=n_slices) {
                        //std::cout<< "warning: slice too high (nsec_p1 <= nsec_p2): " << my_n << " e_p1: " << energy_p2 << std::endl;
                        my_n = n_slices-1;
                    }
                    if (my_n<0) my_n = 0;
                    if (h_after_cut_fit_neutron[my_n]->GetEntries()>0)
                       e_n = h_after_cut_fit_neutron[my_n]->GetRandom();
                    //h_before_cut_nhit_p1.Fill(nhit_p2);
                    //h_before_cut_nhit_p2.Fill(nhit_p1);
                    //h2_before_cut_nhit.Fill(nhit_p2, nhit_p1);
                    h2_before_cut_efit.Fill(energy_p2, energy_p1);
                    energy_p1_corrected = energy_p2+e_rem+e_n;
                    h2_before_cut_efit_corrected.Fill(energy_p1_corrected, energy_p1);
                    h_before_cut_efit_p1_corrected.Fill(energy_p1_corrected);
                    h2_before_cut_efit_neutron.Fill(energy_p2, energy_nu - (energy_p2+e_rem) );
                    h_before_cut_efit_p2.Fill(energy_p1);
                    h_before_cut_time_diff.Fill(time_ns_diff);
                    h_before_cut_position_displacement.Fill(ev_particle_distance);
                    h2_before_cut_time_diff_displacement.Fill(time_ns_diff, ev_particle_distance);
                    h2_before_cut_energy_resolution.Fill((energy_p1_corrected-energy_nu)/energy_nu, (energy_p1-neutron_capture_energy)/neutron_capture_energy);
                    h2_before_cut_position_resolution.Fill((position_p2_x-position_ep_x)/position_ep_x, (position_p1_x-position_n_x)/position_n_x);
                }
            }

            //save data... ev's go into vectors
            //if passed then push events into vectors and mark event to be added to tree
            if (ev_passed){
                output_ev_fit_validity.push_back(ev_fit_validity->at(j));
                output_ev_fit_energy.push_back(ev_fit_energy->at(j));
                output_ev_fit_position_r.push_back(ev_fit_position_r->at(j));
                output_ev_fit_position_x.push_back(ev_fit_position_x->at(j));
                output_ev_fit_position_y.push_back(ev_fit_position_y->at(j));
                output_ev_fit_position_z.push_back(ev_fit_position_z->at(j));
                output_ev_nhit.push_back(ev_nhit->at(j));
                output_ev_time_days.push_back(ev_time_days->at(j));
                output_ev_time_seconds.push_back(ev_time_seconds->at(j));
                output_ev_time_nanoseconds.push_back(ev_time_nanoseconds->at(j));
            }
        }

        // if at least one ev passed, then fill the output tree
        if (any_ev_passed){
            // only push back mc information once (only ever 1 parent particle in mc)
            output_mc_quench.push_back(mc_quench->at(0));
            output_mc_neutrino_energy.push_back(mc_neutrino_energy->at(0));
            output_mc_positron_energy.push_back(mc_positron_energy->at(0));
            output_mc_neutron_energy.push_back(mc_neutron_energy->at(0));
            output_mc_neutrino_position_r.push_back(mc_neutrino_position_r->at(0));
            output_mc_neutrino_position_x.push_back(mc_neutrino_position_x->at(0));
            output_mc_neutrino_position_y.push_back(mc_neutrino_position_y->at(0));
            output_mc_neutrino_position_z.push_back(mc_neutrino_position_z->at(0));
            output_mc_positron_position_r.push_back(mc_positron_position_r->at(0));
            output_mc_positron_position_x.push_back(mc_positron_position_x->at(0));
            output_mc_positron_position_y.push_back(mc_positron_position_y->at(0));
            output_mc_positron_position_z.push_back(mc_positron_position_z->at(0));
            output_mc_neutron_position_r.push_back(mc_neutron_position_r->at(0));
            output_mc_neutron_position_x.push_back(mc_neutron_position_x->at(0));
            output_mc_neutron_position_y.push_back(mc_neutron_position_y->at(0));
            output_mc_neutron_position_z.push_back(mc_neutron_position_z->at(0));
            output_reactor_info_latitude.push_back(reactor_info_latitude->at(0));
            output_reactor_info_longitude.push_back(reactor_info_longitude->at(0));
            output_reactor_info_altitude.push_back(reactor_info_altitude->at(0));
            output_reactor_info_distance.push_back(reactor_info_distance->at(0));

            // fill tree
            tree_output->Fill();
            n_passed++;
        }

        // fill mc event histograms (only one parent per event)
        h_before_cut_emc_nu.Fill(energy_nu);
        h_before_cut_emc.Fill(energy_ep);
        h2_before_cut_emc.Fill(energy_ep, energy_n);
        if (any_ev_passed){ // this event had a trigger which passed
            h_after_cut_emc_nu.Fill(energy_nu);
            h_after_cut_emc.Fill(energy_ep);
            h2_after_cut_emc.Fill(energy_ep, energy_n);
        }
    }

    // cut purity plots (ratio plot)
    TH1D *h_after_cut_emc_nu_ratio = (TH1D*)h_after_cut_emc_nu.Clone("h_after_cut_emc_nu_ratio");
    h_after_cut_emc_nu_ratio->Divide(&h_before_cut_emc_nu);

    TH1D *h_before_cut_efit_ratio = (TH1D*)h_before_cut_efit_p1_corrected.Clone("h_before_cut_efit_ratio");
    h_before_cut_efit_ratio->Divide(&h_before_cut_emc_nu);

    TH1D *h_after_cut_efit_ratio = (TH1D*)h_after_cut_efit_p1_corrected.Clone("h_after_cut_efit_ratio");
    h_after_cut_efit_ratio->Divide(&h_after_cut_emc_nu);

    TH1D *h_after_cut_emc_nu_n = (TH1D*)h_after_cut_emc_nu.Clone("h_after_cut_emc_nu_n");
    TH1D *h_after_cut_efit_p1_corrected_n = (TH1D*)h_after_cut_efit_p1_corrected.Clone("h_after_cut_efit_p1_corrected_n");
    //h_after_cut_emc_nu_n->Scale(1./h_after_cut_emc_nu_n->Integral());
    //h_after_cut_efit_p1_corrected_n->Scale(1./h_after_cut_efit_p1_corrected_n->Integral());
    h_after_cut_emc_nu_n->Scale(1./h_after_cut_emc_nu_n->GetMaximum());
    h_after_cut_efit_p1_corrected_n->Scale(1./h_after_cut_efit_p1_corrected_n->GetMaximum());
    h_after_cut_efit_p1_corrected_n->SetLineColor(kRed);

    // save file
    file_input->Close();
    file_output->cd();

    //h_before_cut_nhit_p1.GetXaxis()->SetTitle("nhit_p1");
    //h_before_cut_nhit_p2.GetXaxis()->SetTitle("nhit_p2");
    //h2_before_cut_nhit.GetXaxis()->SetTitle("nhit_p1");
    //h2_before_cut_nhit.GetYaxis()->SetTitle("nhit_p2");
    h2_before_cut_emc.GetXaxis()->SetTitle("energy_ep");
    h2_before_cut_emc.GetYaxis()->SetTitle("energy_n");
    h2_before_cut_efit.GetXaxis()->SetTitle("energy_p1");
    h2_before_cut_efit.GetYaxis()->SetTitle("energy_p2");
    h2_before_cut_efit_corrected.GetXaxis()->SetTitle("energy_p1+0.784MeV+E_n");
    h2_before_cut_efit_corrected.GetYaxis()->SetTitle("energy_p2");
    h2_before_cut_time_diff_displacement.GetXaxis()->SetTitle("time_ns_diff");
    h2_before_cut_time_diff_displacement.GetYaxis()->SetTitle("ev_particle_distance");
    h2_before_cut_energy_resolution.GetXaxis()->SetTitle("(energy_p1+0.784MeV+E_n-energy_nu)/energy_nu");
    h2_before_cut_energy_resolution.GetYaxis()->SetTitle("(energy_p2-n_capE)/n_capE");
    h2_before_cut_position_resolution.GetXaxis()->SetTitle("(position_p1_x-position_ep_x)/position_ep_x");
    h2_before_cut_position_resolution.GetYaxis()->SetTitle("(position_p2_x-position_n_x)/position_n_x");
    h2_before_cut_efit_neutron.GetXaxis()->SetTitle("energy_p1");
    h2_before_cut_efit_neutron.GetYaxis()->SetTitle("energy_nu - (energy_p1+e_rem)");
    h_before_cut_emc.SetTitle("KE_ep");
    h_before_cut_emc_nu.SetTitle("KE_nu");
    h_before_cut_efit_p1_corrected.SetTitle("energy_p1+0.784MeV+E_n");
    h_before_cut_efit_p2.SetTitle("energy_p2");
    h_before_cut_time_diff.SetTitle("time_ns_diff");
    h_before_cut_position_displacement.SetTitle("ev_particle_distance");
    h_before_cut_emc.GetXaxis()->SetTitle("Energy (MeV)");
    h_before_cut_emc_nu.GetXaxis()->SetTitle("Energy (MeV)");
    h_before_cut_efit_p1_corrected.GetXaxis()->SetTitle("Energy (MeV)");
    h_before_cut_efit_p2.GetXaxis()->SetTitle("Energy (MeV)");
    h_before_cut_time_diff.GetXaxis()->SetTitle("Time between abs(ev_fit_p2 - ev_fit_p1) (ns)");
    h_before_cut_position_displacement.GetXaxis()->SetTitle("Position (ev_fit_p1 - ev_fit_p2) (mm)");

    //h_after_cut_nhit_p1.GetXaxis()->SetTitle("nhit_p1");
    //h_after_cut_nhit_p2.GetXaxis()->SetTitle("nhit_p2");
    //h2_after_cut_nhit.GetXaxis()->SetTitle("nhit_p1");
    //h2_after_cut_nhit.GetYaxis()->SetTitle("nhit_p2");
    h2_after_cut_emc.GetXaxis()->SetTitle("energy_ep");
    h2_after_cut_emc.GetYaxis()->SetTitle("energy_n");
    h2_after_cut_efit.GetXaxis()->SetTitle("energy_p1");
    h2_after_cut_efit.GetYaxis()->SetTitle("energy_p2");
    h2_after_cut_efit_corrected.GetXaxis()->SetTitle("energy_p1+0.784MeV+E_n");
    h2_after_cut_efit_corrected.GetYaxis()->SetTitle("energy_p2");
    h2_after_cut_efit_neutron.GetXaxis()->SetTitle("energy_p1");
    h2_after_cut_efit_neutron.GetYaxis()->SetTitle("energy_nu - (energy_p1+e_rem)");
    h2_after_cut_time_diff_displacement.GetXaxis()->SetTitle("time_ns_diff");
    h2_after_cut_time_diff_displacement.GetYaxis()->SetTitle("ev_particle_distance");
    h2_after_cut_energy_resolution.GetXaxis()->SetTitle("(energy_p1+0.784MeV+E_n-energy_nu)/energy_nu");
    h2_after_cut_energy_resolution.GetYaxis()->SetTitle("(energy_p2+0.784MeV-n_capE)/n_capE");
    h2_after_cut_position_resolution.GetXaxis()->SetTitle("(position_p1_x-position_ep_x)/position_ep_x");
    h2_after_cut_position_resolution.GetYaxis()->SetTitle("(position_p2_x-position_n_x)/position_n_x");
    h_after_cut_emc.SetTitle("KE_ep");
    h_after_cut_emc_nu.SetTitle("KE_nu");
    h_after_cut_emc_nu_n->SetTitle("KE_nu (normalised (maximum = 1)");
    h_after_cut_efit_p1_corrected.SetTitle("energy_p1+0.784MeV+E_n");
    h_after_cut_efit_p2.SetTitle("energy_p2");
    h_after_cut_efit_p1_corrected_n->SetTitle("energy_p1+0.784MeV+E_n");
    h_after_cut_time_diff.SetTitle("time_ns_diff");
    h_after_cut_position_displacement.SetTitle("ev_particle_distance");
    h_after_cut_emc.GetXaxis()->SetTitle("Energy (MeV)");
    h_after_cut_emc_nu.GetXaxis()->SetTitle("Energy (MeV)");
    h_after_cut_emc_nu_n->GetXaxis()->SetTitle("Energy (MeV)");
    h_after_cut_efit_p1_corrected.GetXaxis()->SetTitle("Energy (MeV)");
    h_after_cut_efit_p2.GetXaxis()->SetTitle("Energy (MeV)");
    h_after_cut_efit_p1_corrected_n->GetXaxis()->SetTitle("Energy (MeV)");
    h_after_cut_time_diff.GetXaxis()->SetTitle("Time between abs(ev_fit_p2 - ev_fit_p1) (ns)");
    h_after_cut_position_displacement.GetXaxis()->SetTitle("Position (ev_fit_p1 - ev_fit_p2) (mm)");
    
    h_after_cut_emc_nu_ratio->GetYaxis()->SetTitle("KE_nu after cut / KE_nu before cut");
    h_before_cut_efit_ratio->GetYaxis()->SetTitle("energy_p1+0.784MeV+E_n / KE_nu");
    h_after_cut_efit_ratio->GetYaxis()->SetTitle("energy_p1+0.784MeV+E_n / KE_nu");
    h_after_cut_emc_nu_ratio->GetXaxis()->SetTitle("Energy (MeV)");
    h_before_cut_efit_ratio->GetXaxis()->SetTitle("Energy (MeV)");
    h_after_cut_efit_ratio->GetXaxis()->SetTitle("Energy (MeV)");

    //h_before_cut_nhit_p1.Write();
    //h_before_cut_nhit_p2.Write();
    //h2_before_cut_nhit.Write();
    h2_before_cut_emc.Write();
    h_before_cut_emc_nu.Write();
    h_before_cut_emc.Write();
    h2_before_cut_efit.Write();
    h2_before_cut_efit_corrected.Write();
    h_before_cut_efit_p1_corrected.Write();
    h2_before_cut_efit_neutron.Write();
    h_before_cut_efit_p2.Write();
    h_before_cut_time_diff.Write();
    h_before_cut_position_displacement.Write();
    h2_before_cut_time_diff_displacement.Write();
    h2_before_cut_energy_resolution.Write();
    h2_before_cut_position_resolution.Write();

    //h_after_cut_nhit_p1.Write();
    //h_after_cut_nhit_p1.Write();
    //h2_after_cut_nhit.Write();
    h2_after_cut_emc.Write();
    h_after_cut_emc_nu.Write();
    h_after_cut_emc_nu_n->Write();
    h_after_cut_emc.Write();
    h2_after_cut_efit.Write();
    h2_after_cut_efit_corrected.Write();
    h_after_cut_efit_p1_corrected.Write();
    h_after_cut_efit_p1_corrected_n->Write();
    h2_after_cut_efit_neutron.Write();
    h_after_cut_efit_p2.Write();
    h_after_cut_time_diff.Write();
    h_after_cut_position_displacement.Write();
    h2_after_cut_time_diff_displacement.Write();
    h2_after_cut_energy_resolution.Write();
    h2_after_cut_position_resolution.Write();

    h_after_cut_emc_nu_ratio->Write();
    h_before_cut_efit_ratio->Write();
    h_after_cut_efit_ratio->Write();

    //for (unsigned int ii=0; ii<n_slices; ii++)
    //    h_after_cut_fit_neutron[ii]->Write();

    tree_output->AutoSave();
    file_output->Close();
    std::cout<<"final entries: "<<n_passed<<std::endl;
}

int main(int argc, char *argv[]) {

    if (argc != 11) {
        std::cout<<"Error: 10 arguments expected. Got: "<<argc-1<<std::endl;
        return 1; // return>0 indicates error code
    }
    else {
        TH1::AddDirectory(kFALSE);
        const std::string &filename_input = argv[1];
        const std::string &filename_output = argv[2];
        double energy_ep_min = atof(argv[3]);
        double energy_ep_max = atof(argv[4]);
        double energy_n_min = atof(argv[5]);
        double energy_n_max = atof(argv[6]);
        double deltaT = atof(argv[7]);
        double deltaP = atof(argv[8]);
        double deltaD = atof(argv[9]);
        const bool use_mc = atoi(argv[10]);

        process_cuts(filename_input, filename_output, energy_ep_min, energy_ep_max, energy_n_min, energy_n_max, deltaT, deltaP, deltaD, use_mc);

        return 0; // completed successfully
    }
}
