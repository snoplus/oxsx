#include <TFile.h>
#include <RAT/DB.hh>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <string>
#include <TNtuple.h>
#include <iostream>
#include <RAT/DU/Utility.hh>
#include <TObject.h>

void process_cuts(const std::string filename_input, const std::string filename_output, double energy_ep_min, double energy_ep_max, double energy_n_min, double energy_n_max, double deltaT = 1000000, double deltaP = 6000, double deltaD = 6000, const bool use_mc = false){
    // round 1 checks for event pairs with time difference < 1000micros and both events within 6000mm
    // round 2 checks for event pairs with time difference < 275micros, both events within 5000mm and position difference between events (in pair) < 1250mm
    TFile *file_input = new TFile(filename_input.c_str());
    TTree *tree_input = (TTree*)file_input->Get("tt");
    TFile *file_output = new TFile(filename_output.c_str(), "RECREATE");
    TTree *tree_output = new TTree("tt", "AntinuTree");

    size_t n_entries = tree_input->GetEntries();
    size_t n_passed = 0;
    Int_t days, nextDays;
    Int_t secs, nextSec;
    Int_t nsecs, nextNSec;
    Int_t nhit_p1, nhit_p2;
    bool fit_validity;
    double time_ns_diff;
    //float position_r, position_x, position_y, position_z;
    float position_ep_x;//, position_ep_y, position_ep_z, position_ep_r;
    float position_n_x;//, position_n_y, position_n_z, position_n_r;
    float position_p1_r, position_p1_x, position_p1_y, position_p1_z;
    float position_p2_r, position_p2_x, position_p2_y, position_p2_z;
    float energy_p1, energy_p2, energy_ep, energy_n, ev_particle_distance;
    bool coincidence_pass_previous = false; // assume that the first event is invalidated by an event we didn't catch
    bool all_pass, energy_pass, coincidence_pass, time_window_pass, position_r_pass, particle_distance_pass, ev_passed;
    bool coincidencePartner=false;

    TH1D h_before_cut_nhit("h_before_cut_nhit", "h_before_cut_nhit", 3000, 0, 3000);
    TH1D h_before_cut_efit("h_before_cut_efit", "h_before_cut_efit", 1000, 0, 10);
    TH2D h2_before_cut_nhit("h2_before_cut_nhit", "h2_before_cut_nhit", 3000, 0, 3000, 3000, 0, 3000);
    TH2D h2_before_cut_emc("h2_before_cut_emc", "h2_before_cut_emc", 1000, 0, 1, 10000, 0, 10);
    TH2D h2_before_cut_efit("h2_before_cut_efit", "h2_before_cut_efit", 1000, 0, 10, 1000, 0, 10);
    TH1D h_before_cut_time_diff("h_before_cut_time_diff", "h_before_cut_time_diff", 1000, 0, 5000000);
    TH2D h2_before_cut_energy_resolution("h2_before_cut_energy_resolution", "h2_before_cut_energy_resolution", 200, -2, 2, 201, -2, 2);
    TH2D h2_before_cut_position_resolution("h2_before_cut_position_resolution", "h2_before_cut_position_resolution", 200, -2, 2, 201, -2, 2);
    TH1D h_before_cut_position_displacement("h_before_cut_position_displacement", "h_before_cut_position_displacement", 1000, 0, 10000);
    TH2D h2_before_cut_time_diff_displacement("h2_before_cut_time_diff_displacement", "h2_before_cut_time_diff_displacement", 1000, 0, 5000000, 1000, 0, 10000);

    TH1D h_after_cut_nhit("h_after_cut_nhit", "h_after_cut_nhit", 3000, 0, 3000);
    TH1D h_after_cut_efit("h_after_cut_efit", "h_after_cut_efit", 1000, 0, 10);
    TH2D h2_after_cut_nhit("h2_after_cut_nhit", "h2_after_cut_nhit", 3000, 0, 3000, 3000, 0, 3000);
    TH2D h2_after_cut_emc("h2_after_cut_emc", "h2_after_cut_emc", 1000, 0, 1, 10000, 0, 10);
    TH2D h2_after_cut_efit("h2_after_cut_efit", "h2_after_cut_efit", 1000, 0, 10, 1000, 0, 10);
    TH1D h_after_cut_time_diff("h_after_cut_time_diff", "h_after_cut_time_diff", 1000, 0, 5000000);
    TH2D h2_after_cut_energy_resolution("h2_after_cut_energy_resolution", "h2_after_cut_energy_resolution", 200, -2, 2, 201, -2, 2);
    TH2D h2_after_cut_position_resolution("h2_after_cut_position_resolution", "h2_after_cut_position_resolution", 200, -2, 2, 201, -2, 2);
    TH1D h_after_cut_position_displacement("h_after_cut_position_displacement", "h_after_cut_position_displacement", 1000, 0, 10000);
    TH2D h2_after_cut_time_diff_displacement("h2_after_cut_time_diff_displacement", "h2_after_cut_time_diff_displacement", 1000, 0, 5000000, 1000, 0, 10000);

    //size_t output_ientry=0;
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
    //tree_output->Branch("entry",output_ientry);
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
    //tree_input->SetBranchAddress("entry",&ientry);
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

    std::cout<<"initial entries: "<<n_entries<<std::endl;
    for (size_t i=0; i < (n_entries-1); i++){
        tree_input->GetEntry(i);
        ev_passed = false;

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

        size_t n_ev = ev_fit_validity->size();
        if (n_ev>0) n_ev--; // only minus 1 if there are >0 events. Need to -1 since we are testing pairs together (if i and i+1)
        
        // first test all ev entries, then test between entries
        for (size_t j=0; j < n_ev; j++){
            
            //std::cout<<"i "<<i<<" j "<<j<<" n_ev "<<n_ev<<std::endl;
            
            all_pass = false;
            coincidence_pass = false;
            time_window_pass = false;
            position_r_pass = false;
            energy_pass = false;
            particle_distance_pass = false;

            // this should never print. i.e. in MC data there is only 1 parent.
            if (mc_positron_position_r->size()>1)
                std::cout << i<<"more than one MC parent: "<<mc_positron_position_r->size()<< std::endl;
            
            // get data from vectors
            fit_validity = ev_fit_validity->at(j);
            days = ev_time_days->at(j);
            secs = ev_time_seconds->at(j);
            nsecs = ev_time_nanoseconds->at(j);

            //if (j < (ev_nhit->size()-1)){
            nextDays = ev_time_days->at(j+1);
            nextSec = ev_time_seconds->at(j+1);
            nextNSec = ev_time_nanoseconds->at(j+1);

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

            // take care of 'reversed' nanosecond times
            if (nextNSec>=nsecs){
                position_p1_r = ev_fit_position_r->at(j);
                position_p1_x = ev_fit_position_x->at(j);
                position_p1_y = ev_fit_position_y->at(j);
                position_p1_z = ev_fit_position_z->at(j);
                position_p2_r = ev_fit_position_r->at(j+1);
                position_p2_x = ev_fit_position_x->at(j+1);
                position_p2_y = ev_fit_position_y->at(j+1);
                position_p2_z = ev_fit_position_z->at(j+1);
                energy_p1 = ev_fit_energy->at(j);
                energy_p2 = ev_fit_energy->at(j+1);
                nhit_p1 = ev_nhit->at(j);
                nhit_p2 = ev_nhit->at(j+1);
            }
            else{
                position_p1_r = ev_fit_position_r->at(j+1);
                position_p1_x = ev_fit_position_x->at(j+1);
                position_p1_y = ev_fit_position_y->at(j+1);
                position_p1_z = ev_fit_position_z->at(j+1);
                position_p2_r = ev_fit_position_r->at(j);
                position_p2_x = ev_fit_position_x->at(j);
                position_p2_y = ev_fit_position_y->at(j);
                position_p2_z = ev_fit_position_z->at(j);
                energy_p1 = ev_fit_energy->at(j+1);
                energy_p2 = ev_fit_energy->at(j);
                nhit_p1 = ev_nhit->at(j+1);
                nhit_p2 = ev_nhit->at(j);
            }
            
            // calculate particle pair distance
            //TVector3 position_ep_pos(position_ep_x, position_ep_y, position_ep_z);
            //TVector3 position_n_pos(position_n_x, position_n_y, position_n_z);
            //mc_particle_distance = (position_n_pos - position_ep_pos).Mag(); // this is the first particle step, so is zero.

            TVector3 position_p1_pos(position_p1_x, position_p1_y, position_p1_z);
            TVector3 position_p2_pos(position_p2_x, position_p2_y, position_p2_z);
            ev_particle_distance = (position_p2_pos - position_p1_pos).Mag();
            
            // calculate time difference
            if(std::abs(nextDays - days) > 0){ //begin time tests
                coincidence_pass = false;
                time_ns_diff = std::abs(nextDays - days)*24*60*60*1e9 + std::abs(nextSec - secs)*1e9 + std::abs(nextNSec - nsecs);
            }
            else { // days are equal
                if(std::abs(nextSec - secs) > 0){
                    coincidence_pass = false;
                    time_ns_diff = std::abs(nextSec - secs)*1e9 + std::abs(nextNSec - nsecs);
                }
                else { // seconds are equal
                    time_ns_diff = std::abs(nextNSec - nsecs);
                    if(std::abs(nextNSec - nsecs) < deltaT) // nanosec are to within deltaT
                        coincidence_pass = true;
                    else
                        coincidence_pass = false; // so close..
                }
            }

            // conditions for event to pass cuts
            // both particles to be fitted within radius
            if ((position_p1_r < deltaP)&&(position_p2_r < deltaP))
                position_r_pass = true;
            // inter-particle distance cut
            if (ev_particle_distance < deltaD)
                particle_distance_pass = true;
            // time window cut
            time_window_pass = (coincidence_pass_previous || coincidence_pass);
            // energy cut
            if ((energy_ep>energy_ep_min)&&(energy_ep<energy_ep_max)&&(energy_n>energy_n_min)&&(energy_n<energy_n_max))
                energy_pass = true;
            
            all_pass = fit_validity && time_window_pass && position_r_pass && particle_distance_pass && energy_pass;

            //if passed then push events into vectors and mark event to be added to tree
            if(all_pass){
                if (j < (ev_nhit->size()-1)){
                    output_ev_fit_validity.push_back(fit_validity);
                    output_ev_fit_energy.push_back(ev_fit_energy->at(j));
                    output_ev_fit_position_r.push_back(ev_fit_position_r->at(j));
                    output_ev_fit_position_x.push_back(ev_fit_position_x->at(j));
                    output_ev_fit_position_y.push_back(ev_fit_position_y->at(j));
                    output_ev_fit_position_z.push_back(ev_fit_position_z->at(j));
                    output_ev_nhit.push_back(ev_nhit->at(j));
                    output_ev_time_days.push_back(ev_time_days->at(j));
                    output_ev_time_seconds.push_back(ev_time_seconds->at(j));
                    output_ev_time_nanoseconds.push_back(ev_time_nanoseconds->at(j));
                    ev_passed = true; // if any ev passes, throw a flag to later fill the tree
                }
            }

            //plotting events
            //fill hist with all events and then with those which passed the cuts
            h_before_cut_nhit.Fill(nhit_p1);
            h_before_cut_efit.Fill(energy_p1);
            h2_before_cut_nhit.Fill(nhit_p1, nhit_p2);
            h2_before_cut_emc.Fill(energy_ep, energy_n);
            h2_before_cut_efit.Fill(energy_p1, energy_p2);
            h_before_cut_time_diff.Fill(time_ns_diff);
            h_before_cut_position_displacement.Fill(ev_particle_distance);
            h2_before_cut_time_diff_displacement.Fill(time_ns_diff, ev_particle_distance);
            h2_before_cut_energy_resolution.Fill((energy_p1-energy_ep)/energy_ep, (energy_p2-energy_n)/energy_n);
            h2_before_cut_position_resolution.Fill((position_p1_x-position_ep_x)/position_ep_x, (position_p2_x-position_n_x)/position_n_x);
            
            if (all_pass){ // this individual passed
                h2_after_cut_emc.Fill(energy_ep, energy_n);
            }
            
            
            if(coincidencePartner==true)// only fill if they havent been plotted on the last loop cycle
                coincidencePartner=false; // (and if they were plotted last loop cycle, reset the flag for the next event)
            else{
                if (all_pass){ // this individual passed
                    coincidencePartner = true;

                    h_after_cut_nhit.Fill(nhit_p1);
                    h_after_cut_nhit.Fill(nhit_p2); // using coincidencePartner flag. fill both partners at same time.
                    h_after_cut_efit.Fill(energy_p1);
                    h_after_cut_efit.Fill(energy_p2); // using coincidencePartner flag. fill both partners at same time.
                    h2_after_cut_nhit.Fill(nhit_p1, nhit_p2);
                    h2_after_cut_efit.Fill(energy_p1, energy_p2);
                    h_after_cut_time_diff.Fill(time_ns_diff);
                    h_after_cut_position_displacement.Fill(ev_particle_distance);
                    h2_after_cut_time_diff_displacement.Fill(time_ns_diff, ev_particle_distance);
                    h2_after_cut_energy_resolution.Fill((energy_p1-energy_ep)/energy_ep, (energy_p2-energy_n)/energy_n);
                    h2_after_cut_position_resolution.Fill((position_p1_x-position_ep_x)/position_ep_x, (position_p2_x-position_n_x)/position_n_x);
                }
            }

            // shuffle the flag along, if this event is tagged by it's coincidence with the next (coincidence_pass)
            // the next one is flagged by its coincidence with the one before (coincidence_pass_previous)
            coincidence_pass_previous = coincidence_pass;
        }

        // if at least one ev passed, then fill the output tree
        if (ev_passed){
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

            tree_output->Fill();
            n_passed++;
        }
    }
    file_input->Close();
    file_output->cd();

    h_before_cut_nhit.GetXaxis()->SetTitle("nhit_p1");
    h_before_cut_efit.GetXaxis()->SetTitle("energy_p1");
    h2_before_cut_nhit.GetXaxis()->SetTitle("nhit_p1");
    h2_before_cut_nhit.GetYaxis()->SetTitle("nhit_p2");
    h2_before_cut_emc.GetXaxis()->SetTitle("energy_ep");
    h2_before_cut_emc.GetYaxis()->SetTitle("energy_n");
    h2_before_cut_efit.GetXaxis()->SetTitle("energy_p1");
    h2_before_cut_efit.GetYaxis()->SetTitle("energy_p2");
    h_before_cut_time_diff.GetXaxis()->SetTitle("time_ns_diff");
    h_before_cut_position_displacement.GetXaxis()->SetTitle("ev_particle_distance");
    h2_before_cut_time_diff_displacement.GetXaxis()->SetTitle("time_ns_diff");
    h2_before_cut_time_diff_displacement.GetYaxis()->SetTitle("ev_particle_distance");
    h2_before_cut_energy_resolution.GetXaxis()->SetTitle("(energy_p1-energy_ep)/energy_ep");
    h2_before_cut_energy_resolution.GetYaxis()->SetTitle("(energy_p2-energy_n)/energy_n");
    h2_before_cut_position_resolution.GetXaxis()->SetTitle("(position_p1_r-position_ep_r)/position_ep_r");
    h2_before_cut_position_resolution.GetYaxis()->SetTitle("(position_p2_r-position_n_r)/position_n_r");
    
    h_after_cut_nhit.GetXaxis()->SetTitle("nhit_p1");
    h_after_cut_efit.GetXaxis()->SetTitle("energy_p1");
    h2_after_cut_nhit.GetXaxis()->SetTitle("nhit_p1");
    h2_after_cut_nhit.GetYaxis()->SetTitle("nhit_p2");
    h2_after_cut_emc.GetXaxis()->SetTitle("energy_ep");
    h2_after_cut_emc.GetYaxis()->SetTitle("energy_n");
    h2_after_cut_efit.GetXaxis()->SetTitle("energy_p1");
    h2_after_cut_efit.GetYaxis()->SetTitle("energy_p2");
    h_after_cut_time_diff.GetXaxis()->SetTitle("time_ns_diff");
    h_after_cut_position_displacement.GetXaxis()->SetTitle("ev_particle_distance");
    h2_after_cut_time_diff_displacement.GetXaxis()->SetTitle("time_ns_diff");
    h2_after_cut_time_diff_displacement.GetYaxis()->SetTitle("ev_particle_distance");
    h2_after_cut_energy_resolution.GetXaxis()->SetTitle("(energy_p1-energy_ep)/energy_ep");
    h2_after_cut_energy_resolution.GetYaxis()->SetTitle("(energy_p2-energy_n)/energy_n");
    h2_after_cut_position_resolution.GetXaxis()->SetTitle("(position_p1_r-position_ep_r)/position_ep_r");
    h2_after_cut_position_resolution.GetYaxis()->SetTitle("(position_p2_r-position_n_r)/position_n_r");
    
    h_before_cut_nhit.Write();
    h_before_cut_efit.Write();
    h2_before_cut_nhit.Write();
    h2_before_cut_emc.Write();
    h2_before_cut_efit.Write();
    h_before_cut_time_diff.Write();
    h_before_cut_position_displacement.Write();
    h2_before_cut_time_diff_displacement.Write();
    h2_before_cut_energy_resolution.Write();
    h2_before_cut_position_resolution.Write();

    h_after_cut_nhit.Write();
    h_after_cut_efit.Write();
    h2_after_cut_nhit.Write();
    h2_after_cut_emc.Write();
    h2_after_cut_efit.Write();
    h_after_cut_time_diff.Write();
    h_after_cut_position_displacement.Write();
    h2_after_cut_time_diff_displacement.Write();
    h2_after_cut_energy_resolution.Write();
    h2_after_cut_position_resolution.Write();

    //tree_output->Print();
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
