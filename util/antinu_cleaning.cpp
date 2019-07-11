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

void process_cuts(const std::string filename_input, const std::string filename_output, Double_t energy_ep_min, Double_t energy_ep_max, Double_t energy_n_min, Double_t energy_n_max, Double_t deltaT = 1000000, Double_t deltaP = 6000, Double_t deltaD = 6000){

    // load input file
    TFile *file_input = new TFile(filename_input.c_str());
    TTree *tree_input = (TTree*)file_input->Get("nt");

    // setup output file
    TFile *file_output = new TFile(filename_output.c_str(), "RECREATE");
    TTree *tree_output = tree_input->CloneTree(0);

    char *name = new char[1000];

    TH2D h2_before_cut_emc("h2_before_cut_emc", "h2_before_cut_emc", 10000, 0, 10, 1000, 0, 1);
    TH1D h_before_cut_emc("h_before_cut_emc", "h_before_cut_emc", 100, 0, 10);
    TH1D h_before_cut_emc_nu("h_before_cut_emc_nu", "h_before_cut_emc_nu", 100, 0, 10);

    TH2D h2_after_cut_emc("h2_after_cut_emc", "h2_after_cut_emc", 10000, 0, 10, 1000, 0, 1);
    TH1D h_after_cut_emc("h_after_cut_emc", "h_after_cut_emc", 100, 0, 10);
    TH1D h_after_cut_emc_nu("h_after_cut_emc_nu", "h_after_cut_emc_nu", 100, 0, 10);
    TH2D h2_after_cut_efit("h2_after_cut_efit", "h2_after_cut_efit", 1000, 0, 10, 1000, 0, 10);
    TH1D h_after_cut_efit_p1("h_after_cut_efit_p1", "h_after_cut_efit_p1", 100, 0, 10);
    TH2D h2_after_cut_efit_neutron("h2_after_cut_efit_neutron", "h2_after_cut_efit_neutron", 100, 0, 10, 300, 0, 3);
    TH1D h_after_cut_efit_p2("h_after_cut_efit_p2", "h_after_cut_efit_p2", 100, 0, 10);
    TH1D h_after_cut_time_diff("h_after_cut_time_diff", "h_after_cut_time_diff", 1000, 0, 5000000);
    TH2D h2_after_cut_energy_resolution("h2_after_cut_energy_resolution", "h2_after_cut_energy_resolution", 200, -2, 2, 201, -2, 2);
    TH2D h2_after_cut_position_resolution("h2_after_cut_position_resolution", "h2_after_cut_position_resolution", 200, -2, 2, 201, -2, 2);
    TH1D h_after_cut_position_displacement("h_after_cut_position_displacement", "h_after_cut_position_displacement", 1000, 0, 10000);
    TH2D h2_after_cut_time_diff_displacement("h2_after_cut_time_diff_displacement", "h2_after_cut_time_diff_displacement", 1000, 0, 5000000, 1000, 0, 10000);

    Double_t neutron_capture_energy = 2.2;//1.857;
    Double_t e_rem = 0.784;

    ULong64_t n_passed = 0;
    ULong64_t already_tagged = 0;
    ULong64_t n_parent_entries = 0;
    ULong64_t n_parent_particles = 0;
    Double_t time_ns_diff;
    ULong64_t i_original, j;
    Bool_t all_pass, energy_pass, coincidence_pass, position_r_pass, particle_distance_pass, partner_passed;

    // properties to load from tree
    ULong64_t entry, entry_i;
    Double_t ev_pos_r, ev_pos_x, ev_pos_y, ev_pos_z, ev_particle_distance;
    Double_t ev_pos_r_i, ev_pos_x_i, ev_pos_y_i, ev_pos_z_i;
    Double_t ev_pos_r_p1, ev_pos_x_p1, ev_pos_y_p1, ev_pos_z_p1;
    Double_t ev_pos_r_p2, ev_pos_x_p2, ev_pos_y_p2, ev_pos_z_p2;
    Double_t ev_energy, ev_energy_i, ev_energy_p1, ev_energy_p2;
    UInt_t ev_time_seconds, ev_time_seconds_i, ev_time_days, ev_time_days_i, ev_time_days_p1, ev_time_seconds_p1, ev_time_days_p2, ev_time_seconds_p2;
    Double_t ev_time_nanoseconds, ev_time_nanoseconds_i, ev_time_nanoseconds_p1, ev_time_nanoseconds_p2;
    Int_t ev_nhit;//, ev_nhit_i;
    Bool_t ev_validity;
    ULong64_t ev_index, ev_n_index, mc_n_ev_index, mc_ev_index_ep, mc_ev_index_n, ev_index_p1, ev_index_p2;

    Double_t mc_pos_r_nu, mc_pos_x_nu, mc_pos_y_nu, mc_pos_z_nu;
    Double_t mc_pos_n_r, mc_pos_x_n, mc_pos_y_n, mc_pos_z_n;
    Double_t mc_pos_r_ep, mc_pos_x_ep, mc_pos_y_ep, mc_pos_z_ep;
    Double_t mc_quench_i, mc_energy_nu, mc_energy_n, mc_energy_ep;
    UInt_t mc_time_days, mc_time_seconds;
    Double_t mc_time_nanoseconds;
    Double_t latitude, longitude, altitude, distance;


    // set branches
    tree_input->SetBranchAddress("entry", &entry);
    tree_input->SetBranchAddress("mc_time_days", &mc_time_days);
    tree_input->SetBranchAddress("mc_time_seconds", &mc_time_seconds);
    tree_input->SetBranchAddress("mc_time_nanoseconds", &mc_time_nanoseconds);
    tree_input->SetBranchAddress("mc_quench", &mc_quench_i);
    tree_input->SetBranchAddress("mc_neutrino_energy", &mc_energy_nu);
    tree_input->SetBranchAddress("mc_positron_energy", &mc_energy_ep);
    tree_input->SetBranchAddress("mc_neutron_energy", &mc_energy_n);
    tree_input->SetBranchAddress("mc_neutrino_position_r", &mc_pos_r_nu);
    tree_input->SetBranchAddress("mc_neutrino_position_x", &mc_pos_x_nu);
    tree_input->SetBranchAddress("mc_neutrino_position_y", &mc_pos_y_nu);
    tree_input->SetBranchAddress("mc_neutrino_position_z", &mc_pos_z_nu);
    tree_input->SetBranchAddress("mc_positron_position_r", &mc_pos_r_ep);
    tree_input->SetBranchAddress("mc_positron_position_x", &mc_pos_x_ep);
    tree_input->SetBranchAddress("mc_positron_position_y", &mc_pos_y_ep);
    tree_input->SetBranchAddress("mc_positron_position_z", &mc_pos_z_ep);
    tree_input->SetBranchAddress("mc_neutron_position_r", &mc_pos_n_r);
    tree_input->SetBranchAddress("mc_neutron_position_x", &mc_pos_x_n);
    tree_input->SetBranchAddress("mc_neutron_position_y", &mc_pos_y_n);
    tree_input->SetBranchAddress("mc_neutron_position_z", &mc_pos_z_n);
    tree_input->SetBranchAddress("mc_n_ev_index", &mc_n_ev_index);
    tree_input->SetBranchAddress("mc_ev_index_ep", &mc_ev_index_ep);
    tree_input->SetBranchAddress("mc_ev_index_n", &mc_ev_index_n);
    tree_input->SetBranchAddress("ev_index", &ev_index);
    tree_input->SetBranchAddress("ev_n_index", &ev_n_index);
    tree_input->SetBranchAddress("ev_fit_energy", &ev_energy);
    tree_input->SetBranchAddress("ev_fit_validity", &ev_validity);
    tree_input->SetBranchAddress("ev_fit_position_r", &ev_pos_r);
    tree_input->SetBranchAddress("ev_fit_position_x", &ev_pos_x);
    tree_input->SetBranchAddress("ev_fit_position_y", &ev_pos_y);
    tree_input->SetBranchAddress("ev_fit_position_z", &ev_pos_z);
    tree_input->SetBranchAddress("ev_nhit", &ev_nhit);
    tree_input->SetBranchAddress("ev_time_days", &ev_time_days);
    tree_input->SetBranchAddress("ev_time_seconds", &ev_time_seconds);
    tree_input->SetBranchAddress("ev_time_nanoseconds", &ev_time_nanoseconds);
    tree_input->SetBranchAddress("reactor_info_latitude", &latitude);
    tree_input->SetBranchAddress("reactor_info_longitude", &longitude);
    tree_input->SetBranchAddress("reactor_info_altitude", &altitude);
    tree_input->SetBranchAddress("reactor_info_distance", &distance);

    // values to modify
    tree_output->SetBranchAddress("ev_fit_energy_p1", &ev_energy_p1);
    tree_output->SetBranchAddress("ev_fit_energy_p2", &ev_energy_p2);
    tree_output->SetBranchAddress("ev_index_p1", &ev_index_p1);
    tree_output->SetBranchAddress("ev_index_p2", &ev_index_p2);

    std::vector<TString> tagged_entries;

    ULong64_t n_entries = tree_input->GetEntries();
    ULong64_t percent_interval = n_entries/10; //print every 10%
    if (percent_interval<=0) percent_interval = 1;
    ULong64_t progress_countdown = percent_interval;
    ULong64_t i = 0;
    ULong64_t p2 = 0;
    n_parent_particles = tree_input->GetEntries();

    std::vector<ULong64_t> entries_parents;
    std::vector<ULong64_t> entries_neutrons;
    std::vector<ULong64_t> entries_positrons;
    std::vector<ULong64_t> entries_neutrons2;
    std::vector<ULong64_t> entries_positrons2;
    std::vector<ULong64_t> entries_basic;
    std::vector<ULong64_t> entries_fiducial;
    std::vector<ULong64_t> entries_coincidence_ep;
    std::vector<ULong64_t> entries_coincidence_n;
    std::vector<ULong64_t> entries_distance_ep;
    std::vector<ULong64_t> entries_distance_n;
    std::vector<ULong64_t> entries_passed;
    std::vector<ULong64_t> entries_passed_2;
    std::vector<ULong64_t> entries_passed_2_ep;
    if (n_entries>1){

        // basic cut
        for (ULong64_t ii=0; ii < n_entries; ii++){
            i = ii;
            tree_input->GetEntry(i);

            entries_parents.push_back(entry); // entry for all events
            //std::cout<<"ii: "<<ii<<" entry: "<<entry<<" ev_index: "<<ev_index<<"/"<<ev_n_index<<" ev_energy: "<<ev_energy<<std::endl;

            if ((ev_validity==true)&&(ev_nhit>10)&&(ev_pos_r<=6000)) // basic nhit and position cuts to remove junk
                entries_basic.push_back(i);

            h_before_cut_emc_nu.Fill(mc_energy_nu);
            h_before_cut_emc.Fill(mc_energy_ep);
            h2_before_cut_emc.Fill(mc_energy_ep, mc_energy_n);
        }
        // vector method to get unique entries
        sort( entries_parents.begin(), entries_parents.end() );
        entries_parents.erase( unique( entries_parents.begin(), entries_parents.end() ), entries_parents.end() );
        n_parent_entries = entries_parents.size();
        std::cout<<"parent entries: "<<n_parent_entries<<std::endl;
        std::cout<<"initial particles: "<<n_entries<<std::endl;

        n_entries = entries_basic.size();
        std::cout<<"particles after basic cuts: "<<n_entries<<std::endl;

        // deltaP cut
        for (ULong64_t ii=0; ii < n_entries; ii++){
            i = entries_basic.at(ii);
            tree_input->GetEntry(i);
            if (ev_pos_r <= deltaP) // cut events outside of specified fiducial radius
                entries_fiducial.push_back(i);
        }
        n_entries = entries_fiducial.size();
        std::cout<<"particles after fiducial cut: "<<n_entries<<std::endl;

        // energy cut ** these are fiducial neutron & positron candidates
        for (ULong64_t ii=0; ii < n_entries; ii++){
            i = entries_fiducial.at(ii);
            tree_input->GetEntry(i);
            if ((ev_energy > energy_n_min)&&(ev_energy < energy_n_max)) // cut arround 2.2. MeV gamma of neutron capture
                entries_neutrons.push_back(i);
            if ((ev_energy > energy_ep_min)&&(ev_energy < energy_ep_max)) // cut positron energy spectrum
                entries_positrons.push_back(i);
        }
        std::cout<<"particles after neutron energy cut: "<<entries_neutrons.size()<<std::endl;
        std::cout<<"particles after positron energy cut: "<<entries_positrons.size()<<std::endl;
        n_entries = entries_neutrons.size();
        progress_countdown = percent_interval;

        // ev cut
        for (ULong64_t ii=0; ii < entries_neutrons.size(); ii++){
            i = entries_neutrons.at(ii);
            
            // print progress
            progress_countdown--;
            if (progress_countdown==0){
                progress_countdown = percent_interval;
                printf("ev cut %.0f%% done\n",(Double_t)(ii+1)/n_entries*100);
            }
            
            tree_input->GetEntry(i);
            entry_i = entry;
            ev_index_p2 = 0; // reset to 0 i.e positron like.
            bool passed = false;
            if (ev_index > 0) {// neutron can't be the first fitted particle (i.e. ev=0,1
                ev_index_p2 = ev_index;
                ev_energy_p2 = ev_energy;
                //std::cout<<"n: ev_index_p2: "<<ev_index<<"/"<<ev_n_index<<" entry: "<<entry<<" ev_energy: "<<ev_energy<<std::endl;
                passed = true;
            }

            if (passed){
                // now go through neutron events, find neutron event with same event number, select neutron for ev's > than this
                ULong64_t n_pairs = 0;
                for (ULong64_t iii = 0; iii < entries_positrons.size(); iii++){
                    j = entries_positrons.at(iii);
                    tree_input->GetEntry(j);
                    if ((entry_i == entry)&&(ev_index < ev_index_p2)&&(ev_index < 2)){ // positron must be first of 2 fitted particles
                        if ((int)(ev_n_index-ev_index)>1){ // positron must not be the last fitted particle i.e if n_ev=2 positron can't be ev=1, must be ev=0.
                            n_pairs ++;
                            ev_index_p1 = ev_index;
                            ev_energy_p1 = ev_energy;
                            // pass pairs of particles at a time - only pass neutron if positron was found.
                            entries_neutrons2.push_back(i); // but this is how we get duplicates in the neutron list.
                            entries_positrons2.push_back(j);
                            //std::cout<<"e: ev_index_p1: "<<ev_index_p1<<"/"<<ev_n_index<<" entry: "<<entry<<" ev_energy: "<<ev_energy<<" n_pairs: "<<n_pairs<<std::endl;
                            //std::cout<<"n: ev_index_p2: "<<ev_index_p2<<"/"<<ev_n_index<<" entry: "<<entry<<" ev_energy: "<<ev_energy_p2<<" n_pairs: "<<n_pairs<<std::endl;
                            if (ev_index == 1) std::cout<<"n: ev_index_p2: "<<ev_index_p2<<"/"<<ev_n_index<<" entry: "<<entry<<" ev_energy: "<<ev_energy_p2<<" n_pairs: "<<n_pairs<<std::endl;
                            
                        }
                    }
                }
                if (n_pairs>1)
                    std::cout<<"\n***\n***\n***** warning n_pairs: *****\n***\n***\n"<<n_pairs<<std::endl;
            }
        }
        std::cout<<"neutron after ev cut: "<<entries_neutrons2.size()<<std::endl;
        std::cout<<"positron after ev cut: "<<entries_positrons2.size()<<std::endl;
        n_entries = entries_neutrons2.size();
        progress_countdown = percent_interval;


        // time window cut
        for (ULong64_t ii=0; ii < n_entries; ii++){
            i = entries_neutrons2.at(ii);
            // print progress
            progress_countdown--;
            if (progress_countdown==0){
                progress_countdown = percent_interval;
                printf("coincident %.0f%% done\n",(Double_t)(ii+1)/n_entries*100);
            }
            tree_input->GetEntry(i);
            entry_i = entry;
            ev_time_days_p2 = ev_time_days;
            ev_time_seconds_p2 = ev_time_seconds;
            ev_time_nanoseconds_p2 = ev_time_nanoseconds;

            // find closest fiducial positron
            ULong64_t jj=0;
            for (jj = 0; jj < entries_positrons2.size(); jj++){
                j = entries_positrons2.at(jj);
                if (i <= j) break; //jj is now the index of the positron with the same particle number as the neutron
            }
            if (jj > entries_positrons2.size()) jj = entries_positrons2.size();
            //std::cout<<"\ti:"<<i<<" ii:"<<ii<<" jj:"<<entries_positrons2.at(jj)<<std::endl;

            // go backward through events to find positron, but ensure we go backwards without going beyond 0.
            ULong64_t jj_limit = jj-5; //set limit since jj will change value
            while ((jj > jj_limit)&&(jj > 0)&&(jj < entries_positrons2.size())){
                j = entries_positrons2.at(jj);
                std::cout<<"jj"<<jj<<" j:"<<j<<std::endl;

                tree_input->GetEntry(j);

                if ((i!=j)&&(j<i)){ //don't tag particle with itself and positron must be before neutron
                    if (entry_i==entry){ // for MC generated data, we know the positron and neutron are in the same event number **not true for data**

                        if(std::abs(ev_time_days_p2 - ev_time_days) > 0){ //begin time tests
                            time_ns_diff = std::abs(ev_time_days_p2 - ev_time_days)*24*60*60*1e9 + std::abs(ev_time_seconds_p2 - ev_time_seconds)*1e9 + std::abs(ev_time_nanoseconds_p2 - ev_time_nanoseconds);
                        }
                        else { // days are equal
                            if(std::abs(ev_time_seconds_p2 - ev_time_seconds) > 0){
                                time_ns_diff = std::abs(ev_time_seconds_p2 - ev_time_seconds)*1e9 + std::abs(ev_time_nanoseconds_p2 - ev_time_nanoseconds);
                            }
                            else { // seconds are equal
                                time_ns_diff = std::abs(ev_time_nanoseconds_p2 - ev_time_nanoseconds);
                                if(std::abs(ev_time_nanoseconds_p2 - ev_time_nanoseconds) < deltaT){ // nanosec are to within deltaT
                                    entries_coincidence_n.push_back(i);
                                    entries_coincidence_ep.push_back(j);
                                }
                            }
                        }
                    }
                }
                jj--;
            }
        }
        n_entries = entries_coincidence_n.size();
        std::cout<<"particles after coincident cut: "<<entries_coincidence_n.size()+entries_coincidence_ep.size()<<" (n:"<<entries_coincidence_n.size()<<", ep:"<<entries_coincidence_ep.size()<<")"<<std::endl;

        progress_countdown = percent_interval;
        // Inter-particle distance cut
        for (ULong64_t ii=0; ii < n_entries; ii++){
            i = entries_coincidence_n.at(ii);
            
            // print progress
            progress_countdown--;
            if (progress_countdown==0){
                progress_countdown = percent_interval;
                printf("distance %.0f%% done\n",(Double_t)(ii+1)/n_entries*100);
            }
            tree_input->GetEntry(i);
            TVector3 position_p2_pos(ev_pos_x, ev_pos_y, ev_pos_z);

            j = entries_coincidence_ep.at(ii);
            tree_input->GetEntry(j);
            TVector3 position_p1_pos(ev_pos_x, ev_pos_y, ev_pos_z);
            //TVector3 position_p1_pos(mc_pos_x_ep, mc_pos_y_ep, mc_pos_z_ep);

            ev_particle_distance = (position_p2_pos - position_p1_pos).Mag();
            if ((ev_particle_distance >= 1)&&(ev_particle_distance < deltaD)){
                entries_distance_n.push_back(i);
                entries_distance_ep.push_back(j);
                entries_passed.push_back(i);
                entries_passed.push_back(j);
            }
        }
        sort( entries_passed.begin(), entries_passed.end() );
        n_entries = entries_passed.size();
        std::cout<<"particles after inter-particle distance cut: "<<n_entries<<" (n:"<<entries_distance_n.size()<<", ep:"<<entries_distance_ep.size()<<")"<<std::endl;

        std::vector<ULong64_t>::iterator it1 = std::unique( entries_distance_n.begin(), entries_distance_n.end() );
        bool wasUnique_n = (it1 == entries_distance_n.end() );
        std::vector<ULong64_t>::iterator it2 = std::unique( entries_distance_ep.begin(), entries_distance_ep.end() );
        bool wasUnique_ep = (it2 == entries_distance_ep.end() );
        std::cout<<"uniqness n:"<<wasUnique_n<<" ep:"<<wasUnique_ep<<std::endl;

        // get final entry numbers
        for (ULong64_t ii=0; ii < entries_distance_n.size(); ii++){
            i = entries_distance_n.at(ii);
            tree_input->GetEntry(i);
            entries_passed_2.push_back(entry);

        }
        for (ULong64_t ii=0; ii < entries_distance_ep.size(); ii++){
            i = entries_distance_ep.at(ii);
            tree_input->GetEntry(i);
            entries_passed_2_ep.push_back(entry);
        }

        // for (ULong64_t ii=0; ii < entries_distance_n.size(); ii++){
            // bool same = 0;
            // if (entries_passed_2.at(ii)==entries_passed_2_ep.at(ii)) same=true;
            // std::cout<<"i:"<<ii<<" particle_j:"<<entries_distance_ep.at(ii)<<" particle_i:"<<entries_distance_n.at(ii)<<" entry_j:"<<entries_passed_2_ep.at(ii)<<" entry_i:"<<entries_passed_2.at(ii)<<":"<<same<<std::endl;
        // }

        // fill events which passed
        n_passed = n_entries;
        for (ULong64_t ii=0; ii < entries_distance_n.size(); ii++){
            i = entries_distance_n.at(ii);
            tree_input->GetEntry(i);
            ev_energy_p2 = ev_energy;
            ev_index_p2 = ev_index;
            ev_time_nanoseconds_p2 = ev_time_nanoseconds;
            ev_pos_x_p2 = ev_pos_x;
            ev_pos_y_p2 = ev_pos_y;
            ev_pos_z_p2 = ev_pos_z;
            TVector3 position_p2_pos(ev_pos_x_p2, ev_pos_y_p2, ev_pos_z_p2);

            j = entries_distance_ep.at(ii);
            tree_input->GetEntry(j);
            ev_energy_p1 = ev_energy;
            ev_index_p1 = ev_index;

            time_ns_diff = std::abs(ev_time_nanoseconds_p2 - ev_time_nanoseconds);
            TVector3 position_p1_pos(ev_pos_x, ev_pos_y, ev_pos_z);
            ev_particle_distance = (position_p2_pos - position_p1_pos).Mag();

            // // fill histograms with all triggers which pass cuts
            h2_after_cut_efit.Fill(ev_energy_p1, ev_energy_p2);
            h_after_cut_efit_p1.Fill(ev_energy_p1);
            h2_after_cut_efit_neutron.Fill(ev_energy_p1, mc_energy_nu - (ev_energy_p1+e_rem) );
            h_after_cut_efit_p2.Fill(ev_energy_p2);
            h_after_cut_time_diff.Fill(time_ns_diff);
            h_after_cut_position_displacement.Fill(ev_particle_distance);
            h2_after_cut_time_diff_displacement.Fill(time_ns_diff, ev_particle_distance);
            h2_after_cut_energy_resolution.Fill((ev_energy_p1-mc_energy_nu)/mc_energy_nu, (ev_energy_p2-neutron_capture_energy)/neutron_capture_energy);
            h2_after_cut_position_resolution.Fill((ev_pos_x-mc_pos_x_ep)/mc_pos_x_ep, (ev_pos_x_p2-mc_pos_x_n)/mc_pos_x_n);

            h_after_cut_emc_nu.Fill(mc_energy_nu);
            h_after_cut_emc.Fill(mc_energy_ep);
            h2_after_cut_emc.Fill(mc_energy_ep, mc_energy_n);

            tree_output->Fill();
        }

        std::vector<ULong64_t>::iterator it3 = std::unique( entries_passed.begin(), entries_passed.end() );
        bool wasUnique = (it3 == entries_passed.end() );
        std::cout<<"uniqness n:"<<wasUnique<<std::endl;

        p2 = entries_passed.size();

    }
    std::cout<<"passed entries: "<<p2<<std::endl;

    // cut purity plots (ratio plot)
    TH1D *h_after_cut_emc_nu_ratio = (TH1D*)h_after_cut_emc_nu.Clone("h_after_cut_emc_nu_ratio");
    h_after_cut_emc_nu_ratio->Divide(&h_before_cut_emc_nu);

    TH1D *h_after_cut_efit_ratio = (TH1D*)h_after_cut_efit_p1.Clone("h_after_cut_efit_ratio");
    h_after_cut_efit_ratio->Divide(&h_after_cut_emc_nu);

    TH1D *h_after_cut_emc_nu_n = (TH1D*)h_after_cut_emc_nu.Clone("h_after_cut_emc_nu_n");
    TH1D *h_after_cut_efit_p1_n = (TH1D*)h_after_cut_efit_p1.Clone("h_after_cut_efit_p1_n");
    h_after_cut_emc_nu_n->Scale(1./h_after_cut_emc_nu_n->GetMaximum());
    h_after_cut_efit_p1_n->SetLineColor(kRed);
    h_after_cut_efit_p1_n->Scale(1./h_after_cut_efit_p1_n->GetMaximum());

    // save file
    file_input->Close();
    file_output->cd();

    // set axes...
    h2_before_cut_emc.GetXaxis()->SetTitle("mc_energy_ep");
    h2_before_cut_emc.GetYaxis()->SetTitle("mc_energy_n");
    h_before_cut_emc.SetTitle("KE_ep");
    h_before_cut_emc_nu.SetTitle("KE_nu");
    h_before_cut_emc.GetXaxis()->SetTitle("Energy (MeV)");
    h_before_cut_emc_nu.GetXaxis()->SetTitle("Energy (MeV)");

    h2_after_cut_emc.GetXaxis()->SetTitle("mc_energy_ep");
    h2_after_cut_emc.GetYaxis()->SetTitle("mc_energy_n");
    h2_after_cut_efit.GetXaxis()->SetTitle("ev_energy_p1");
    h2_after_cut_efit.GetYaxis()->SetTitle("ev_energy_p2");
    h2_after_cut_efit_neutron.GetXaxis()->SetTitle("ev_energy_p1");
    h2_after_cut_efit_neutron.GetYaxis()->SetTitle("mc_energy_nu - (ev_energy_p1+e_rem)");
    h2_after_cut_time_diff_displacement.GetXaxis()->SetTitle("time_ns_diff");
    h2_after_cut_time_diff_displacement.GetYaxis()->SetTitle("ev_particle_distance");
    h2_after_cut_energy_resolution.GetXaxis()->SetTitle("(ev_energy_p1+0.784MeV-mc_energy_nu)/mc_energy_nu");
    h2_after_cut_energy_resolution.GetYaxis()->SetTitle("(ev_energy_p2+0.784MeV-n_capE)/n_capE");
    h2_after_cut_position_resolution.GetXaxis()->SetTitle("(ev_pos_x_p1-mc_pos_x_ep)/mc_pos_x_ep");
    h2_after_cut_position_resolution.GetYaxis()->SetTitle("(ev_pos_x_p2-mc_pos_x_n)/mc_pos_x_n");
    h_after_cut_emc.SetTitle("KE_ep");
    h_after_cut_emc_nu.SetTitle("KE_nu");
    h_after_cut_emc_nu_n->SetTitle("KE_nu (normalised (maximum = 1)");
    h_after_cut_efit_p1.SetTitle("ev_energy_p1");
    h_after_cut_efit_p1_n->SetTitle("ev_energy_p1_n");
    h_after_cut_efit_p2.SetTitle("ev_energy_p2");
    h_after_cut_time_diff.SetTitle("time_ns_diff");
    h_after_cut_position_displacement.SetTitle("ev_particle_distance");
    h_after_cut_emc.GetXaxis()->SetTitle("Energy (MeV)");
    h_after_cut_emc_nu.GetXaxis()->SetTitle("Energy (MeV)");
    h_after_cut_emc_nu_n->GetXaxis()->SetTitle("Energy (MeV)");
    h_after_cut_efit_p1.GetXaxis()->SetTitle("Energy (MeV)");
    h_after_cut_efit_p1_n->GetXaxis()->SetTitle("Energy (MeV)");
    h_after_cut_efit_p2.GetXaxis()->SetTitle("Energy (MeV)");
    h_after_cut_time_diff.GetXaxis()->SetTitle("Time between abs(ev_fit_p2 - ev_fit_p1) (ns)");
    h_after_cut_position_displacement.GetXaxis()->SetTitle("Position (ev_fit_p1 - ev_fit_p2) (mm)");

    h_after_cut_emc_nu_ratio->GetYaxis()->SetTitle("KE_nu after cut / KE_nu before cut");
    h_after_cut_efit_ratio->GetYaxis()->SetTitle("ev_energy_p1+0.784MeV / KE_nu");
    h_after_cut_emc_nu_ratio->GetXaxis()->SetTitle("Energy (MeV)");
    h_after_cut_efit_ratio->GetXaxis()->SetTitle("Energy (MeV)");

    // write objects...
    h2_before_cut_emc.Write();
    h_before_cut_emc_nu.Write();
    h_before_cut_emc.Write();

    h2_after_cut_emc.Write();
    h_after_cut_emc_nu.Write();
    h_after_cut_emc_nu_n->Write();
    h_after_cut_emc.Write();
    h2_after_cut_efit.Write();
    h_after_cut_efit_p1.Write();
    h_after_cut_efit_p1_n->Write();
    h2_after_cut_efit_neutron.Write();
    h_after_cut_efit_p2.Write();
    h_after_cut_time_diff.Write();
    h_after_cut_position_displacement.Write();
    h2_after_cut_time_diff_displacement.Write();
    h2_after_cut_energy_resolution.Write();
    h2_after_cut_position_resolution.Write();

    h_after_cut_emc_nu_ratio->Write();
    h_after_cut_efit_ratio->Write();

    tree_output->AutoSave();
    file_output->Close();
    std::cout<<"Summary: "<<std::endl;
    std::cout<<"\tparent entries: "<<n_parent_entries<<std::endl;
    std::cout<<"\tfinal entries: "<<p2<<std::endl;
    std::cout<<"\tpass percent: "<<(double)p2/(double_t)n_parent_entries*100<<std::endl;


    std::cout<<"n_min: "<<energy_n_min<<" n_max: "<<energy_n_max<<std::endl;

    //write csv output file with event numbers
    sprintf(name, "%s.csv",filename_output.c_str());
    FILE *fOut = fopen(name,"w");
    fprintf(fOut,"filename_input,filename_output,energy_ep_min,energy_ep_max,energy_n_min,energy_n_max,deltaT,deltaP,deltaD,initial_entries,final_entries,already_tagged,finished\n");
    fprintf(fOut,"%s,%s,%f,%f,%f,%f,%f,%f,%f,%i,%llu,%llu,%llu,%i\n", filename_input.c_str(), filename_output.c_str(), energy_ep_min, energy_ep_max, energy_n_min, energy_n_max, deltaT, deltaP, deltaD, n_parent_entries, n_passed, already_tagged,1);
    fclose(fOut);
}

Int_t main(Int_t argc, char *argv[]) {

    if (argc != 10) {
        std::cout<<"Error: 9 arguments expected. Got: "<<argc-1<<std::endl;
        return 1; // return>0 indicates error code
    }
    else {
        //TH1::AddDirectory(kFALSE);
        const std::string &filename_input = argv[1];
        const std::string &filename_output = argv[2];
        double energy_ep_min = atof(argv[3]);
        double energy_ep_max = atof(argv[4]);
        double energy_n_min = atof(argv[5]);
        double energy_n_max = atof(argv[6]);
        double deltaT = atof(argv[7]);
        double deltaP = atof(argv[8]);
        double deltaD = atof(argv[9]);

        //write csv output file to show process has begun (values filled upon completion)
        char *name = new char[1000];
        sprintf(name, "%s.csv",filename_output.c_str());
        FILE *fOut = fopen(name,"w");
        fprintf(fOut,"filename_input,filename_output,energy_ep_min,energy_ep_max,energy_n_min,energy_n_max,deltaT,deltaP,deltaD,initial_entries,final_entries,already_tagged,finished\n");
        fprintf(fOut,"%s,%s,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i\n", filename_input.c_str(), filename_output.c_str(), -9000, -9000, -9000, -9000, -9000, -9000, -9000, -9000, -9000, -9000, -9000,0);
        fclose(fOut);

        process_cuts(filename_input, filename_output, energy_ep_min, energy_ep_max, energy_n_min, energy_n_max, deltaT, deltaP, deltaD);

        return 0; // completed successfully
    }
}
