#include <TFile.h>
#include <RAT/DB.hh>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <string>
#include <TTree.h>
#include <iostream>
#include <RAT/DU/Utility.hh>
#include <TObject.h>
#include <math.h>

void process_cuts(const std::string filename_input, const std::string filename_output, Double_t energy_ep_min, Double_t energy_ep_max, Double_t energy_n_min, Double_t energy_n_max, Double_t deltaT = 1000000, Double_t deltaP = 6000, Double_t deltaD = 6000){

    // load input file
    TFile *file_input = new TFile(filename_input.c_str());
    TTree *tree_input = (TTree*)file_input->Get("nt");

    // load neutron pdf
    TFile *file_input_n = new TFile("/data/snoplus/lidgard/OXSX/flux100/myNeutron_flux100.root");
    TH2D *h2_after_cut_fit_neutron = (TH2D*)file_input_n->Get("h2_after_cut_efit_neutron");

    // setup output file
    TFile *file_output = new TFile(filename_output.c_str(), "RECREATE");
    TTree *tree_output = tree_input->CloneTree(0);

    const ULong64_t n_slices = 100;
    char *name = new char[1000];
    TH1D *h_after_cut_fit_neutron[n_slices];
    //TF1 *f_after_cut_fit_neutron[n_slices];
    for (ULong64_t i=0; i<n_slices; i++){
        sprintf(name, "h_after_cut_fit_neutron%llu",i);
        h_after_cut_fit_neutron[i] = h2_after_cut_fit_neutron->ProjectionY(name,i*1-3,i*1+3);
        sprintf(name, "h_after_cut_fit_neutron%llu ",i);
        h_after_cut_fit_neutron[i]->SetTitle(name);
        //if (h_after_cut_fit_neutron[i]->GetEntries()>0){ // to use TF1 to fit slices (no better)
        //    sprintf(name, "f_after_cut_fit_neutron%d",i);
            //f_after_cut_fit_neutron[i] = new TF1(name,"gaus",0,3);
            //h_after_cut_fit_neutron[i]->Fit(f_after_cut_fit_neutron[i], "RQN");
        //}
    }

    TH2D h2_before_cut_emc("h2_before_cut_emc", "h2_before_cut_emc", 10000, 0, 10, 1000, 0, 1);
    TH1D h_before_cut_emc("h_before_cut_emc", "h_before_cut_emc", 1000, 0, 10);
    TH1D h_before_cut_emc_nu("h_before_cut_emc_nu", "h_before_cut_emc_nu", 1000, 0, 10);

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

    Double_t neutron_capture_energy = 1.857;
    Double_t e_rem = 0.784;

    ULong64_t n_passed = 0;
    ULong64_t n_parent_entries = 0;
    Double_t time_ns_diff;
    ULong64_t j;
    Bool_t all_pass, energy_pass, coincidence_pass, position_r_pass, particle_distance_pass, any_ev_passed;

    // properties to load from tree
    ULong64_t entry;
    Double_t ev_pos_r, ev_pos_x, ev_pos_y, ev_pos_z, ev_particle_distance;
    Double_t ev_pos_r_i, ev_pos_x_i, ev_pos_y_i, ev_pos_z_i;
    Double_t ev_pos_r_p1, ev_pos_x_p1, ev_pos_y_p1, ev_pos_z_p1;
    Double_t ev_pos_r_p2, ev_pos_x_p2, ev_pos_y_p2, ev_pos_z_p2;
    Double_t ev_energy, ev_energy_i, ev_energy_p1, ev_energy_corrected_p1, ev_energy_p2;
    UInt_t ev_time_seconds, ev_time_seconds_i, ev_time_days, ev_time_days_i, ev_time_days_p1, ev_time_seconds_p1, ev_time_days_p2, ev_time_seconds_p2;
    Double_t ev_time_nanoseconds, ev_time_nanoseconds_i, ev_time_nanoseconds_p1, ev_time_nanoseconds_p2;
    Int_t ev_nhit;//, ev_nhit_i;
    Bool_t ev_validity;
    ULong64_t ev_index;

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
    tree_input->SetBranchAddress("ev_index", &ev_index);
    tree_input->SetBranchAddress("ev_fit_validity", &ev_validity);
    tree_input->SetBranchAddress("ev_fit_energy", &ev_energy);
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

    std::vector<ULong64_t> tagged_entries;

    ULong64_t n_entries = tree_input->GetEntries();
    Int_t i_percent = (n_entries/100)*10;
    std::cout<<"initial entries: "<<n_entries<<std::endl;
    for (ULong64_t i=0; i < (n_entries-1); i++){ // since we're comparing partners we go to entry n-1

        // print progress
        if(i_percent && !(i%i_percent))
            std::cout << i/i_percent*10 << "% done " << std::endl;

        tree_input->GetEntry(i);

        if (ev_index==0) n_parent_entries++; //count parent entries (not second triggers)

        any_ev_passed = false;

        // // basic cuts to remove junk
        // if ((ev_validity==false)||(ev_nhit<=20)||(ev_pos_r>6000))
            // continue;

        // for events which pass most basic cut, compare events within 1s
        // first test all ev entries, then test between entries
        ev_time_days_i = ev_time_days;
        ev_time_seconds_i = ev_time_seconds;
        ev_time_nanoseconds_i = ev_time_nanoseconds;
        ev_pos_r_i = ev_pos_r;
        ev_pos_x_i = ev_pos_x;
        ev_pos_y_i = ev_pos_y;
        ev_pos_z_i = ev_pos_z;
        ev_energy_i = ev_energy;

        j = i; // go forward from i'th event
        while ((j <= (i+100)) && (j < (n_entries-1))){ // look through events (next 100)

            j++; // test the next trigger
            tree_input->GetEntry(j);

            // reset pass flags
            all_pass = false;
            coincidence_pass = false;
            position_r_pass = false;
            energy_pass = false;
            particle_distance_pass = false;

            // basic cuts to remove junk
            //if ((ev_validity==false)||(ev_nhit<=20)||(ev_pos_r>6000))
            //    continue;

            // get data from vectors for this trigger
            // p1 is first (prompt particle i.e. assumed positron)
            // check to see if nanoseconds are reversed (and switch if needed)
            if (ev_time_nanoseconds_i <= ev_time_nanoseconds){ //
                ev_time_days_p1 = ev_time_days_i;
                ev_time_seconds_p1 = ev_time_seconds_i;
                ev_time_nanoseconds_p1 = ev_time_nanoseconds_i;
                ev_pos_r_p1 = ev_pos_r_i;
                ev_pos_x_p1 = ev_pos_x_i;
                ev_pos_y_p1 = ev_pos_y_i;
                ev_pos_z_p1 = ev_pos_z_i;
                ev_energy_p1 = ev_energy_i;

                ev_time_days_p2 = ev_time_days;
                ev_time_seconds_p2 = ev_time_seconds;
                ev_time_nanoseconds_p2 = ev_time_nanoseconds;
                ev_pos_r_p2 = ev_pos_r;
                ev_pos_x_p2 = ev_pos_x;
                ev_pos_y_p2 = ev_pos_y;
                ev_pos_z_p2 = ev_pos_z;
                ev_energy_p2 = ev_energy;
            }
            else{
                ev_time_days_p1 = ev_time_days;
                ev_time_seconds_p1 = ev_time_seconds;
                ev_time_nanoseconds_p1 = ev_time_nanoseconds;
                ev_pos_r_p1 = ev_pos_r;
                ev_pos_x_p1 = ev_pos_x;
                ev_pos_y_p1 = ev_pos_y;
                ev_pos_z_p1 = ev_pos_z;
                ev_energy_p1 = ev_energy;

                ev_time_days_p2 = ev_time_days_i;
                ev_time_seconds_p2 = ev_time_seconds_i;
                ev_time_nanoseconds_p2 = ev_time_nanoseconds_i;
                ev_pos_r_p2 = ev_pos_r_i;
                ev_pos_x_p2 = ev_pos_x_i;
                ev_pos_y_p2 = ev_pos_y_i;
                ev_pos_z_p2 = ev_pos_z_i;
                ev_energy_p2 = ev_energy_i;
            }

            // find 'other' energy and correction
            Double_t e_n = 0;
            UInt_t my_n = round(ev_energy_p1*10);
            if (my_n>=n_slices) {
                //std::cout<< "warning: slice too high (ev_time_nanoseconds_p1 <= ev_time_nanoseconds_p2): " << my_n << " e_p1: " << ev_energy_p1 << std::endl;
                my_n = n_slices-1;
            }
            if (my_n<0) my_n = 0;
            if (h_after_cut_fit_neutron[my_n]->GetEntries()>0)
                e_n = h_after_cut_fit_neutron[my_n]->GetRandom();
            ev_energy_corrected_p1 = ev_energy_p1+e_rem+e_n;

            // time window cut
            // calculate time difference
            if(std::abs(ev_time_days_p2 - ev_time_days_p1) > 0){ //begin time tests
                coincidence_pass = false;
                time_ns_diff = std::abs(ev_time_days_p2 - ev_time_days_p1)*24*60*60*1e9 + std::abs(ev_time_seconds_p2 - ev_time_seconds_p1)*1e9 + std::abs(ev_time_nanoseconds_p2 - ev_time_nanoseconds_p1);
            }
            else { // days are equal
                if(std::abs(ev_time_seconds_p2 - ev_time_seconds_p1) > 0){
                    coincidence_pass = false;
                    time_ns_diff = std::abs(ev_time_seconds_p2 - ev_time_seconds_p1)*1e9 + std::abs(ev_time_nanoseconds_p2 - ev_time_nanoseconds_p1);
                }
                else { // seconds are equal
                    time_ns_diff = std::abs(ev_time_nanoseconds_p2 - ev_time_nanoseconds_p1);
                    if(std::abs(ev_time_nanoseconds_p2 - ev_time_nanoseconds_p1) < deltaT) // nanosec are to within deltaT
                        coincidence_pass = true;
                    else{
                        coincidence_pass = false; // so close..
                        //continue; // if this cut failed, don't go any further..
                    }
                }
            }

            // both particles to be fitted within radius
            if ((ev_pos_r_p1 < deltaP)&&(ev_pos_r_p2 < deltaP))
                position_r_pass = true;
            //else
            //    continue; // if this cut failed, don't go any further..

            // Inter-particle distance cut
            TVector3 position_p1_pos(ev_pos_x_p1, ev_pos_y_p1, ev_pos_z_p1);
            TVector3 position_p2_pos(ev_pos_x_p2, ev_pos_y_p2, ev_pos_z_p2);
            ev_particle_distance = (position_p2_pos - position_p1_pos).Mag();
            if (ev_particle_distance < deltaD)
                particle_distance_pass = true;
            //else
            //    continue; // if this cut failed, don't go any further..

            // energy cut
            if ((ev_energy_corrected_p1 > energy_ep_min) && (ev_energy_corrected_p1 < energy_ep_max) && (ev_energy_p2 > energy_n_min) && (ev_energy_p2 < energy_n_max))
                energy_pass = true;
            //else
            //    continue; // if this cut failed, don't go any further..

            // reached master cut flag
            all_pass = ev_validity && coincidence_pass && position_r_pass && particle_distance_pass && energy_pass;

            if (all_pass){
                // mark this i'th trigger as passed (it and at least 1 partner passed cuts)
                any_ev_passed = true;

                //tagged_entries.push_back(j);

                // fill histograms with all triggers which pass cuts
                h2_after_cut_efit.Fill(ev_energy_p1, ev_energy_p2);
                h2_after_cut_efit_corrected.Fill(ev_energy_corrected_p1, ev_energy_p2);
                h_after_cut_efit_p1_corrected.Fill(ev_energy_corrected_p1);
                h2_after_cut_efit_neutron.Fill(ev_energy_p1, mc_energy_nu - (ev_energy_p1+e_rem) );
                h_after_cut_efit_p2.Fill(ev_energy_p2);
                h_after_cut_time_diff.Fill(time_ns_diff);
                h_after_cut_position_displacement.Fill(ev_particle_distance);
                h2_after_cut_time_diff_displacement.Fill(time_ns_diff, ev_particle_distance);
                h2_after_cut_energy_resolution.Fill((ev_energy_corrected_p1-mc_energy_nu)/mc_energy_nu, (ev_energy_p2-neutron_capture_energy)/neutron_capture_energy);
                h2_after_cut_position_resolution.Fill((ev_pos_x_p1-mc_pos_x_ep)/mc_pos_x_ep, (ev_pos_x_p2-mc_pos_x_n)/mc_pos_x_n);
            }
        }
        // has this particle been tagged as a partner of another particle previously?
        //for (ULong64_t ii = 0; ii < tagged_entries.size(); ii++){
            //if (i == tagged_entries.at(ii))
            //    std::cout << "WARNING! already tagged i" << ii << " ev: " << i << std::endl;
        //}
        
        // load the original event (not the j'th)
        tree_input->GetEntry(i);

        // if at least one ev passed, then fill the output tree
        if (any_ev_passed){
            // fill tree
            tree_output->Fill();
            n_passed++;
        }

        // fill mc event histograms
        // fill with all events (ev_index == 0 and in this outer loop since there's only one parent per event)
        if (ev_index==0){
            h_before_cut_emc_nu.Fill(mc_energy_nu);
            h_before_cut_emc.Fill(mc_energy_ep);
            h2_before_cut_emc.Fill(mc_energy_ep, mc_energy_n);

            if (any_ev_passed){ // fill with only those which had a trigger which passed cuts
                h_after_cut_emc_nu.Fill(mc_energy_nu);
                h_after_cut_emc.Fill(mc_energy_ep);
                h2_after_cut_emc.Fill(mc_energy_ep, mc_energy_n);
            }
        }
    }

    // if (tagged_entries.size()>0){ // if the doubly-tagged entry vector has any entries, then print the vector.
        // for (ULong64_t ii = 0; ii < tagged_entries.size(); ii++)
            // std::cout << "i(" << ii << ") entry_i: " << tagged_entries.at(ii) << std::endl;
    // }

    //HISTOGRAMS
    // cut purity plots (ratio plot)
    TH1D *h_after_cut_emc_nu_ratio = (TH1D*)h_after_cut_emc_nu.Clone("h_after_cut_emc_nu_ratio");
    h_after_cut_emc_nu_ratio->Divide(&h_before_cut_emc_nu);

    //TH1D *h_before_cut_efit_ratio = (TH1D*)h_before_cut_efit_p1_corrected.Clone("h_before_cut_efit_ratio");
    //h_before_cut_efit_ratio->Divide(&h_before_cut_emc_nu);

    TH1D *h_after_cut_efit_ratio = (TH1D*)h_after_cut_efit_p1_corrected.Clone("h_after_cut_efit_ratio");
    h_after_cut_efit_ratio->Divide(&h_after_cut_emc_nu);

    TH1D *h_after_cut_emc_nu_n = (TH1D*)h_after_cut_emc_nu.Clone("h_after_cut_emc_nu_n");
    TH1D *h_after_cut_efit_p1_corrected_n = (TH1D*)h_after_cut_efit_p1_corrected.Clone("h_after_cut_efit_p1_corrected_n");
    h_after_cut_emc_nu_n->Scale(1./h_after_cut_emc_nu_n->GetMaximum());
    h_after_cut_efit_p1_corrected_n->Scale(1./h_after_cut_efit_p1_corrected_n->GetMaximum());
    h_after_cut_efit_p1_corrected.SetLineColor(kRed);
    h_after_cut_efit_p1_corrected_n->SetLineColor(kRed);

    // save file
    file_input->Close();
    file_input_n->Close();
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
    h2_after_cut_efit_corrected.GetXaxis()->SetTitle("ev_energy_p1+0.784MeV+E_n");
    h2_after_cut_efit_corrected.GetYaxis()->SetTitle("ev_energy_p2");
    h2_after_cut_efit_neutron.GetXaxis()->SetTitle("ev_energy_p1");
    h2_after_cut_efit_neutron.GetYaxis()->SetTitle("mc_energy_nu - (ev_energy_p1+e_rem)");
    h2_after_cut_time_diff_displacement.GetXaxis()->SetTitle("time_ns_diff");
    h2_after_cut_time_diff_displacement.GetYaxis()->SetTitle("ev_particle_distance");
    h2_after_cut_energy_resolution.GetXaxis()->SetTitle("(ev_energy_p1+0.784MeV+E_n-mc_energy_nu)/mc_energy_nu");
    h2_after_cut_energy_resolution.GetYaxis()->SetTitle("(ev_energy_p2+0.784MeV-n_capE)/n_capE");
    h2_after_cut_position_resolution.GetXaxis()->SetTitle("(ev_pos_x_p1-mc_pos_x_ep)/mc_pos_x_ep");
    h2_after_cut_position_resolution.GetYaxis()->SetTitle("(ev_pos_x_p2-mc_pos_x_n)/mc_pos_x_n");
    h_after_cut_emc.SetTitle("KE_ep");
    h_after_cut_emc_nu.SetTitle("KE_nu");
    h_after_cut_emc_nu_n->SetTitle("KE_nu (normalised (maximum = 1)");
    h_after_cut_efit_p1_corrected.SetTitle("ev_energy_p1+0.784MeV+E_n");
    h_after_cut_efit_p2.SetTitle("ev_energy_p2");
    h_after_cut_efit_p1_corrected_n->SetTitle("ev_energy_p1+0.784MeV+E_n");
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
    h_after_cut_efit_ratio->GetYaxis()->SetTitle("ev_energy_p1+0.784MeV+E_n / KE_nu");
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
    h_after_cut_efit_ratio->Write();

    //for (ULong64_t ii=0; ii<n_slices; ii++)
    //    h_after_cut_fit_neutron[ii]->Write();

    tree_output->AutoSave();
    file_output->Close();
    std::cout<<"number of (parent) entries: "<<n_parent_entries<<std::endl;
    std::cout<<"final entries: "<<n_passed<<std::endl;
}

Int_t main(Int_t argc, char *argv[]) {

    if (argc != 10) {
        std::cout<<"Error: 9 arguments expected. Got: "<<argc-1<<std::endl;
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

        process_cuts(filename_input, filename_output, energy_ep_min, energy_ep_max, energy_n_min, energy_n_max, deltaT, deltaP, deltaD);

        return 0; // completed successfully
    }
}
