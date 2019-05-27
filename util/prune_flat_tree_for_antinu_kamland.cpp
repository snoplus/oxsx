#include <TFile.h>
#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/FitResult.hh>
#include <RAT/DB.hh>
#include <TMath.h>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <TTree.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <RAT/DU/Utility.hh>
#include <RAT/DB.hh>
//#include <TObject.h>

//const TVector3 SNO_ECEF_coord_ = TVector3(672.87, -4347.18, 4600.51);
const TVector3 SNO_ECEF_coord_ = TVector3(-3777.14425893, 3483.58137383, 3766.0181443); // Kamland

Double_t CalculateDistance(TVector3 point1, TVector3 point2) {
    return (point2 - point1).Mag();
}

TVector3 LLAtoECEF(Double_t latitude, Double_t longitude, Double_t altitude) {
    // reference http://www.mathworks.co.uk/help/aeroblks/llatoecefposition.html
    const Double_t toRad = TMath::Pi()/180.;
    const Double_t Earthradius = 6378137.0; //Radius of the Earth (in meters)
    const Double_t f_f = 1./298.257223563; //Flattening factor WGS84 Model
    Double_t x, y, z;
    latitude = latitude*toRad;
    longitude = longitude*toRad;
    const Double_t f_l = TMath::ATan(TMath::Power((1. - f_f),2)*TMath::Tan(latitude))/toRad;
    const Double_t f_rs = TMath::Sqrt( TMath::Power(Earthradius,2)/(1. + (1./TMath::Power((1. - f_f),2) - 1.)*TMath::Power(TMath::Sin(f_l*toRad),2)));
    x = (f_rs*TMath::Cos(f_l*toRad)*TMath::Cos(longitude) + altitude*TMath::Cos(latitude)*TMath::Cos(longitude))/1000.; // in km
    y = (f_rs*TMath::Cos(f_l*toRad)*TMath::Sin(longitude) + altitude*TMath::Cos(latitude)*TMath::Sin(longitude))/1000.; // in km
    z = (f_rs*TMath::Sin(f_l*toRad) + altitude*TMath::Sin(latitude))/1000.; // in km
    TVector3 ECEF(x,y,z);
    //std::cout << " LLAtoECEF: lat:" << latitude/toRad << " long:" << longitude/toRad << " alt:" << altitude << std::endl;
    //std::cout << " LLAtoECEF: f_f:" << f_f << " f_l:" << f_l << " f_s:" << f_rs/10000000. << std::endl;
    //std::cout << " LLAtoECEF: x:" << x << " y:" << y << " z:" << z << std::endl;
    return ECEF;
}

Double_t GetReactorDistanceLLA(Double_t latitude, Double_t longitude, Double_t altitude) {
    return CalculateDistance(SNO_ECEF_coord_,LLAtoECEF(latitude, longitude, altitude));
}

void ntload(std::string input_filename, std::string output_filename, std::string fitter_name) {

    TTree *tree_output = new TTree("nt","Anti-neutrino processed tree");
    ULong64_t n_total_events = 0;
    std::string reactorName="";
    std::string reactorcoreStr="";
    RAT::DB *db = RAT::DB::Get();
    RAT::DBLinkPtr fLink;
    std::vector<Double_t> latitude;
    std::vector<Double_t> longitude;
    std::vector<Double_t> altitude;

    // variables for branches
    ULong64_t entry;
    Double_t mc_energy_nu, mc_pos_x_nu, mc_pos_y_nu, mc_pos_z_nu, mc_pos_r_nu, mc_quench_i;
    Double_t mc_energy_n, mc_pos_x_n, mc_pos_y_n, mc_pos_z_n, mc_pos_r_n;
    Double_t mc_energy_ep, mc_pos_x_ep, mc_pos_y_ep, mc_pos_z_ep, mc_pos_r_ep;
    Double_t latitude_i, longitude_i, altitude_i;
    Double_t distance_i;
    UInt_t mc_time_days, mc_time_seconds;
    Double_t mc_time_nanoseconds;
    Int_t ev_nhit;
    UInt_t ev_time_days;
    Double_t ev_time_nanoseconds;
    UInt_t ev_time_seconds;
    Double_t ev_energy, ev_pos_x, ev_pos_y, ev_pos_z, ev_pos_r, ev_energy_p1, ev_energy_p2;
    Bool_t ev_validity;
    ULong64_t ev_index, mc_ev_index_ep, mc_ev_index_n, ev_index_p1, ev_index_p2;

    // set branches
    tree_output->Branch("entry", &entry, "entry/l");
    tree_output->Branch("mc_time_days", &mc_time_days, "mc_time_days/i");
    tree_output->Branch("mc_time_seconds", &mc_time_seconds, "mc_time_seconds/i");
    tree_output->Branch("mc_time_nanoseconds", &mc_time_nanoseconds, "mc_time_nanoseconds/D");
    tree_output->Branch("mc_quench", &mc_quench_i, "mc_quench_i/D");
    tree_output->Branch("mc_neutrino_energy", &mc_energy_nu, "mc_energy_nu/D");
    tree_output->Branch("mc_positron_energy", &mc_energy_ep, "mc_energy_ep/D");
    tree_output->Branch("mc_neutron_energy", &mc_energy_n, "mc_energy_n/D");
    tree_output->Branch("mc_neutrino_position_r", &mc_pos_r_nu, "mc_pos_r_nu/D");
    tree_output->Branch("mc_neutrino_position_x", &mc_pos_x_nu, "mc_pos_x_nu/D");
    tree_output->Branch("mc_neutrino_position_y", &mc_pos_y_nu, "mc_pos_y_nu/D");
    tree_output->Branch("mc_neutrino_position_z", &mc_pos_z_nu, "mc_pos_z_nu/D");
    tree_output->Branch("mc_positron_position_r", &mc_pos_r_ep, "mc_pos_r_ep/D");
    tree_output->Branch("mc_positron_position_x", &mc_pos_x_ep, "mc_pos_x_ep/D");
    tree_output->Branch("mc_positron_position_y", &mc_pos_y_ep, "mc_pos_y_ep/D");
    tree_output->Branch("mc_positron_position_z", &mc_pos_z_ep, "mc_pos_z_ep/D");
    tree_output->Branch("mc_neutron_position_r", &mc_pos_r_n, "mc_pos_r_n/D");
    tree_output->Branch("mc_neutron_position_x", &mc_pos_x_n, "mc_pos_x_n/D");
    tree_output->Branch("mc_neutron_position_y", &mc_pos_y_n, "mc_pos_y_n/D");
    tree_output->Branch("mc_neutron_position_z", &mc_pos_z_n, "mc_pos_z_n/D");
    tree_output->Branch("mc_ev_index_ep", &mc_ev_index_ep, "mc_ev_index_ep/l");
    tree_output->Branch("mc_ev_index_n", &mc_ev_index_n, "mc_ev_index_n/l");
    tree_output->Branch("ev_index", &ev_index, "ev_index/l");
    tree_output->Branch("ev_fit_validity", &ev_validity, "ev_validity/O");
    tree_output->Branch("ev_fit_energy", &ev_energy, "ev_energy/D");
    tree_output->Branch("ev_fit_position_r", &ev_pos_r, "ev_pos_r/D");
    tree_output->Branch("ev_fit_position_x", &ev_pos_x, "ev_pos_x/D");
    tree_output->Branch("ev_fit_position_y", &ev_pos_y, "ev_pos_y/D");
    tree_output->Branch("ev_fit_position_z", &ev_pos_z, "ev_pos_z/D");
    tree_output->Branch("ev_nhit", &ev_nhit, "ev_nhit/I");
    tree_output->Branch("ev_time_days", &ev_time_days, "ev_time_days/i");
    tree_output->Branch("ev_time_seconds", &ev_time_seconds, "ev_time_seconds/i");
    tree_output->Branch("ev_time_nanoseconds", &ev_time_nanoseconds, "ev_time_nanoseconds/D");
    tree_output->Branch("reactor_info_latitude", &latitude_i, "latitude_i/D");
    tree_output->Branch("reactor_info_longitude", &longitude_i, "longitude_i/D");
    tree_output->Branch("reactor_info_altitude", &altitude_i, "altitude_i/D");
    tree_output->Branch("reactor_info_distance", &distance_i, "distance_i/D");

    // values to modify
    tree_output->Branch("ev_fit_energy_p1", &ev_energy_p1, "ev_energy_p1/D");
    tree_output->Branch("ev_fit_energy_p2", &ev_energy_p2, "ev_energy_p2/D");
    tree_output->Branch("ev_index_p1", &ev_index_p1, "ev_index_p1/l");
    tree_output->Branch("ev_index_p2", &ev_index_p2, "ev_index_p2/l");

    // load ds
    RAT::DU::DSReader ds_reader(input_filename.c_str());
    ULong64_t percent_interval = ds_reader.GetEntryCount()/10; //print every 10%
    if (percent_interval<=0) percent_interval = 1;
    ULong64_t progress_countdown = percent_interval;
    std::cout << "Number of events: " << ds_reader.GetEntryCount() << std::endl;

    for (ULong64_t i_entry = 0; i_entry < ds_reader.GetEntryCount(); i_entry++) {

        // entry index
        const RAT::DS::Entry& ds_entry = ds_reader.GetEntry(i_entry);
        entry = i_entry;

        // print progress
        progress_countdown--;
        if (progress_countdown==0){
            progress_countdown = percent_interval;
            printf("%.0f%% done\n",(Double_t)(i_entry+1)/ds_reader.GetEntryCount()*100);
        }

        //*** MC entries
        // pdg code -12 = anti electron neutrino, see http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
        // pdg code -11 = positron
        // pdg code 2112 = neutron

        const RAT::DS::MC &rMC = ds_entry.GetMC();
        ULong64_t n_mcparent = 0; // only ever 1 parent //for(ULong64_t i_mcparent = 0; i_mcparent < rMC.GetMCParentCount(); i_mcparent++) {

        // set all variables to some unphysical value
        mc_energy_nu = -9000;
        mc_pos_x_nu = -9000;
        mc_pos_y_nu = -9000;
        mc_pos_z_nu = -9000;
        mc_pos_r_nu = -9000;
        mc_quench_i = -9000;
        mc_energy_n = -9000;
        mc_pos_x_n = -9000;
        mc_pos_y_n = -9000;
        mc_pos_z_n = -9000;
        mc_pos_r_n = -9000;
        mc_energy_ep = -9000;
        mc_pos_x_ep = -9000;
        mc_pos_y_ep = -9000;
        mc_pos_z_ep = -9000;
        mc_pos_r_ep = -9000;
        latitude_i = -9000;
        longitude_i = -9000;
        altitude_i = -9000;
        distance_i = -9000;
        mc_ev_index_ep = -9000;
        mc_ev_index_n = -9000;
        ev_energy_p1 = -9000;
        ev_energy_p2 = -9000;
        ev_index_p1 = -9000;
        ev_index_p2 = -9000;

        if ((rMC.GetMCParticleCount()==2)&&(rMC.GetMCParent(n_mcparent).GetPDGCode()==-12)){ // check the parent is an anti-neutrino and there are two child particles

            const RAT::DS::UniversalTime &time_day_sec_ns = rMC.GetUniversalTime();
            mc_time_days = time_day_sec_ns.GetDays();
            mc_time_seconds = time_day_sec_ns.GetSeconds();
            mc_time_nanoseconds = time_day_sec_ns.GetNanoSeconds();

            mc_quench_i = rMC.GetScintQuenchedEnergyDeposit();

            // parent (anti-neutrino) particle properties
            const RAT::DS::MCParticle &mc_parent = rMC.GetMCParent(0);
            mc_energy_nu = mc_parent.GetKineticEnergy();
            mc_pos_r_nu = mc_parent.GetPosition().Mag();
            mc_pos_x_nu = mc_parent.GetPosition().X();
            mc_pos_y_nu = mc_parent.GetPosition().Y();
            mc_pos_z_nu = mc_parent.GetPosition().Z();

            // child particle properties
            for(ULong64_t i_mcparticle = 0; i_mcparticle < rMC.GetMCParticleCount(); i_mcparticle++) {
                const RAT::DS::MCParticle &mc_particle = rMC.GetMCParticle(i_mcparticle);
                if (mc_particle.GetPDGCode()==-11){  // positron properties
                    mc_energy_ep = mc_particle.GetKineticEnergy();
                    mc_pos_r_ep = mc_particle.GetPosition().Mag();
                    mc_pos_x_ep = mc_particle.GetPosition().X();
                    mc_pos_y_ep = mc_particle.GetPosition().Y();
                    mc_pos_z_ep = mc_particle.GetPosition().Z();
                    mc_ev_index_ep = i_mcparticle;
                }
                if (mc_particle.GetPDGCode()==2112){  // neutron properties
                    mc_energy_n = mc_particle.GetKineticEnergy();
                    mc_pos_r_n = mc_particle.GetPosition().Mag();
                    mc_pos_x_n = mc_particle.GetPosition().X();
                    mc_pos_y_n = mc_particle.GetPosition().Y();
                    mc_pos_z_n = mc_particle.GetPosition().Z();
                    mc_ev_index_n = i_mcparticle;
                }
            }
            //*** DB entries
            // lat and long from REACTORS.ratdb
            // (split reactor name into name and core components)
            reactorcoreStr = mc_parent.GetMetaInfo();
            Int_t pos = reactorcoreStr.find_last_of(" ");
            reactorName = reactorcoreStr.substr(0, pos);
            Int_t coreNumber = atoi(reactorcoreStr.substr(pos+1, reactorcoreStr.size()).c_str());
            fLink = db->GetLink("REACTOR", reactorName);
            latitude  = fLink->GetDArray("latitude");
            longitude = fLink->GetDArray("longitude");
            altitude  = fLink->GetDArray("altitude");
            latitude_i = latitude[coreNumber];
            longitude_i = longitude[coreNumber];
            altitude_i = altitude[coreNumber];
            distance_i = GetReactorDistanceLLA(latitude_i, longitude_i, altitude_i);
        }

        //*** EV tree entries
        for(ULong64_t i_ev = 0; i_ev < ds_entry.GetEVCount(); i_ev++) {

            const RAT::DS::EV &rEV = ds_entry.GetEV(i_ev);
            ev_nhit = rEV.GetNhitsCleaned();
            const RAT::DS::UniversalTime &time_day_sec_ns = rEV.GetUniversalTime();
            ev_time_days = time_day_sec_ns.GetDays();
            ev_time_seconds = time_day_sec_ns.GetSeconds();
            ev_time_nanoseconds = time_day_sec_ns.GetNanoSeconds();
            ev_validity = false;
            ev_index = i_ev;

            // set all to some unphysical value, maintain correspondence of number of events in tree to events in root file
            ev_energy = -9000;
            ev_pos_x = -9000;
            ev_pos_y = -9000;
            ev_pos_z = -9000;
            ev_pos_r = -9000;

            const vector<std::string> fit_names = rEV.GetFitNames();
            for (ULong64_t i = 0; i < fit_names.size(); i++){ //ensure fit exists by going through all fits and ..
                if (fit_names.at(i) == fitter_name.c_str()) { // ..selecting the desired name
                    try {
                        const RAT::DS::FitResult &fitResult = rEV.GetFitResult(fitter_name.c_str());
                        const RAT::DS::FitVertex &fitVertex = fitResult.GetVertex(0);
                        ev_validity = fitResult.GetValid() && fitVertex.ContainsPosition() && fitVertex.ValidPosition() && fitVertex.ContainsEnergy() && fitVertex.ValidEnergy();
                        ev_energy = fitVertex.GetEnergy();
                        ev_pos_x = fitVertex.GetPosition().X();
                        ev_pos_y = fitVertex.GetPosition().Y();
                        ev_pos_z = fitVertex.GetPosition().Z();
                        ev_pos_r = fitVertex.GetPosition().Mag();
                    }
                    catch (RAT::DS::FitResult::NoVertexError &e){;}
                    catch (RAT::DS::FitResult::NoFitResultError &e){;}
                    catch (RAT::DS::FitVertex::NoValueError &e){;}
                    break;
                }
            }
            tree_output->Fill();
            n_total_events++;
        }
    }
    std::cout << " Processed events: " << ds_reader.GetEntryCount() << std::endl;
    std::cout << " Total number of events (all ev's): " << n_total_events << std::endl;

    // write output ntuple
    TFile *output_file = new TFile(output_filename.c_str(),"RECREATE");
    tree_output->Write();
    output_file->Close();
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cout<<"Error: 3 arguments expected. Got: "<<argc-1<<std::endl;
        return 1; // return>0 indicates error code
    }
    else{
        std::string input_filename = argv[1];
        std::string output_filename = argv[2];
        std::string fitter_name = argv[3];
        ntload(input_filename, output_filename, fitter_name);
        return 0; // return=0 indicates no error
    }
}
