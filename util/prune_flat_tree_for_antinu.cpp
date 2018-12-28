#include <TFile.h>
#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/FitResult.hh>
#include <RAT/DB.hh>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TMath.h>
#include <TCanvas.h>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <TNtuple.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <RAT/DU/Utility.hh>
#include <RAT/TrackNav.hh>
#include <RAT/TrackCursor.hh>
#include <RAT/TrackNode.hh>
#include "TView.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include <RAT/DB.hh>
#include <TObject.h>

const TVector3 SNO_LLA_coord_ = TVector3(-81.2014, 46.4753, -1766.0);
const TVector3 SNO_ECEF_coord_ = TVector3(672.87, -4347.18, 4600.51);

double CalculateDistance(TVector3 point1, TVector3 point2) {
    return (point2 - point1).Mag();
}

TVector3 LLAtoECEF(double longitude, double latitude, double altitude) {
    // reference http://www.mathworks.co.uk/help/aeroblks/llatoecefposition.html
    static double toRad = TMath::Pi()/180.;
    static double Earthradius = 6378137.0; //Radius of the Earth (in meters)
    static double f = 1./298.257223563; //Flattening factor WGS84 Model
    static double L, rs, x, y, z;
    L = atan( pow((1. - f),2)*TMath::Tan(latitude*toRad))*180./TMath::Pi();
    rs = TMath::Sqrt( pow(Earthradius,2)/(1. + (1./pow((1. - f),2) - 1.)*pow(TMath::Sin(L*toRad),2)));
    x = (rs*TMath::Cos(L*toRad)*TMath::Cos(longitude*toRad) + altitude*TMath::Cos(latitude*toRad)*TMath::Cos(longitude*toRad))/1000; // in km
    y = (rs*TMath::Cos(L*toRad)*TMath::Sin(longitude*toRad) + altitude*TMath::Cos(latitude*toRad)*TMath::Sin(longitude*toRad))/1000; // in km
    z = (rs*TMath::Sin(L*toRad) + altitude*TMath::Sin(latitude*toRad))/1000; // in km

    TVector3 ECEF(x,y,z);

    return ECEF;
}

double GetReactorDistanceLLA(double longitude, double latitude, double altitude) {
    return CalculateDistance(SNO_ECEF_coord_,LLAtoECEF(longitude, latitude,altitude));
}

void ntload(std::string input_filename, std::string output_filename) {

    //TNtuple* nt = new TNtuple("nt","nt","entry:MCQuench:ParPosr:ParKE:Part1Posr:Part1KE:ReactorLatitude:ReactorLongitude:ReactorDistance");
    TTree *tt = new TTree("tt", "AntinuTree");
    //Float_t v[9];

    size_t ientry=0;
    std::vector<float> MCenergy;
    std::vector<float> MCposition;
    std::vector<float> EVenergy;
    std::vector<float> EVposition;
    std::vector<float> EVnhit;
    std::vector<float> EVtime;
    std::vector<float> reactorInfo;
    tt->Branch("entry",ientry);
    tt->Branch("MCenergy",&MCenergy);
    tt->Branch("MCposition",&MCposition);
    tt->Branch("EVenergy",&EVenergy);
    tt->Branch("EVposition",&EVposition);
    tt->Branch("EVnhit",&EVnhit);
    tt->Branch("EVtime",&EVtime);
    tt->Branch("reactorInfo",&reactorInfo);

    std::string reactorName="";
    std::string reactorcoreStr="";
    RAT::DB *db = RAT::DB::Get();
    RAT::DBLinkPtr fLink;
    std::vector<double> fLatitude;
    std::vector<double> fLongitude;
    std::vector<double> fAltitude;

    // load ds
    RAT::DU::DSReader ds_reader(input_filename.c_str());

    for (size_t i_entry=0; i_entry < ds_reader.GetEntryCount(); ++i_entry) {

        // entry index
        const RAT::DS::Entry& ds_entry = ds_reader.GetEntry(i_entry);

        //*** MC entries
        // pdg code -12 = anti electron neutrino, see http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
        // pdg code -11 = positron
        // pdg code 2112 = neutron

        const RAT::DS::MC &rMC = ds_entry.GetMC();
        for(size_t i_mcparent = 0; i_mcparent < rMC.GetMCParentCount(); i_mcparent++ ) {

            Double_t mc_energy, mc_x, mc_y, mc_z, mc_mag, mc_quench;
            Double_t mc_n_energy, mc_n_x, mc_n_y, mc_n_z, mc_n_mag;
            Double_t mc_ep_energy, mc_ep_x, mc_ep_y, mc_ep_z, mc_ep_mag;
            double latitude, longitude, altitude, distance;

            // set all variables to some unphysical value
            mc_energy = -9000;
            mc_x = -9000;
            mc_y = -9000;
            mc_z = -9000;
            mc_mag = -9000;
            mc_quench = -9000;
            mc_n_energy = -9000;
            mc_n_x = -9000;
            mc_n_y = -9000;
            mc_n_z = -9000;
            mc_n_mag = -9000;
            mc_ep_energy = -9000;
            mc_ep_x = -9000;
            mc_ep_y = -9000;
            mc_ep_z = -9000;
            mc_ep_mag = -9000;
            latitude = -9000;
            longitude = -9000;
            altitude = -9000;
            distance = -9000;

            if ((rMC.GetMCParticleCount()==2)&&(rMC.GetMCParent(i_mcparent).GetPDGCode()==-12)){ // check the parent is an anti-neutrino and there are two child particles

                mc_quench = rMC.GetScintQuenchedEnergyDeposit();

                // parent (antineutrino) particle properties
                const RAT::DS::MCParticle &mc_parent = rMC.GetMCParent(0);
                mc_energy = mc_parent.GetKineticEnergy();
                mc_mag = mc_parent.GetPosition().Mag();
                mc_x = mc_parent.GetPosition().X();
                mc_y = mc_parent.GetPosition().Y();
                mc_z = mc_parent.GetPosition().Z();

                // child particle properties
                for( size_t i_mcparticle = 0; i_mcparticle < rMC.GetMCParticleCount(); i_mcparticle++ ) {
                    const RAT::DS::MCParticle &mc_particle = rMC.GetMCParticle(i_mcparticle);
                    if (mc_particle.GetPDGCode()==2112){  // positron properties
                        mc_ep_energy = mc_particle.GetKineticEnergy();
                        mc_ep_mag = mc_particle.GetPosition().Mag();
                        mc_ep_x = mc_particle.GetPosition().X();
                        mc_ep_y = mc_particle.GetPosition().Y();
                        mc_ep_z = mc_particle.GetPosition().Z();
                    }
                    if (mc_particle.GetPDGCode()==-11){  // neutron properties
                        mc_n_energy = mc_particle.GetKineticEnergy();
                        mc_n_mag = mc_particle.GetPosition().Mag();
                        mc_n_x = mc_particle.GetPosition().X();
                        mc_n_y = mc_particle.GetPosition().Y();
                        mc_n_z = mc_particle.GetPosition().Z();
                    }
                }
                //*** DB entries
                // lat and long from REACTORS.ratdb
                // (split reactor name into name and core components)
                reactorcoreStr = mc_parent.GetMetaInfo();
                int pos = reactorcoreStr.find_last_of(" ");
                reactorName = reactorcoreStr.substr(0, pos);
                int coreNumber = atoi(reactorcoreStr.substr(pos+1, reactorcoreStr.size()).c_str());
                fLink = db->GetLink("REACTOR", reactorName);
                fLatitude  = fLink->GetDArray("latitude");
                fLongitude = fLink->GetDArray("longitude");
                fAltitude  = fLink->GetDArray("altitude");
                latitude = fLatitude[coreNumber];
                longitude = fLongitude[coreNumber];
                altitude = fAltitude[coreNumber];
                distance = GetReactorDistanceLLA(latitude, longitude, altitude);
            }

            //add to tree
            ientry = i_entry;
            //particle order convention: antineutrino, positron, neutron.
            //energy order convention: energy, quench.
            MCenergy.push_back(mc_energy);
            MCenergy.push_back(mc_ep_energy);
            MCenergy.push_back(mc_n_energy);
            MCenergy.push_back(mc_quench);
            //position particle order convention: mag, x, y, z.
            MCposition.push_back(mc_mag);
            MCposition.push_back(mc_x);
            MCposition.push_back(mc_y);
            MCposition.push_back(mc_z);
            MCposition.push_back(mc_ep_mag);
            MCposition.push_back(mc_ep_x);
            MCposition.push_back(mc_ep_y);
            MCposition.push_back(mc_ep_z);
            MCposition.push_back(mc_n_mag);
            MCposition.push_back(mc_n_x);
            MCposition.push_back(mc_n_y);
            MCposition.push_back(mc_n_z);
            //info order convention: latitude, longitude, altitude, distance.
            reactorInfo.push_back(latitude);
            reactorInfo.push_back(longitude);
            reactorInfo.push_back(altitude);
            reactorInfo.push_back(distance);
            //v[0] = (Float_t)(i_entry);
            //v[1] = (Float_t)(mc_quench);
            //v[2] = (Float_t)(mc_mag);
            //v[3] = (Float_t)(mc_energy);
            //v[4] = (Float_t)(mc_ep_mag);
            //v[5] = (Float_t)(mc_ep_energy);
            //v[6] = latitude;
            //v[7] = longitude;
            //v[8] = distance;
        }

        //*** EV tree entries
        for(size_t i_ev = 0; i_ev < ds_entry.GetEVCount(); i_ev++) {

            const RAT::DS::EV &rEV = ds_entry.GetEV(i_ev);

            const int nHit = rEV.GetNhits();
            const RAT::DS::UniversalTime &time_day_sec_ns = rEV.GetUniversalTime();
            const unsigned int time_days = time_day_sec_ns.GetDays();
            const unsigned int time_seconds = time_day_sec_ns.GetSeconds();
            const Double_t time_nanoSeconds = time_day_sec_ns.GetNanoSeconds();
            Double_t vertex_energy, vertex_x, vertex_y, vertex_z, vertex_mag;
            bool vertex_validity = false;

            // set all to some unphysical value, maintain correspondence of number of events in tree to events in root file
            vertex_energy = -9000;
            vertex_x = -9000;
            vertex_y = -9000;
            vertex_z = -9000;
            vertex_mag = -9000;

            const unsigned long nFits = rEV.GetFitNames().size(); //use this to avoid trivial NoFitResultError exceptions
            if (nFits>0){
                try {
                    const RAT::DS::FitResult &fitResult = rEV.GetFitResult("scintFitter");
                    const RAT::DS::FitVertex &fitVertex = fitResult.GetVertex(0);
                    vertex_validity = fitResult.GetValid();
                    vertex_energy = fitVertex.GetEnergy();
                    vertex_x = fitVertex.GetPosition().X();
                    vertex_y = fitVertex.GetPosition().Y();
                    vertex_z = fitVertex.GetPosition().Z();
                    vertex_mag = fitVertex.GetPosition().Mag();
                }
                catch (RAT::DS::FitResult::NoFitResultError &e){;}
                catch (RAT::DS::FitVertex::NoValueError &e){;}
            }
            else{
                vertex_validity = false;
            }

            //add all events to tree
            EVenergy.push_back((int)vertex_validity);
            EVenergy.push_back(vertex_energy);
            //position particle order convention: mag, x, y, z.
            EVposition.push_back(vertex_x);
            EVposition.push_back(vertex_y);
            EVposition.push_back(vertex_z);
            EVposition.push_back(vertex_mag);
            EVnhit.push_back(nHit);
            EVtime.push_back(time_days);
            EVtime.push_back(time_seconds);
            EVtime.push_back(time_nanoSeconds);
        }

        if (i_entry % 1000 == 0) std::cout <<  "    Processed events: " << i_entry << " (" << (double)i_entry/ds_reader.GetEntryCount()*100. << "%) " << std::endl;

        // fill
        tt->Fill();
        //nt->Fill(v);
    }
    std::cout <<  " Processed events: " << ds_reader.GetEntryCount() << std::endl;

    // write output ntuple
    TFile *output_file = new TFile(output_filename.c_str(),"RECREATE");
    //nt->Write();
    tt->Write();
    output_file->Close();
}

int main(int argc, char* argv[]) {
  std::string input_filename = argv[1];
  std::string output_filename = argv[2];
  ntload(input_filename, output_filename);
}
