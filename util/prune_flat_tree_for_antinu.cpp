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

const TVector3 SNO_LLA_coord_ = TVector3(-81.2014, 46.4753,-1766.0);
const TVector3 SNO_ECEF_coord_ = TVector3(672.87,-4347.18,4600.51);

double CalculateDistance( TVector3 point1, TVector3 point2) {
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

TNtuple* ntload(const char* fname, const char* nout) {

  TNtuple* nt = new TNtuple("nt","nt","entry:MCQuench:ParPosr:ParKE:Part1Posr:Part1KE:ReactorLatitude:ReactorLongitude:ReactorDistance");
  Float_t v[9];

  RAT::DU::DSReader ds(fname);
  
  std::string rName="";
  int coreNumber;
  std::string reactorName="";
  std::string reactorcoreStr="";
  int pos;

  RAT::DB *db = RAT::DB::Get();
  RAT::DBLinkPtr fLink;
  std::vector<double> fLatitude;
  std::vector<double> fLongitude;
  std::vector<double> fAltitude;
  double fDistance = 0.;
  
  size_t ie = 0;
  size_t Nevents = ds.GetEntryCount();

  for (ie=0; ie<Nevents; ++ie) {
      
    // number of entries
    const RAT::DS::Entry& entry = ds.GetEntry(ie);
    v[0] = (Float_t)(ie);
    
    // quenched energy
    const RAT::DS::MC& rMC = entry.GetMC();
    v[1] = (Float_t)(rMC.GetScintQuenchedEnergyDeposit());
    
    // particle position and energy
    const RAT::DS::MCParticle& rmcparent = rMC.GetMCParent(0);
    const RAT::DS::MCParticle& rmcparticle1 = rMC.GetMCParticle(0);
    v[2] = (Float_t)((rmcparent.GetPosition()).Mag());
    v[3] = (Float_t)(rmcparent.GetKineticEnergy());
    v[4] = (Float_t)((rmcparticle1.GetPosition()).Mag());
    v[5] = (Float_t)(rmcparticle1.GetKineticEnergy());

    // lat and long from REACTORS.ratdb
    // (split reactor name into name and core components)
    reactorcoreStr = rmcparent.GetMetaInfo();
    pos = reactorcoreStr.find_last_of(" ");
    reactorName = reactorcoreStr.substr(0, pos);
    coreNumber = atoi(reactorcoreStr.substr(pos+1, reactorcoreStr.size()).c_str());
    fLink = db->GetLink("REACTOR", reactorName);
    fLatitude  = fLink->GetDArray("latitude");
    fLongitude = fLink->GetDArray("longitude");
    fAltitude  = fLink->GetDArray("altitude");
    v[6] = fLatitude[coreNumber];
    v[7] = fLongitude[coreNumber];
    
    // distance
    fDistance = GetReactorDistanceLLA( fLongitude[coreNumber],fLatitude[coreNumber],fAltitude[coreNumber] );
    v[8] = fDistance;
    //std::cout <<  "Name=" << reactorName << " core=" << coreNumber+1 << " longitude=" << v[6] << " latitude=" << v[7] << " distance=" << v[8] << std::endl;

    if (ie % 1000 == 0) std::cout <<  "    Processed events: " << ie << " (" << (double)ie/Nevents*100. << "%) " << std::endl;

    // fill ntuple
    nt->Fill(v);
  }
  
  // 
  std::cout <<  " Processed events: " << Nevents << std::endl;
  
  // write output ntuple
  std::stringstream outname;
  outname<<nout;
  TFile *ntout = new TFile(outname.str().c_str(),"RECREATE");
  nt->Write();
  ntout->Close();
  return nt;
}

int main(int argc, char* argv[])
{
  const char* rootin = argv[1];
  const char* ntout = argv[2];
  ntload(rootin, ntout);
}
