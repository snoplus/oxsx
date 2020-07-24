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
#include <algorithm>
#include <TRandom3.h>

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
  //TAKING ONLYEVENTS WITH FITRESTULTS AND VALID FITTED POSITION
  TNtuple* nt = new TNtuple("nt","nt","MCentry:Evindex:ParKE:Part1KE:Part2KE:ReactorLatitude:ReactorLongitude:ReactorDistance:Days:Seconds:Nanoseconds:nhits:Energy:posx:posy:posz:Fitted");
  int Entries = 17;
  Float_t v[Entries];
  int Fit = 0;
  int NoFit = 0;
  RAT::DU::DSReader ds(fname);
  
  std::string rName="";
  int coreNumber;
  std::string reactorName="";
  std::vector<std::string> reactorNames;
  std::vector<int> reactorNorms;
  std::string reactorcoreStr="";
  int pos;

  RAT::DB *db = RAT::DB::Get();
  RAT::DBLinkPtr fLink;
  std::vector<double> fLatitude;
  std::vector<double> fLongitude;
  std::vector<double> fAltitude;
  double fDistance = 0.;

  TRandom3 *r1 = new TRandom3();
  r1->SetSeed(0);
  
  int evcount = 0;
  //int sec = 1;
  int Nevents = abs(ds.GetEntryCount());
  for (int ie=0; ie<Nevents; ie++) {  
    
    // number of entries
    const RAT::DS::Entry& rds = ds.GetEntry(ie);
    const RAT::DS::MC& rMC = rds.GetMC();
    const RAT::DS::MCParticle& rmcparent = rMC.GetMCParent(0);
    //int parentpdg = rmcparent.GetPDGCode();
    //std::cout<<parentpdg<<std::endl;
    
    const RAT::DS::MCParticle& rmcparticle1 = rMC.GetMCParticle(0);
    const RAT::DS::MCParticle& rmcparticle2 = rMC.GetMCParticle(1);
      
    reactorcoreStr = rmcparent.GetMetaInfo();
    pos = reactorcoreStr.find_last_of(" ");
    reactorName = reactorcoreStr.substr(0, pos);
      
    coreNumber = atoi(reactorcoreStr.substr(pos+1, reactorcoreStr.size()).c_str());
    fLink = db->GetLink("REACTOR", reactorName);
    fLatitude  = fLink->GetDArray("latitude");
    fLongitude = fLink->GetDArray("longitude");
    fAltitude  = fLink->GetDArray("altitude");    
    // distance
    fDistance = GetReactorDistanceLLA( fLongitude[coreNumber],fLatitude[coreNumber],fAltitude[coreNumber] );
      
    int EVcount = rds.GetEVCount();
    for( int iEV = 0; iEV < EVcount; iEV++ ){
      const RAT::DS::EV& rEV = rds.GetEV(iEV);
      evcount += 1;
      v[0] = (Float_t)(ie);
      v[1] = (Float_t)(iEV);
      v[2] = (Float_t)(rmcparent.GetKineticEnergy());
      v[3] = (Float_t)(rmcparticle1.GetKineticEnergy());	
      v[4] = (Float_t)(rmcparticle2.GetKineticEnergy());	
      v[5] = (Float_t)(fLatitude[coreNumber]);
      v[6] = (Float_t)(fLongitude[coreNumber]);
      v[7] = (Float_t)(fDistance);	
	
      bool BadFit = true;
      if( rEV.FitResultExists("scintFitter")){
	if( rEV.GetFitResult("scintFitter").GetVertex(0).ContainsPosition() ){
	  if( rEV.GetFitResult("scintFitter").GetVertex(0).ValidPosition() ) {
	    if( rEV.GetFitResult("scintFitter").GetVertex(0).ContainsEnergy() ) {
	      if( rEV.GetFitResult("scintFitter").GetVertex(0).ValidEnergy() ) {
		BadFit = false;
		//std::cout<<"MCentry: "<<ie<<" Evindex "<<iEV<<" Fitted +1"<<std::endl;
		RAT::DS::FitResult fResult = rEV.GetFitResult("scintFitter");
		RAT::DS::FitVertex fVertex = fResult.GetVertex(0);
		v[8] = (Float_t)(rEV.GetUniversalTime().GetDays());
		v[9] = (Float_t)(rEV.GetUniversalTime().GetSeconds());
		v[10] = (Float_t)(rEV.GetUniversalTime().GetNanoSeconds());
		v[11] = (Float_t)(rEV.GetNhits());
		v[12] = (Float_t)(fVertex.GetEnergy());
		v[13] = (Float_t)(fVertex.GetPosition().X());
		v[14] = (Float_t)(fVertex.GetPosition().Y());
		v[15] = (Float_t)(fVertex.GetPosition().Z());
		v[16] = (Float_t)(1);
		Fit +=1;
	      }
	    }
	  }
	}
      }
	
      if (BadFit){
	//std::cout<<"MCentry: "<<ie<<" Evindex "<<iEV<<" No Fit: 0"<<std::endl;
	for (int i = 8; i < Entries; i++)
	  v[i] = (Float_t)(-1);
	NoFit += 1;
      }
      nt->Fill(v);
    }
    
    if (ie % 100 == 0) std::cout <<  "    Processed events: " << ie << " (" << (double)ie/Nevents*100. << "%) " << std::endl;
  }
  
  std::cout <<  " Processed events: " << Nevents << std::endl;
  std::cout <<  " Fits: " <<Fit<<" NoFits: "<<NoFit<< std::endl;
  
  // write output ntuple
  std::stringstream outname;
  outname<<nout;
  TFile *ntout = new TFile(outname.str().c_str(),"RECREATE");
  std::cout<<outname.str()<<std::endl;
  nt->Write();
  ntout->Close();
  return nt;
}

int main(int argc, char* argv[])
{
  const char* rootin= argv[1];
  const char* ntout = argv[2];
  if (argc != 3) {
    std::cout<<"2 arguments expected: \n 1: location of input root file \n 2: location/filename\
 of output pruned ntuple"<<std::endl;
  }
  else{
  ntload(rootin, ntout);
  }
}
