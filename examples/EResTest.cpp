// A simple fit in energy for signal and a background
#include <string>
#include <vector>
#include <math.h>
#include <Rand.h>
#include <fstream>
#include <iostream>

#include <TH1D.h>
#include <TH2D.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <ROOTNtuple.h>
#include <TFile.h>
#include <TMath.h>
#include <string>
#include <TNtuple.h>
#include <iostream>
#include <TObject.h>
#include <TRandom3.h>
#include <TVector3.h>

#include <BinnedED.h>
#include <BinnedEDGenerator.h>
#include <SystematicManager.h>
#include <BinnedNLLH.h>
#include <FitResult.h>
#include <Minuit.h>
#include <DistTools.h>
#include <Minuit.h>
#include <Convolution.h>
#include <Scale.h>
#include <BoolCut.h>
#include <BoxCut.h>
#include <ParameterDict.h>
#include <ContainerTools.hpp>
#include <ROOTMultiPlot.h>
#include <TRandom3.h>
#include <Gaussian.h>
#include <GaussianERes.h>

// flag which is used for repeating fits that weren't good
bool badfit = true;

int unoscPdfint;
int oscPdfint;

//Surv Probablity function for antinu of KE nuE and distance traveled of baseline.
double probability(double nuE, double baseline, double delmsqr21, double sinsqrtheta12, double sinsqrtheta13) {
  double fSSqr2Theta12 = pow(sin(2.0 * TMath::ASin(sqrt(sinsqrtheta12))), 2.0);
  double fS4 = pow(sinsqrtheta13, 2.0);
  double fC4 = pow(1.0-sinsqrtheta13, 2.0);
  double scale = 1.267e3; // for nuE in [MeV] and baseline in [km]
  double sSqrDmBE = pow(sin(scale * delmsqr21 * baseline / nuE), 2.0);
  double fOscProb = (fC4 * (1.0 - fSSqr2Theta12 * sSqrDmBE) + fS4);
  return fOscProb;
}

//Coincidence tagging - function to apply cuts to isolate positron and neutron events
// At present, looking at the MC event number and EV index number to ensure no pairs are missed
void OscPromptE (const std::string infile, const std::string outfile,  double delmsqr21, double sinsqrtheta12, double sinsqrtheta13, double Distance) {
  
  unoscPdfint = 0;
  oscPdfint = 0;
  
  TFile *fin = TFile::Open(infile.c_str());
  TNtuple *T = (TNtuple *) fin->Get("nt");

  TFile *fout = new TFile(outfile.c_str(), "RECREATE");
  TNtuple* ntout = new TNtuple("nt","nt","E1");

  float parke, Energy;
  
  T->SetBranchAddress("ParKE", &parke);
  T->SetBranchAddress("E1", &Energy);
  //T->SetBranchAddress("ReactorDistance", &ReactorDistance); // Average distance taken from info file, not event by event (though it can be)

  double survprob;
  TRandom3 *r1 = new TRandom3();
  r1->SetSeed(0);    

  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    
    survprob = probability(parke, Distance, delmsqr21, sinsqrtheta12, sinsqrtheta13);
    unoscPdfint +=1;
    if (survprob > r1->Rndm()){
      oscPdfint +=1;
      ntout->Fill(Energy);
    }
  }

  ntout->Write();
  fin->Close();
  fout->Close();
  delete fin;
  delete fout;
}
  

void readInfoFile(const std::string &runInfoFileName, std::vector<std::string> &reactorNames, std::vector<double> &distances, std::vector<std::string> &reactorTypes, std::vector<int> &nCores, std::vector<double> &powers ) {
  // Read couchDB run-info text file
  std::ifstream in;
  in.open(runInfoFileName.c_str());
  std::cout << "opening file: " << runInfoFileName.c_str() << std::endl;

  std::fill(reactorNames.begin(), reactorNames.end(), "");
  std::fill(distances.begin(), distances.end(), 0.);
  std::fill(reactorTypes.begin(), reactorTypes.end(), "");
  std::fill(nCores.begin(), nCores.end(), 0);
  std::fill(powers.begin(), powers.end(), 0.);

  std::string reactorName,distance,reactorType,nCore,power;
  int lineNo = 0;

  // read until end of file.
  while(in.good()){
    std::getline(in,reactorName,',');
    std::getline(in,distance,',');
    std::getline(in,reactorType,',');
    std::getline(in,nCore,',');
    std::getline(in,power,'\n');

    if (lineNo>0){ //skip csv header
      if (strcmp(reactorName.c_str(),"")!=0) {
	reactorNames.push_back(reactorName);
	distances.push_back(atof(distance.c_str()));
	reactorTypes.push_back(reactorType.c_str());
	nCores.push_back(atoi(nCore.c_str()));
	powers.push_back(atof(power.c_str()));

	std::cout << "v: reactorName: " << reactorNames[lineNo-1] << ", distance: " << distances[lineNo-1] << ", reactorType: " << reactorTypes[lineNo-1] << ", nCore: " << nCores[lineNo-1] << ", power: " << powers[lineNo-1] << std::endl; //debug check ('-1' for header)
      }
    }
    lineNo++;
  }
  in.close();
}

double LHFit(const std::string UnOscfile, const std::string tempFile, AxisCollection axes, BinnedED dataSetPdf, int numPdfs, std::vector<std::string> reactorNames, std::vector<double> reactorDistances, double d21fix, double s12fix,std::vector<double> normconstraints,int Normmin, int Normmax, std::vector<double> normconstraintsuncert){
  
  char name[100];
  TRandom3 *r1 = new TRandom3();
  r1->SetSeed(0);

  ////////////////////
  // 1. Set Up PDFs //
  ////////////////////
  
  
  ParameterDict minima;
  ParameterDict maxima;
  ParameterDict initialval;
  ParameterDict initialerr;
 
  BinnedNLLH lhFunction;
  lhFunction.SetBufferAsOverflow(true);
  int Buff = 1;
  lhFunction.SetBuffer(0,Buff,Buff);  
  lhFunction.SetDataDist(dataSetPdf); // initialise withe the data set
  double rand = r1->Rndm();

  // allowing for oscillation of constraints, calculating the generator's predicted oscillated total integral
  double oscconstraintint = 0.;
  for (int i = 0; i < numPdfs; i++){
    OscPromptE(UnOscfile,tempFile ,d21fix,s12fix,0.,reactorDistances[i]); //or s13 = 0.0215
    double osc_loss = oscPdfint/(double)unoscPdfint;      // calculate loss in pdf integral
    std::cout<<" osc loss factor: "<<osc_loss<<std::endl;
    normconstraints[i] = normconstraints[i]*osc_loss;
    oscconstraintint += normconstraints[i];
  }
  
  //renomalising oscillated constraints, to KL Oscillated spectrum integral
  double DataInt = dataSetPdf.Integral();
  std::cout<<" Data Integral: "<<DataInt<<std::endl;
  double oscrenormInt = 0.;
  for (int i = 0; i < numPdfs; i++){
    normconstraints[i] = ((normconstraints[i]/oscconstraintint)*DataInt); 
    oscrenormInt += normconstraints[i];
    std::cout<<reactorNames[i]<<"   mean: "<<normconstraints[i]<<"  sigma: "<<(normconstraintsuncert[i]*normconstraints[i])<<std::endl;
  }
  
  std::cout<<"Integral of measured (oscillated) KL Data: "<<DataInt<<std::endl;
  std::cout<<"Integral of oscillated renormalised constraints: "<<oscrenormInt<<std::endl;
  


  //Adding convolution
  /*ObsSet  obsSet(0);
  Convolution* conv = new Convolution("conv");
  Gaussian* gaus = new Gaussian(0,1,"gaus"); 
  gaus->RenameParameter("means_0","mean_");
  gaus->RenameParameter("stddevs_0","stddev_");
  conv->SetFunction(gaus);
  conv->SetAxes(axes);
  conv->SetTransformationObs(obsSet);
  conv->SetDistributionObs(obsSet);
  conv->Construct();
  
  minima["mean_"] = -15;
  maxima["mean_"] = 15;    
  minima["stddev_"] = 0;  
  maxima["stddev_"] = 1;     
  initialval["mean_"] = 0;
  initialerr["mean_"] = 10;    
  initialval["stddev_"] = 0.15;  
  initialerr["stddev_"] = 0.5;     
      
  lhFunction.AddSystematic(conv);
  lhFunction.SetConstraint("mean_", 0, 0.005);
  lhFunction.SetConstraint("stddev_", 0.15, 0.01);
  */

  ObsSet  obsSet(0);
  Convolution* conv = new Convolution("conv");
  GaussianERes* gaus = new GaussianERes(0.07,"gaus"); 
  gaus->RenameParameter("eres_0","eres_");
  conv->SetEResolution(gaus);
  conv->SetAxes(axes);
  conv->SetTransformationObs(obsSet);
  conv->SetDistributionObs(obsSet);
  conv->Construct();

  std::cout<<"\n \n ------------ finished constructing convolution -----------\n "<<std::endl;
  
  minima["eres_"] = 0;  
  maxima["eres_"] = 1;     
  initialval["eres_"] = 0.07;  
  initialerr["eres_"] = 0.01;     
      
  lhFunction.AddSystematic(conv);
  lhFunction.SetConstraint("stddev_", 0.07, 0.01);

  return 0.;
    
    /*
  for (int i = 0; i < numPdfs; i++){
    OscPromptE(UnOscfile,tempFile ,d21fix,s12fix,0.,reactorDistances[i]);   //oscillate each input reactor KE pdfs -> outputs a pruned ntuple with an oscillated EPrompt positron spectrum called tempFile, to be used below to fill a reactorPdf
    
    BinnedED *reactorPdf = new BinnedED(reactorNames[i],axes);
    reactorPdf->SetObservables(0);
    
    ROOTNtuple reactorNtp(tempFile, "nt");
    for(size_t j = 0; j < reactorNtp.GetNEntries(); j++)
      reactorPdf->Fill(reactorNtp.GetEntry(j));
    //std::cout<<"pre osc Int: "<<unoscPdfint<<std::endl;   //integral of input unoscillated pdfs
    //int postoscint = (int)(reactorPdf->Integral());
    //std::cout<<"post osc Int: "<<oscPdfint<<std::endl;  //integral of pdf after KE->EPrompt + oscillation
    double osc_loss = oscPdfint/(double)unoscPdfint;      // calculate loss in pdf integral
    std::cout<<" osc loss factor: "<<osc_loss<<std::endl;
    reactorPdf->Normalise();

    // Setting optimisation limits
    sprintf(name,"%s_norm",reactorNames[i].c_str());
    minima[name] = Normmin;
    maxima[name] = Normmax;
    
    double sigma = (normconstraintsuncert[i]*normconstraints[i]);
    
    initialval[name] = normconstraints[i];//(rand*(maxima[name]-minima[name]))+minima[name];
    initialerr[name] = sigma;//0.5*initialval[name];
    
    lhFunction.AddDist(*reactorPdf);
    //double constraint = normconstraints[i] * osc_loss; //vector of unoscillated norm contraints made in main() Jeff you can probs improve this
    //lhFunction.SetConstraint(name, constraint, (normconstraintsuncert[i]*constraint));  // I just randomly picked 10% for contraint uncertainty 
    lhFunction.SetConstraint(name, normconstraints[i], sigma);
    std::cout<<reactorNames[i]<<"   m: "<<normconstraints[i]<<"    sigma: "<<sigma<<"\n"<<std::endl;
  }
  
  ////////////
  // 4. Fit //
  ////////////
  Minuit min;
  min.SetMethod("Migrad");
  min.SetMaxCalls(10000000);
  min.SetTolerance(0.001);
  min.SetMinima(minima);
  min.SetMaxima(maxima);
  min.SetInitialValues(initialval);
  min.SetInitialErrors(initialerr);
  
  /////////////////////////////////////////////
  ////////        Fit Result        ///////////
  /////////////////////////////////////////////
  FitResult fitResult = min.Optimise(&lhFunction);
  ParameterDict bestFit = fitResult.GetBestFit();

  bool fitValid = fitResult.GetValid();
  if (fitValid){
    fitResult.Print();
    badfit = false;
    std::cout<<"goodfit!"<<std::endl;
  }else{
    std::cout<<"INVALID FIT!! \n Trying Again:"<<std::endl;
    badfit = true;
  }
  
  lhFunction.SetParameters(bestFit);
  double lhval = (-1.)*lhFunction.Evaluate();//(-1)*lhFunction.Evaluate();
  
  printf("\n-------------------------------------------------------------------------------------\n");
  std::cout<<"LH val at best fit: "<<lhval<<std::endl;
  printf("-------------------------------------------------------------------------------------\n");
  */


  // Gaus for results:
  /*Convolution* fitconv = new Convolution("fitconv");
  Gaussian* fitgaus = new Gaussian(bestFit.at("mean_"),bestFit.at("stddev_"),"fitgaus"); 
  fitconv->SetFunction(fitgaus);
  fitconv->SetAxes(axes);
  fitconv->SetTransformationObs(obsSet);
  fitconv->SetDistributionObs(obsSet);
  fitconv->Construct();
  
  /*
  //If Want to plot for testing purposes:  
  BinnedED* Result = new BinnedED("Result",axes);
  BinnedED TotalResult("TotalResult",axes);
  BinnedED ConstrainedTotalResult("ConstrainedTotalResult",axes);
  TotalResult.SetObservables(0);
  ConstrainedTotalResult.SetObservables(0);
  ROOTMultiPlot* Plot = new ROOTMultiPlot;
  
  TPaveText pt(0.75,0.35,1.0,0.65,"NDC");
  for (int i = 0; i< numPdfs; i++){    
    OscPromptE(UnOscfile,tempFile ,d21fix,s12fix,0.0215,reactorDistances[i]);
    
    BinnedED *reactorPdf = new BinnedED(reactorNames[i],axes);
    reactorPdf->SetObservables(0);
    
    ROOTNtuple reactorNtp(tempFile, "nt");
    for(size_t j = 0; j < reactorNtp.GetNEntries(); j++)
      reactorPdf->Fill(reactorNtp.GetEntry(j));
    
    //std::cout<<"pre osc Int: "<<unoscPdfint<<std::endl;   //integral of input unoscillated pdfs
    //std::cout<<"post osc Int: "<<oscPdfint<<std::endl;  //integral of pdf after KE->EPrompt + oscillation
    double osc_loss = oscPdfint/(double)unoscPdfint;      // calculate loss in pdf integral
    //std::cout<<" osc loss factor: "<<osc_loss<<std::endl;
    reactorPdf->Normalise();
    
    //Result->Add(fitconv(*reactorPdf),1);
    Result->Add(*reactorPdf,1);
    Result->Normalise();
      
    sprintf(name,"%s_norm",reactorNames[i].c_str());
    std::cout<<"Best Fit Norm: "<<bestFit.at(name)<<std::endl;
    Result->Scale(bestFit.at(name));
    
    TotalResult.Add(*Result,1);
    
    Result->Empty();

    //pt.AddText(Form("norm = %.5f" ,bestFit[name]));
    std::stringstream Name;
    Name<<reactorNames[i]<<"_norm = "<<bestFit[name];
    pt.AddText(Name.str().c_str());

    OscPromptE(UnOscfile,tempFile ,d21fix,s12fix,0.0215,reactorDistances[i]);
    
    BinnedED *ConstrainedreactorPdf = new BinnedED(reactorNames[i],axes);
    ConstrainedreactorPdf->SetObservables(0);
    
    ROOTNtuple ConstrainedreactorNtp(tempFile, "nt");
    for(size_t j = 0; j < ConstrainedreactorNtp.GetNEntries(); j++)
      ConstrainedreactorPdf->Fill(ConstrainedreactorNtp.GetEntry(j));
    ConstrainedreactorPdf->Normalise();

    std::cout<<reactorNames[i]<<"                "<<normconstraints[i]<<std::endl;
    ConstrainedreactorPdf->Scale(normconstraints[i]);

    ConstrainedTotalResult.Add(*ConstrainedreactorPdf,1);
  }
  
  pt.SetFillColor(kWhite);
  pt.SetShadowColor(kWhite);
 
  TH1D DataHist;
  TH1D FitHist;
  TH1D ConstrainedHist;
  
  DataHist = DistTools::ToTH1D(dataSetPdf);
  FitHist = DistTools::ToTH1D(TotalResult);
  ConstrainedHist = DistTools::ToTH1D(ConstrainedTotalResult);

  std::cout<<"Data Int: "<<DataHist.Integral()<<std::endl;
  std::cout<<"Best Fit Pdf Int: "<<FitHist.Integral()<<std::endl;
  std::cout<<"Constrained No Fit Pdf Int: "<<FitHist.Integral()<<std::endl;
  
  TH1D FullFit("FullFit","",FitHist.GetNbinsX(),FitHist.GetXaxis()->GetXmin(),FitHist.GetXaxis()->GetXmax());
  FullFit.Add(&FitHist);
  TH1D ConstrainedFullFit("ConstrainedFullFit","",ConstrainedHist.GetNbinsX(),ConstrainedHist.GetXaxis()->GetXmin(),ConstrainedHist.GetXaxis()->GetXmax());
  ConstrainedFullFit.Add(&ConstrainedHist);
  DataHist.Sumw2();

  TLegend* leg = new TLegend(0.75,0.8,1.0,1.0);
  leg->AddEntry(&DataHist,"Data","lf");
  leg->AddEntry(&FitHist,"Fit Result","lf");
  leg->AddEntry(&ConstrainedHist,"Constrained NoFit Result","lf");

  TCanvas* c1 = new TCanvas("c1");
  c1->cd();  
  DataHist.SetTitle("Data to Fit");  
  DataHist.GetYaxis()->SetTitle(Form("Counts"));
  DataHist.Draw("same");
  FitHist.SetLineColor(kRed);
  FitHist.SetLineWidth(3);
  FitHist.Draw("same");//FitHist.Draw("same e");
  ConstrainedHist.SetLineColor(kGreen);
  ConstrainedHist.SetLineWidth(3);
  ConstrainedHist.Draw("same");  //ConstrainedHist.Draw("same e");
  leg->Draw();
  
  pt.Draw();
  c1->cd();
  
  TFile * fitout = new TFile("TestFit.root","RECREATE");
  c1->Write();
    
  fitout->Close();
  */    

  //return lhval;
  
}

int main(int argc, char *argv[]) {

  int FakeDataArgc = 14 + 1;
  int NoFakeDataArgc = 12 + 1;

  if (argc < NoFakeDataArgc) {
    std::cout<<"12 (14) arguments expected: \n 1: location of UnOsc pruned Ntuple file (needed for oscillating of prompt energy event spectrum) \n 2: filename of rootfile containting constraints and Data (no .root!) \n 3: location/filename for reactor info file \n \n 4: s12 \n 5: d21 \n 6: s12bin \n 7: d21bin \n \n 8: Emin \n 9: Emax \n 10: numbins \n 11: result text file name/loc \n 12: temporary root file to fill ntuples (can delete afterwards) \n \n (13): Fake data at d21 \n (14):Fake data at s12"<<std::endl;
  }
  else{
    int Normmin = 0;
    int Normmax = 1000;
    
    std::stringstream data_constraintsstream;
    const std::string &data_constraintsFile= argv[2];
    
    if (argc == NoFakeDataArgc)
      data_constraintsstream<<data_constraintsFile;
    else if(argc == FakeDataArgc){
      double FakeDatad21 = atof(argv[13]);
      double FakeDatas12 = atof(argv[14]);
      
      data_constraintsstream<<data_constraintsFile<<"ds21_"<<FakeDatad21<<"s12_"<<FakeDatas12;
    }
    else{
      std::cout<<"Wrong number of arguments!!"<<std::endl;
      return 0;
    }	

    TFile *fin = TFile::Open((data_constraintsstream.str()+".root").c_str(),"READ");
    
    //Data
    TNtuple *datacontentNT = (TNtuple*) fin->Get("databincontents");
    
    float Count;
    datacontentNT->SetBranchAddress("counts", &Count);
    //const std::vector<double> Datacontent;
    std::vector<double> Datacontent;
    
    std::cout<<"\n  Digitised KL Paper binc contents (2nd time): "<<std::endl;
    for(int i = 0; i < datacontentNT->GetEntries(); i++){
      datacontentNT->GetEntry(i);
      std::cout<<"i="<<i<<"  "<<Count<<std::endl;
      Datacontent.push_back(Count);
    }
    
    double Emin = atof(argv[8]);
    double Emax = atof(argv[9]);
    int numbins = atoi(argv[10]);
    ObsSet dataRep(0);
    AxisCollection axes;
    axes.AddAxis(BinAxis("EPrompt (MeV)", Emin, Emax, numbins));
    
    BinnedED dataSetPdf("dataSetPdf",axes);
    dataSetPdf.SetObservables(dataRep);

    dataSetPdf.SetBinContents(Datacontent);

    //norm constraints
    TNtuple* normconstraintsNT = (TNtuple*) fin->Get("normconstraints");
    std::vector<double> normconstraints;
    std::vector<double> normconstraintsuncert;
    float norm;
    float normuncert;
    normconstraintsNT->SetBranchAddress("norm", &norm);
    normconstraintsNT->SetBranchAddress("normuncert", &normuncert);
    for(int i = 0; i < normconstraintsNT->GetEntries(); i++){
      normconstraintsNT->GetEntry(i);
      normconstraints.push_back(norm);
      normconstraintsuncert.push_back(normuncert);
    }
    
    
    //Fit, Initialisation of lhfunction
    const std::string &UnOscfile = argv[1];
    const std::string &infoFile = argv[3];
    
    std::vector<std::string> reactorNames;
    std::vector<double> distances;
    std::vector<std::string> reactorTypes;
    std::vector<int> nCores;
    std::vector<double> powers;
    std::cout<<"starting filling info....... \n"<<std::endl;
    readInfoFile(infoFile, reactorNames, distances, reactorTypes, nCores, powers); // get reactor information
    std::cout<<"\n finished filling info"<<std::endl;
    int numPdfs = reactorNames.size();

    std::cout<<"Number of reactors: "<<numPdfs<<std::endl;
    
    for (int i = 0; i< numPdfs; i++){
      std::cout<<"Reactor name: "<<reactorNames[i]<<" Distance:"<<distances[i]<<std::endl;
    }

    double s12 = atof(argv[4]);
    double d21 = atof(argv[5]);

    const std::string &tempFile = argv[12];

    std::cout<<"s12: "<<s12<<"      d21: "<<d21<<std::endl;
    double LHval;
    int fitattempts = 0;
    badfit = true;
    while (badfit){
      fitattempts += 1;
      std::cout<<"fit attempt: "<<fitattempts<<std::endl;
      LHval = LHFit(UnOscfile,tempFile,axes,dataSetPdf,numPdfs,reactorNames,distances,d21,s12,normconstraints,Normmin,Normmax,normconstraintsuncert);
    }
    
    //TH2D *h_2d = (TH2D*) fin->Get("h_2d");
    int bin_i = atoi(argv[6]);
    int bin_j = atoi(argv[7]);

    std::stringstream sstxt;
    const std::string &textfilename = argv[11];
    //sstxt<<textfilename<<"binx_"<<bin_i<<"biny_"<<bin_j<<".txt"
    ofstream textfile;
    textfile.open (textfilename.c_str(),std::ios_base::app);
    textfile <<"binx_"<<bin_i<<"biny_"<<bin_j<<"LHval_"<<LHval<<"\n";
    textfile.close();
    //h_2d->SetBinContent(bin_i,bin_j,LHval);

    //TH1D DataHist;
    //DataHist = DistTools::ToTH1D(dataSetPdf);
    //DataHist.SetName("DataHist");

    //std::cout<<" saving h2 file: "<<data_constraintsstream.str()<<"_.root"<<std::endl;
    
    //TFile *fileOut = new TFile((data_constraintsstream.str()+"_.root").c_str(),"RECREATE");

    //h_2d->Write();
    //datacontentNT->Write();
    //DataHist.Write();
    //normconstraintsNT->Write();
    fin->Close();
    //fileOut->Close();
    
    //const char *command1 = ("rm " + data_constraintsstream.str() + ".root").c_str();
    //const char *command2 = ("mv " + data_constraintsstream.str() + "_.root " + data_constraintsstream.str() + ".root").c_str();
    //system(command1);
    //system(command2);

  }
}
