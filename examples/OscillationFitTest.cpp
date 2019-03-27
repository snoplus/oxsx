// A simple fit in energy for signal and a background
#include <string>
#include <vector>
#include <math.h>
#include <Rand.h>
#include <fstream>

#include <TH1D.h>
#include <TH2D.h>
#include <TGraph2D.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <ROOTNtuple.h>

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
#include <Gaussian.h>
#include <ParameterDict.h>
#include <ContainerTools.hpp>
#include <NuOsc.h>
#include <SurvProb.h>
#include <TH1D.h>
#include <ROOTMultiPlot.h>
#include <TRandom3.h>
#include <TVector.h>


#include <ctime>
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

  std::string reactorName,distance,reactorType,nCore,power,norm;
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

double lhmax = -1000000000;
double lhmin = 1000000000;

void LHFit(const std::string UnOscfile, const std::string dataFile, int numPdfs, std::vector<std::string> reactorNames, std::vector<double> reactorDistances, double flux, int numbins,double Emin, double Emax,const std::string plotname){
  char name[100];
  
  TRandom3 *r1 = new TRandom3();
  r1->SetSeed(0);
  //double rand = r1->Rndm();
  ////////////////////
  // 1. Set Up PDFs //
  ////////////////////

  // Only interested in first bit of data ntuple
  ObsSet dataRep(0);

  // Set up binning
  AxisCollection axes;
  axes.AddAxis(BinAxis("ParKE", Emin, Emax, numbins));

  BinnedED dataSetPdf("dataSetPdf",axes);
  dataSetPdf.SetObservables(dataRep);

  ROOTNtuple dataNtp(dataFile, "nt");
  for(size_t i = 0; i < dataNtp.GetNEntries(); i++)
    dataSetPdf.Fill(dataNtp.GetEntry(i));  
  
  ROOTNtuple reactorNtp(UnOscfile, "nt");
  NuOsc *reactorSystematic;
  ParameterDict minima;
  ParameterDict maxima;
  ParameterDict initialval;
  ParameterDict initialerr;

  BinnedED * reactorPdf0 = new BinnedED(name,axes);
  reactorPdf0->SetObservables(0);
  for(size_t j = 0; j < reactorNtp.GetNEntries(); j++)
    reactorPdf0->Fill(reactorNtp.GetEntry(j));
  reactorPdf0->Normalise();
  
  BinnedNLLH lhFunction;
  lhFunction.SetBufferAsOverflow(true);
  int Buff = 1;
  lhFunction.SetBuffer(1,Buff,Buff);
  lhFunction.SetDataDist(dataSetPdf); // initialise withe the data set
  
  minima["d21"] = 5e-5;
  minima["s12"] = 0.2;
  maxima["d21"] = 9e-5;
  maxima["s12"] = 0.4;
  initialval["d21"] = 7.0e-5;//(r1->Uniform(minima["d21"],maxima["d21"]));//5.5e-5;//6.5e-5;
  initialval["s12"] = 0.35;//(r1->Uniform(minima["s12"],maxima["s12"]));//0.3;
  std::cout<<"\n Initial d21:  "<<initialval["d21"]<<"\n"<<std::endl;
  std::cout<<" Initial s12:  "<<initialval["s12"]<<"\n"<<std::endl;
  initialerr["d21"] = 0.1*initialval["d21"];
  initialerr["s12"] = 0.1*initialval["s12"];
  
  /*
  minima["mean"] = 2;
  minima["stdDev"] = 1;
  maxima["mean"] = 6;
  maxima["stdDev"] = 4;
  initialval["mean"] = (r1->Uniform(minima["mean"],maxima["mean"]));//6.5e-5;
  initialval["stdDev"] = (r1->Uniform(minima["stdDev"],maxima["stdDev"]));//0.3;
  std::cout<<"\n Initial mean:  "<<initialval["mean"]<<"\n"<<std::endl;
  std::cout<<" Initial stdDev:  "<<initialval["stdDev"]<<"\n"<<std::endl;
  initialerr["mean"] = 0.1*initialval["mean"];
  initialerr["stdDev"] = 0.1*initialval["stdDev"];
  */
  /*Rand::SetSeed(0);
  Gaussian gaus1(5, 1);
  BinnedED BG("BG", DistTools::ToHist(gaus1, axes));
  BG.Normalise();
  BG.Scale(200000);
  */
    
  for (int i = 0; i< numPdfs; i++){
    double rand = r1->Rndm();
    BinnedED * reactorPdf = new BinnedED(reactorNames[i],axes);
    reactorPdf->SetObservables(0);
    reactorPdf->Add(*reactorPdf0,1);
    
    sprintf(name,"%s_Survival",reactorNames[i].c_str());
    SurvProb* survprob = new SurvProb(0.1,0.1,reactorDistances[i],name);
    survprob->Setsinsqrtheta13s(0.0215);
    survprob->RenameParameter("delmsqr21_0","d21");
    survprob->RenameParameter("sinsqrtheta12_0","s12");
    sprintf(name,"%s_Systematic",reactorNames[i].c_str());
    reactorSystematic = new NuOsc(name);
    reactorSystematic->SetFunction(survprob);
    reactorSystematic->SetAxes(axes);
    reactorSystematic->SetTransformationObs(dataRep);
    reactorSystematic->SetDistributionObs(dataRep);
    
    // Setting optimisation limits
    sprintf(name,"%s_norm",reactorNames[i].c_str());
    minima[name] = 0;//Normmin;
    maxima[name] = 50000;//Normmax;
    initialval[name] = (rand*(maxima[name]-minima[name]))+minima[name];//flux*2000;
    initialerr[name] = 0.5*initialval[name];

    sprintf(name,"osc_%s",reactorNames[i].c_str());
    std::cout<<name<<std::endl;
    //lhFunction.AddSystematic(reactorSystematic,name,true);
    lhFunction.AddSystematic(reactorSystematic,name);
    lhFunction.AddDist(*reactorPdf,std::vector<std::string>(1,name),true);  
    /*
    if (i == 0){
      Convolution* convsys = new Convolution("conv_Systematic");
      Gaussian* gaussian = new Gaussian(0,1,"gaus");
      gaussian->RenameParameter("means_0","mean");
      gaussian->RenameParameter("stddevs_0","stdDev");
      convsys->SetFunction(gaussian);
      convsys->SetAxes(axes);
      convsys->SetTransformationObs(dataRep);
      convsys->SetDistributionObs(dataRep);
      convsys->Construct();

      
      sprintf(name,"osc_%s",reactorNames[i].c_str());
      std::cout<<name<<std::endl;
      lhFunction.AddSystematic(reactorSystematic,name,true);
      lhFunction.AddSystematic(convsys,"Conv1",false);
      std::vector<std::string> groups;
      groups.push_back(name);
      groups.push_back("Conv1");
      lhFunction.AddDist(*reactorPdf,groups);
    }else{
      sprintf(name,"osc_%s",reactorNames[i].c_str());
      std::cout<<name<<std::endl;
      lhFunction.AddSystematic(reactorSystematic,name,true);
      lhFunction.AddDist(*reactorPdf,std::vector<std::string>(1,name));  
    
    }
    */
  }

  lhFunction.SetConstraint("BRUCE_norm",25000,500);
  lhFunction.SetConstraint("PICKERING_norm",5600,120);
  lhFunction.SetConstraint("DARLINGTON_norm",5500,110);
  //lhFunction.SetConstraint("ReactorPdf1_norm",9450,2000);1652 D
  //lhFunction.SetConstraint("ReactorPdf2_norm",9500,2000);1728 P
  //lhFunction.SetConstraint("ReactorPdf3_norm",500,1000);

  //lhFunction.SetConstraint("d21",7.37e-5,1.6e-6);
  lhFunction.SetConstraint("s12",0.297,0.016);

  std::cout << "Built LH function " << std::endl;
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
  //min.Fix("s13");
  
  /////////////////////////////////////////////
  ////////        Fit Result        ///////////
  /////////////////////////////////////////////
  
  FitResult fitResult = min.Optimise(&lhFunction);
  ParameterDict bestFit = fitResult.GetBestFit();
  fitResult.SetPrintPrecision(6);
  //fitResult.Print();
  bool fitValid = fitResult.GetValid();
  if (fitValid)
    fitResult.Print();
  else
    std::cout<<"INVALID FIT!!"<<std::endl;
  // fitvalid flag
  
  BinnedED * Result = new BinnedED("Result",axes);
  Result->SetObservables(0);
  BinnedED TotalResult("TotalResult",axes);
  TotalResult.SetObservables(0);
  ROOTMultiPlot* Plot = new ROOTMultiPlot;
  
  TPaveText pt(0.75,0.35,1.0,0.65,"NDC");
  for (int i = 0; i< numPdfs; i++){    
    NuOsc OscResult("OscResult");
    SurvProb* survprob = new SurvProb(bestFit.at("d21"),bestFit.at("s12"),reactorDistances[i]);
    survprob->Setsinsqrtheta13s(0.0215);
    OscResult.SetFunction(survprob);
    OscResult.SetAxes(axes);
    OscResult.SetTransformationObs(dataRep);
    OscResult.SetDistributionObs(dataRep);
    OscResult.Construct();
    std::cout<<"reactorpdf0 integral: "<<reactorPdf0->Integral()<<std::endl;
    Result->Add(OscResult(*reactorPdf0),1);
    //sprintf(name,"ReactorPdf%d_norm",i);
    sprintf(name,"%s_norm",reactorNames[i].c_str());
    std::cout<<"Best Fit Norm: "<<bestFit.at(name)<<std::endl;
    Result->Scale(bestFit.at(name));

    TotalResult.Add(*Result,1);
    
    //pt.AddText(Form("norm = %.5f" ,bestFit[name]));
    std::stringstream Name;
    Name<<reactorNames[i]<<"_norm = "<<bestFit[name];
    pt.AddText(Name.str().c_str());
    Result->Empty();
  }
  
  pt.AddText(Form("#Delta m_{21} = %.6f",bestFit["d21"]));
  pt.AddText(Form("#theta_{12} = %.3f",bestFit["s12"]));
  pt.SetFillColor(kWhite);
  pt.SetShadowColor(kWhite);
 
  TH1D DataHist;
  TH1D FitHist;
  
  DataHist = DistTools::ToTH1D(dataSetPdf);
  FitHist = DistTools::ToTH1D(TotalResult);
  
  TH1D FullFit("FullFit","",FitHist.GetNbinsX(),FitHist.GetXaxis()->GetXmin(),FitHist.GetXaxis()->GetXmax());
  FullFit.Add(&FitHist);
  DataHist.Sumw2();

  TLegend* leg = new TLegend(0.75,0.8,1.0,1.0);
  leg->AddEntry(&DataHist,"Data","lf");
  leg->AddEntry(&FitHist,"Fit Result","lf");

  TCanvas* c1 = new TCanvas("c1");
  c1->cd();  
  DataHist.SetTitle("Data to Fit");  
  DataHist.GetYaxis()->SetTitle(Form("Counts"));
  DataHist.Draw();
  FitHist.SetLineColor(kRed);
  FitHist.SetLineWidth(3);
  FitHist.Draw("same e");
  leg->Draw();
  
  pt.Draw();
  c1->cd();
  
  TFile * fitout = new TFile(plotname.c_str(),"RECREATE");
  c1->Write();
  
  fitout->Close();

  /*lhFunction.SetParameters(bestFit);
  double lhval =(-1)*lhFunction.Evaluate();
  if (lhval > lhmax)
    lhmax = lhval;
  if (lhval < lhmin)
    lhmin = lhval;
  std::cout<<"LH val at best fit: "<<lhval<<std::endl;
  return lhval;*/
}

int main(int argc, char *argv[]){
  
  if (argc != 5) {
    std::cout<<"4 arguments expected: \n 1: location of UnOscfile \n 2: location/filename\
 for data spectrum to fit \n 3: location/filename for reactor info file \n 4: filename of output\
 plot"<<std::endl;
  }
  else{
    clock_t begin = clock();
    double flux = 1;
    int numbins = 30;
    
    double Emin = 2;
    double Emax = 8;
    
    std::stringstream argParser;
    const std::string &UnOscfile = argv[1];
    const std::string &dataFile = argv[2];
    const std::string &infoFile = argv[3];
    const std::string &outFile = argv[4];

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

    LHFit(UnOscfile,dataFile,numPdfs,reactorNames,distances,flux,numbins,Emin,Emax,outFile.c_str());

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout<<"elapsed time: "<<elapsed_secs<<std::endl;
    return 0;
  }
}
