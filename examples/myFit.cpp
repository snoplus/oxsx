#ifndef __myFit__
#define __myFit__

//includes from g++
#include <vector>
#include <iostream>

//includes from oxo
#include <AxisCollection.h>
#include <BinnedED.h>
#include <Rand.h>
#include <DistTools.h>
#include <IO.h>
#include <BinnedNLLH.h>
#include <ParameterDict.h>
#include <Minuit.h> 
#include <GridSearch.h>
#include <BinnedEDManager.h>

//includes from root
#include <TMath.h>
#include <TFile.h>
#include <TH1.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLegend.h>

void Drawing(TString filename1, TString filename2, TString filename3, TString filename4, TString nbins, TString filename5 = "Nothing"){
  TFile* file1 = new TFile(filename1+".root","READ");
  TFile* file2 = new TFile(filename2+".root","READ");
  TFile* file3 = new TFile(filename3+".root","READ");
  TFile* file4 = new TFile(filename4+".root","READ");
  TFile* file5;
  if(filename5 != "Nothing"){
    file5 = new TFile(filename5+".root","READ");
  }else{
    file5 = file4;
  }

  TH1D* dat = (TH1D*) file1->Get(filename1);
  TH1D* Pa = (TH1D*) file2->Get(filename2);
  TH1D* pa = (TH1D*) file3->Get(filename3);
  TH1D* tl = (TH1D*) file4->Get(filename4);
  TH1D* elem;
  if(filename5 != "Nothing"){
    elem = (TH1D*) file5->Get(filename5);
  }else{
    elem = tl;
  }
  TFile* output = new TFile("FitOutput.root","RECREATE");

  dat->SetTitle(filename5+" Gold+Bronze Fit "+nbins+" bins");
  if(filename5 == "Nothing") dat->SetTitle("PaTl Gold+Bronze Fit "+nbins+" bins");
  dat->Write("data");
  Pa->Write("Pa234m");
  pa->Write("Pa234mFit");
  tl->Write("Tl208");
  elem->Write(filename5);
  
  TCanvas* C = new TCanvas("C","Fit",800,600);
  C->SetLogy();
  TLegend* L = new TLegend(0.6,0.7,0.9,0.9);
  dat->SetLineColor(1);
  Pa->SetLineColor(4);
  pa->SetLineColor(3);
  tl->SetLineColor(2);
  L->AddEntry(dat,"Data Gold + Bronze");
  L->AddEntry(Pa,"Pa234m");
  L->AddEntry(pa,"Pa234mFit");
  L->AddEntry(tl,"Tl208");
  dat->Draw();
  Pa->Draw("same");
  pa->Draw("same");
  tl->Draw("same");
  TH1D* tot = (TH1D*) Pa->Clone();
  tot->Add(tl);
  
  
  elem->SetLineColor(95);
  elem->Draw("same");
  L->AddEntry(elem,filename5);
  if(filename5 != "Nothing") tot->Add(elem);
  
  tot->SetLineColor(6);
  tot->Draw("same");
  tot->Write("total");
  L->AddEntry(tot,"Total");
  L->Draw("same");
  output->Write();
  C->SaveAs(filename5+"Fit"+nbins+"bins.pdf");
  if(filename5 == "Nothing") C->SaveAs("PaTlFit"+nbins+"bins.pdf");
  
  
  file1->Close();
  file2->Close();
  file3->Close();
  file4->Close();
  output->Close();

  delete file1;
  delete file2;
  delete file3;
  delete file4;
  delete output;
}



Histogram cut(Histogram &hist, double lmin, double lmax){
  Histogram cuttedH = hist;
  std::cout << "hist NBins "<< hist.GetNBins() << std::endl;
  int minbin, maxbin;
  std::vector <double> lm;
  std::vector <double> LM;
  lm.push_back(lmin);
  LM.push_back(lmax);
  minbin = hist.FindBin(lm);
  maxbin = hist.FindBin(LM);
  std::cout << "minbin = " << minbin << std::endl;
  std::cout << "maxbin = " << maxbin << std::endl;
  //std::cout << "sm " << cuttedH.GetBinContent(200) << std::endl;
  for(int i = 1; i <= hist.GetNBins(); i++){
    if(i < minbin  ||  i > maxbin){
      cuttedH.SetBinContent(i,0.01);
    }
  }
  return cuttedH;
}

Histogram Rebin(Histogram &hist, int NBins){
  
  AxisCollection ACo = hist.GetAxes();
  std::vector<std::string> axisnames;
  axisnames = ACo.GetAxisNames();
  BinAxis Ax = ACo.GetAxis(axisnames[0]);
  
  BinAxis BA(axisnames[0], Ax.GetMin(), Ax.GetMax(), NBins, Ax.GetName());
  AxisCollection AC;
  AC.AddAxis(BA);
  
  Histogram RebinH(AC);
  double sum=0;
  int j=1, MaxBins=hist.GetNBins();
  double divider = MaxBins/NBins;
  int div = std::ceil(divider);
  std::cout << "Rebinning" << std::endl;
  for(int i=0; i < MaxBins; i++){
    sum += hist.GetBinContent(i);
    if(i >= j*div || i == MaxBins-1){
      RebinH.SetBinContent(j,sum);
      j++;
      sum=0;
    }
  }
  return RebinH;
}

//INITIALIZING FITTING PARAMETERS-----------------
ParameterDict minima,maxima,scale,errors;
 


int main(int argc, char* argv[]){
  std::cout << "----------------------------------------------------------------------------------------" << std::endl;

  //RETRIEVING HISTOGRAMS FROM FILES
  if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
  }

  //OPEN FILES
  TString filename = argv[1];
  std::cout << "Filename provided: " << filename << std::endl;

  TFile *file1 = TFile::Open(filename,"READ");
  if (!file1 || file1->IsZombie()) {
    std::cerr << "Error: Could not open ROOT file " << filename << std::endl;
    return 1;
  }

  std::cout << "File opened sucessfully" << std::endl;

  //CHOOSE ON THE TERMINAL THE BINNING AND THE ELEMENT
  //ALLWAYS RUN ./program file1Data.root ... file1MC.root ... ELEM_NAME NBINS
  //---------------------------------------------------------
  int NBINS = std::stoi(argv[argc-1]);
  std::string elem;
  if (argc > 6){
    elem = argv[argc-2];
  }else{
    elem = "Nothing";
  }
  
  //////////////////////////////////////////////////////////////////////////////
  //CHOOSE HERE THE FIT RANGE
  double MIN,MAX;
  MIN = 1.22;//1.22
  MAX = 6.0;
  
  //GET DATA GOLD HISTOGRAM FOR FIT
  TH1D* DataG = (TH1D*) file1->Get("hsum_llh");
  std::cout << "data gold File obtained" << std::endl;
  Histogram dataG;
  dataG = DistTools::ToHist(*DataG);
  dataG = Rebin(dataG,NBINS);
  std::cout << "Data Gold Nbins = " << dataG.GetNBins()  << std::endl;
  
  for(int i=1; i<=dataG.GetNBins(); i++){
    if(dataG.GetBinContent(i) == 0) dataG.SetBinContent(i,0.1);
  }
  
  //--------------GET DATA BRONZE --------------------------
  //OPEN FILES
  filename = argv[2];
  std::cout << "Filename provided: " << filename << std::endl;

  TFile *file2 = TFile::Open(filename,"READ");
  if (!file2 || file2->IsZombie()) {
    std::cerr << "Error: Could not open ROOT file " << filename << std::endl;
    return 1;
  }

  //GET GOLD/BRONZE HISTOGRAM FOR FIT
  TH1D* DataB = (TH1D*) file2->Get("hsum_rho");
  std::cout << "data bronze File obtained" << std::endl;
  Histogram dataB;
  dataB = DistTools::ToHist(*DataB);
  dataB = Rebin(dataB,NBINS);
  std::cout << "Bronze Gold Nbins = " << dataB.GetNBins()  << std::endl;
  
  for(int i=1; i<=dataB.GetNBins(); i++){
    if(dataB.GetBinContent(i) == 0) dataB.SetBinContent(i,0.1);
  }

  bool b = true;//variable to control use of data bronze+gold vs only gold
  Histogram NCdataB = dataB;
  Histogram NCdata = dataG;
  if(b) NCdata.Add(NCdataB,1);
  dataB = cut(dataG,MIN,MAX);
  dataG = cut(dataG,MIN,MAX);
  std::cout<<"argc " << argc << std::endl;
  IO::SaveHistogramROOT(dataB,"DataB.root","dataBronze");
  IO::SaveHistogramROOT(dataG,"DataG.root","dataGold");
  BinnedED DATAB("bronze_data",dataB);
  BinnedED data("energy_data",dataG);
  if(b) data.Add(DATAB,1);
  std::cout << "data integral = " << data.Integral() << std::endl;
  //data.Normalise();
  BinnedNLLH logl;

  //Adding the data first
  logl.SetDataDist(data);

 
  Histogram Pa234;
  double AreaPa, AreaEl;

  Histogram NotcutPa234;
  //GET PA234 HIST FOR FIT
  if(argc > 3){

    filename = argv[3];
    std::cout << "Filename provided: " << filename << std::endl;
    TFile *file3 = TFile::Open(filename,"READ");
    if (!file3 || file3->IsZombie()) {
      std::cerr << "Error: Could not open ROOT file " << filename << std::endl;
      return 1;
    }

    TH1D* Pa234m = (TH1D*) file3->Get("Energy_LLHCut");
    std::cout << "Pa234m File obtained" << std::endl;
   
    Pa234 = DistTools::ToHist(*Pa234m);
    Pa234 = Rebin(Pa234,NBINS);
    std::cout << "Pa234m Nbins = " << Pa234.GetNBins()  << std::endl;


    for(int i=0; i < Pa234.GetNBins(); i++){
      if(Pa234.GetBinContent(i) == 0) Pa234.SetBinContent(i,0.01);
    }
    
    IO::SaveHistogramROOT(Pa234,"Pa234m.root","Pa234m");
    std::cout << "Setted to non zero" << std::endl;
    NotcutPa234 = Pa234;
    Pa234 = cut(Pa234, MIN, MAX);
    AreaPa = Pa234.Integral(); //to later get the whole histogram and not the fit region
    //Adding the PDFs
    Pa234.Normalise();
    
    BinnedED Pa234M("Pa234M",Pa234);
    Pa234M.Normalise();

    std::cout << "Pa234M integral = " << Pa234M.Integral() << std::endl;
    logl.AddPdf(Pa234M);

    //SETTING MINIMA/MAXIMA FOR THE FIT
    maxima["Pa234M"] = 2.e8;
    minima["Pa234M"] = 3.;
    scale["Pa234M"] = 0.5;
    errors["Pa234M"] = 1.e-4;
  }

  Histogram Tl208;

  if(argc > 4){
    filename = argv[4];
    std::cout << "Filename provided: " << filename << std::endl;
    TFile *file4 = TFile::Open(filename,"READ");
    if (!file4 || file4->IsZombie()) {
      std::cerr << "Error: Could not open ROOT file " << filename << std::endl;
      return 1;
    }

    TH1D* Tl = (TH1D*) file4->Get("Energy_LLHCut");
    std::cout << "Tl208 File obtained" << std::endl;
   
    Tl208 = DistTools::ToHist(*Tl);
    Tl208 = Rebin(Tl208,NBINS);
    std::cout << "Tl208 Nbins = " << Tl208.GetNBins()  << std::endl;
  
    for(int i=1; i<=Tl208.GetNBins(); i++){
      if(Tl208.GetBinContent(i) == 0) Tl208.SetBinContent(i,0.01);
    }

    IO::SaveHistogramROOT(Tl208,"Tl208.root","Tl208");
    std::cout << "Setted to non zero" << std::endl;
    //Adding the PDFs
    Tl208.Normalise();
    
    BinnedED TL208("TL208",Tl208);
    TL208.Normalise();

    std::cout << "Tl208 integral = " << TL208.Integral() << std::endl;
    logl.AddPdf(TL208);

    //SETTING MINIMA/MAXIMA FOR THE FIT
    maxima["TL208"] = 2.e8;
    minima["TL208"] = 3.e-1;
    scale["TL208"] = 0.005;
    errors["TL208"] = 1.e-4;
  }


  //ELEMENT HYPOTHYSIS
  Histogram Elem, NCElem;
  if(argc > 6){
    filename = argv[5];
    std::cout << "Filename provided: " << filename << std::endl;
    TFile *file5 = TFile::Open(filename,"READ");
    if (!file5 || file5->IsZombie()) {
      std::cerr << "Error: Could not open ROOT file " << filename << std::endl;
      return 1;
    }

    TH1D* el = (TH1D*) file5->Get("Energy_LLHCut");
    std::cout << "element File obtained" << std::endl;
    Elem = DistTools::ToHist(*el);
    Elem = Rebin(Elem,NBINS);
    std::cout << "Elem Nbins = " << Elem.GetNBins()  << std::endl;
  
    for(int i=1; i<=Elem.GetNBins(); i++){
      if(Elem.GetBinContent(i) == 0) Elem.SetBinContent(i,0.01);
    }

    IO::SaveHistogramROOT(Elem,elem+"ff.root",elem);
    std::cout << "Setted to non zero" << std::endl;
    //Adding the PDFs
    NCElem = Elem;
    Elem = cut(Elem,MIN,MAX);
    AreaEl = NCElem.Integral();
    Elem.Normalise();
    
    BinnedED ELEM("ELEM",Elem);
    ELEM.Normalise();

    std::cout << "Elem integral = " << ELEM.Integral() << std::endl;
    logl.AddPdf(ELEM);

    //SETTING MINIMA/MAXIMA FOR THE FIT
    maxima["ELEM"] = 2.e8;
    minima["ELEM"] = 3.;
    scale["ELEM"] = 0.05;
    errors["ELEM"] = 1.e-4;
  }
 

  //ParameterDict steps;
  //  steps["Pa234M"] = 0.1;

  
  //logl.SetDebugMode(true);


  //logl.Evaluate();

  ////////////
  // 4. Fit //
  ////////////
  
  //-------GRIDSEARCH-------------------------------------------
  /*GridSearch gSearch;

  // Set optimisation parameters.
  gSearch.SetMinima(minima);
  gSearch.SetMaxima(maxima);
  gSearch.SetStepSizes(steps);

  */
  //-------MINUIT-----------------------------------
  Minuit gSearch;
  std::cout << "Fitting0" << std::endl;
  gSearch.SetMethod("Migrad");
  gSearch.SetMaxCalls(10000);

  std::cout << "setting Fit Parameters" << std::endl;
  
  gSearch.SetInitialValues(scale);
  gSearch.SetInitialErrors(errors);
  
  gSearch.SetMinima(minima);
  gSearch.SetMaxima(maxima);
  
  std::cout << "maxtol" << std::endl;
  gSearch.SetTolerance(0.0001);
  
  std::cout << "Fitting" << std::endl;
  //IO::SaveHistogramROOT(Pa234,"Pa234mFitF.root","Pafit234");
  
  FitResult result = gSearch.Optimise(&logl);

  std::cout << "Fit Done" << std::endl;
  ParameterDict fit = result.GetBestFit();
  result.Print();
  result.SaveAs("simpleFit_result.txt");
  
  
  //----------GRAPH PLOTTING---------------
  //Pa234 = Pa234M.GetHistogram();
  IO::SaveHistogramROOT(NCdata,"NCdata.root","NCdata");
  IO::SaveHistogramROOT(dataG,"dataG.root","dataG");
  IO::SaveHistogramROOT(Pa234,"PaRaw.root","PaRaw");

    
  std::cout << "Pa234m = "<< fit["Pa234M"] << std::endl;
  std::cout << "Elem = "<< fit["ELEM"] << std::endl;
  Pa234.Scale(fit["Pa234M"]);
  Elem.Scale(fit["ELEM"]);
  std::cout << "-------------------------------" << std::endl;
  std::cout << "Num of Data Events = " << NCdata.Integral() << std::endl;

  double areaPa = Pa234.Integral();
  NotcutPa234.Scale(areaPa/AreaPa);
  std::cout << "Num of Pa234 Events = " << NotcutPa234.Integral() << std::endl;

  Tl208.Scale(fit["TL208"]);
  std::cout << "Num of Tl208 Events = " << Tl208.Integral() << std::endl;

  double areaEl = Elem.Integral();
  NCElem.Scale(areaEl/AreaEl);
  std::cout << "Num of " << elem << " Events = " << NCElem.Integral() << std::endl;

  IO::SaveHistogramROOT(NotcutPa234,"Pa234m.root","Pa234m");
  IO::SaveHistogramROOT(Pa234,"Pa234mFit.root","Pa234mFit");
  IO::SaveHistogramROOT(Tl208,"Tl208.root","Tl208");
  if(argc>6) IO::SaveHistogramROOT(Elem,elem+"Fit.root",elem);
  if(argc>6) IO::SaveHistogramROOT(NCElem,elem+".root",elem);
  std::cout << "Drawing " << std::endl;
  Drawing("NCdata","Pa234m","Pa234mFit","Tl208", argv[argc-1], elem);
  
  std::cout << "Finished " << std::endl;
  return 0;
}


#endif
