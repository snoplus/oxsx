// A fit in energy for signal and a background
#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include <Rand.h>
#include <fstream>
#include <iostream>

#include <TFile.h>
#include <TCanvas.h>
#include <ROOTNtuple.h>
#include <TRandom3.h>
#include <TH1D.h>
#include <TLegend.h>

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
#include <GaussianERes.h>
#include <ParameterDict.h>
#include <ContainerTools.hpp>
#include <NuOsc.h>
#include <SurvProb.h>
#include "AntinuUtils.cpp"
#include "../util/oscillate_util_eres.cpp"

int main(int argc, char *argv[]) {
//int main() {
  /*if (argc != 25){
    std::cout<<"Error: 24 arguments expected."<<std::endl;
      return 1; // return>0 indicates error code
      }
  else{*/
    std::string data_path = "/data/snoplus/PrunedBiPos/output/combined_out_z_257669_259063_norm.root";
    std::string mc_path = "/data/snoplus/PrunedBiPos/output/combined_mc_out_Bipo214_-13p5_-23_-98p5_55_335_115_alpha_-13_-43_-650_58_25_17_alphabirk_0p076_norm.root";
    std::string outname = "out_scaled_nhit_r5000.root";

    double e_scaling_estimate = 1.0904;//1.175;
    //double e_scaling_estimate_sigma = 0.05;

    std::string hist_name = "h_nhit_Bi_r5000_norm";//"h_nhit_Bi_r6000_norm";
    printf("Begin--------------------------------------\n");
    
    TFile * f = new TFile(mc_path.c_str());
    //f->ls();
    TH1D* h_mc_nhit_Bi = (TH1D*)f->Get(hist_name.c_str());

    const size_t nhit_min = h_mc_nhit_Bi->GetXaxis()->GetXmin();//500;
    const size_t nhit_max = h_mc_nhit_Bi->GetXaxis()->GetXmax();//900;
    const size_t n_bins = h_mc_nhit_Bi->GetXaxis()->GetNbins();//(nhit_max-nhit_min)/8.;

    ObsSet  obsSet(0);

    AxisCollection axes;
    axes.AddAxis(BinAxis("nhits", nhit_min, nhit_max, n_bins/4));

    BinnedED *nhit_mc_pdf = new BinnedED("mc", axes);
    nhit_mc_pdf->SetObservables(0);
    
    for(size_t j = 0; j < h_mc_nhit_Bi->GetXaxis()->GetNbins(); j++){
      double centre = h_mc_nhit_Bi->GetXaxis()->GetBinCenter(j+1);
      double content = h_mc_nhit_Bi->GetBinContent(j+1);
      nhit_mc_pdf->Fill(centre,content);
    }

    //nhit_mc_pdf->Normalise();
    f->Close();
    
    TFile * f_data = new TFile(data_path.c_str());
    //f->ls();
    TH1D* h_data_nhit_Bi = (TH1D*)f_data->Get(hist_name.c_str());

    BinnedED *nhit_data_pdf = new BinnedED("data", axes);
    nhit_data_pdf->SetObservables(0);
      
    for(size_t j = 0; j < h_data_nhit_Bi->GetXaxis()->GetNbins() ; j++){
      double centre = h_data_nhit_Bi->GetXaxis()->GetBinCenter(j+1);
      double content = h_data_nhit_Bi->GetBinContent(j+1);
      nhit_data_pdf->Fill(centre,content);
    }

    std::cout<<"int: "<<nhit_data_pdf->Integral()<<"  "<<nhit_mc_pdf->Integral()<<std::endl;
    //scaling
    Scale* scale = new Scale("scale");
    scale->SetScaleFactor(e_scaling_estimate);
    scale->SetAxes(axes);
    scale->SetTransformationObs(obsSet);
    scale->SetDistributionObs(obsSet);
    scale->Construct();
  
    nhit_mc_pdf->Normalise();
    nhit_data_pdf->Normalise();
    
    BinnedED pdf;
    pdf = *nhit_mc_pdf;
    pdf = scale->operator()(pdf);
    pdf.Normalise();

    TFile * fout = new TFile(outname.c_str(),"RECREATE");
    
    // data set
    TH1D data_set_hist = DistTools::ToTH1D(*nhit_data_pdf);
    data_set_hist.SetName("data_set_hist");
    //data_set_hist.GetYaxis()->SetTitle("Counts");
    data_set_hist.GetXaxis()->SetTitle("nhits");
    data_set_hist.SetLineColor(kBlue);
    data_set_hist.Write();

    TH1D mc_set_hist = DistTools::ToTH1D(pdf);
    mc_set_hist.SetName("mc_set_hist");
    //mc_set_hist.GetYaxis()->SetTitle("Counts");
    mc_set_hist.GetXaxis()->SetTitle("nhits");
    mc_set_hist.SetLineColor(kRed);
    mc_set_hist.Write();

    mc_set_hist.Draw("same");
    data_set_hist.Draw("same");


    TCanvas* c_scaled = new TCanvas("c_scaled");
    data_set_hist.SetStats(0);
    c_scaled->SetLogy();
    data_set_hist.Draw("same");
    mc_set_hist.Draw("same");

    TLegend * legend = new TLegend(0.74,0.68,0.88,0.88);
    legend->AddEntry(&data_set_hist,"Data","l");
    legend->AddEntry(&mc_set_hist,"MC","l");
    legend->Draw("same");

    c_scaled->Write();
    // close unoscillated reactor file
    f_data->Close();
    fout->Close();

    printf("End--------------------------------------\n");
    return 0; // completed successfully
}
