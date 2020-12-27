// A fit in energy for signal and a background
#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include <Rand.h>
#include <fstream>
#include <iostream>

#include <TLegend.h>
#include <TCanvas.h>
#include <ROOTNtuple.h>
#include <TRandom3.h>
#include <TH1D.h>

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
#include "AntinuUtils.cpp"
#include "../util/oscillate_util.cpp"

void Make_Fake_Data(BinnedED &data_set_pdf, const std::string &spectrum_phwr_unosc_filepath,
    const std::string &spectrum_pwr_unosc_filepath,
    const std::string &spectrum_geo_uraniumthorium_unosc_filepath,
    const std::string &spectrum_bkg_alphan_unosc_filepath_1,
    const std::string &spectrum_bkg_alphan_unosc_filepath_2,
    std::vector<std::string> &reactor_names, std::vector<std::string> &reactor_types,
    std::vector<Double_t> &distances,
    std::vector<Double_t> &constraint_means, std::vector<Double_t> &constraint_sigmas,
    std::vector<int> &constrain_on_vec,
    TFile *file_out,
    Double_t param_d21, Double_t param_s12, Double_t param_s13,
    const double e_min, const double e_max, const size_t n_bins,
    const double flux_data, const double mc_scale_factor,
    const double geo_uth_custom_flux, const double bkg_alphan_custom_flux,
    BinnedED &reactors_data_set_pdf, BinnedED &geos_data_set_pdf, BinnedED &alpha_n_data_set_pdf,
    BinnedED &alpha_n_1st_data_set_pdf,  BinnedED &alpha_n_2nd_data_set_pdf,
    const std::string &split_pdf_params_file ){

    printf("Beginning Fake data prep--------------------------------------\n");
    printf("Fake Data:: del_21:%.9f, sin2_12:%.7f, sin2_13:%.7f\n", param_d21, param_s12, param_s13);

    char name[1000];
    TRandom3 *random_generator = new TRandom3();
    const ULong64_t n_pdf = reactor_names.size();

    // set up binning
    AxisCollection axes;
    axes.AddAxis(BinAxis("ev_prompt_fit", e_min, e_max, n_bins));


    BinnedED **reactor_unosc_pdf = new BinnedED*[n_pdf];
    BinnedED **reactor_osc_pdf = new BinnedED*[n_pdf];

    Double_t constraint_osc_mean_total = 0.;
    //Double_t data_set_pdf_integral = data_set_pdf.Integral();

    bool geos_included = false;
    bool alphans_included = false;

    double reactor_total = 0.;
    double geos_total = 0.;
    double alphan_13c_total = 0.;
  
    for (ULong64_t i = 0; i < n_pdf; i++){

        sprintf(name, "%s_unosc", reactor_names[i].c_str());
        reactor_unosc_pdf[i] = new BinnedED(name, axes);
        reactor_unosc_pdf[i]->SetObservables(0);
        reactor_osc_pdf[i] = new BinnedED(reactor_names[i], axes);
        reactor_osc_pdf[i]->SetObservables(0);

        bool apply_oscillation = false;
        bool is_further_reactors = false;
        bool split_alpha_n_pdf_1st = false;
        bool split_alpha_n_pdf_2nd = false;
  
        if ((reactor_types[i]=="PWR")||(reactor_types[i]=="BWR")){
          sprintf(name, "%s", spectrum_pwr_unosc_filepath.c_str());
          apply_oscillation = true;
        }else if (reactor_types[i]=="PHWR"){
          sprintf(name, "%s", spectrum_phwr_unosc_filepath.c_str());
          apply_oscillation = true;
        }else if (reactor_types[i]=="further_reactors"){
          sprintf(name, "%s", spectrum_pwr_unosc_filepath.c_str());
          is_further_reactors = true;
        }else if (reactor_names[i]=="uraniumthorium"){
          sprintf(name, "%s", spectrum_geo_uraniumthorium_unosc_filepath.c_str());
          geos_included = true;
          //apply_geo_oscillation = true;
        }else if (reactor_names[i]=="alpha_n_lab_13c"){
          sprintf(name, "%s", spectrum_bkg_alphan_unosc_filepath_1.c_str());
          alphans_included = true;
        }else if (reactor_names[i]=="alpha_n_lab_13c_1st"){
          sprintf(name, "%s", spectrum_bkg_alphan_unosc_filepath_1.c_str());
          alphans_included = true;
          split_alpha_n_pdf_1st = true;
        }else if (reactor_names[i]=="alpha_n_lab_13c_2nd"){
          sprintf(name, "%s", spectrum_bkg_alphan_unosc_filepath_2.c_str());
          alphans_included = true;
          split_alpha_n_pdf_2nd = true;
        }else{
          printf("Throw: Reactor doesn't match any loaded type...\n");
          exit(0); // throw std::exception(); //continue;
        }

        // load unoscillated reactor file (to oscillate, and to plot)
        //ROOTNtuple reactor_unosc_ntp(spectrum_unosc_filepath.c_str(), "nt"); // this would be made easier if this worked for specific branches!!
        //std::cout<<reactor_names[i]<<" pdf name: "<<name<<std::endl;
        TFile *f_in = new TFile(name);
        file_out->cd();
        TTree *reactor_unosc_ntp = (TTree*)f_in->Get("nt");
        TNtuple *reactor_osc_ntp = new TNtuple("nt", "Oscillated Prompt Energy", "ev_fit_energy_p1");

        // oscillate tree
        if (apply_oscillation)
          ntOscillate_pruned(reactor_unosc_ntp, reactor_osc_ntp, param_d21, param_s12, param_s13, distances[i]);
        else if (is_further_reactors) //|| apply_geo_oscillation)
          ntOscillate_pruned_geo(reactor_unosc_ntp, reactor_osc_ntp, param_s12);
        else
          ntNoOscillate_pruned(reactor_unosc_ntp, reactor_osc_ntp);

	// reset branch addresses after oscillating in function (otherwise crash before setting again below..)
        reactor_unosc_ntp->SetBranchStatus("*", 0);
        reactor_unosc_ntp->SetBranchStatus("ev_fit_energy_p1", 1); // (re-enable all branches in use)

        // fill unoscillated pdf
        Double_t ev_unosc_energy_p1;
        reactor_unosc_ntp->SetBranchAddress("ev_fit_energy_p1", &ev_unosc_energy_p1);
        for(size_t j = 0; j < reactor_unosc_ntp->GetEntries(); j++){
            reactor_unosc_ntp->GetEntry(j);
            reactor_unosc_pdf[i]->Fill(ev_unosc_energy_p1);
        }

        // fill oscillated pdf
        Float_t ev_osc_energy_p1;
        reactor_osc_ntp->SetBranchAddress("ev_fit_energy_p1", &ev_osc_energy_p1);
        for(size_t j = 0; j < reactor_osc_ntp->GetEntries(); j++){
            reactor_osc_ntp->GetEntry(j);
            reactor_osc_pdf[i]->Fill(ev_osc_energy_p1);
        }

        // close unoscillated reactor file
        f_in->Close();

        if (apply_oscillation || is_further_reactors){
          // work out total oscillated integral of constraints
          Double_t normalisation_unosc = reactor_unosc_pdf[i]->Integral();
          Double_t normalisation_reactor = reactor_osc_pdf[i]->Integral();
          Double_t osc_loss = normalisation_reactor/normalisation_unosc;

          Double_t constraint_osc_mean = constraint_means[i]*osc_loss*mc_scale_factor;
          Double_t constraint_osc_sigma = (constraint_sigmas[i]/constraint_means[i])*constraint_osc_mean;
          reactor_osc_pdf[i]->Normalise(); 
          reactor_osc_pdf[i]->Scale(constraint_osc_mean);
	  
          data_set_pdf.Add(*reactor_osc_pdf[i]);
          //for plotting only:
          reactors_data_set_pdf.Add(*reactor_osc_pdf[i]);

          printf("  added reactor %d/%d: %s, osc_survival: %.3f, norm_constraint: %.3f err: %.3f data_int:%.3f\n", i+1, n_pdf, reactor_names[i].c_str(), osc_loss, reactor_osc_pdf[i]->Integral(), constraint_osc_sigma, data_set_pdf.Integral());
	  
          reactor_total += constraint_osc_mean;
        }else if (geos_included){
          Double_t normalisation_unosc = reactor_unosc_pdf[i]->Integral();
          Double_t normalisation_reactor = reactor_osc_pdf[i]->Integral();
          Double_t osc_loss = normalisation_reactor/normalisation_unosc;

          Double_t constraint_osc_mean = constraint_means[i]*osc_loss*mc_scale_factor;
          
          // geo custom flux = coefficient times 1 years predicted geos
          reactor_osc_pdf[i]->Normalise();
          reactor_osc_pdf[i]->Scale(constraint_osc_mean*geo_uth_custom_flux);//geo_uth_custom_flux/(double)flux_data);

          data_set_pdf.Add(*reactor_osc_pdf[i]);
          //for plotting only:
          geos_data_set_pdf.Add(*reactor_osc_pdf[i]);

          printf("  added geo %d/%d: %s, osc_survival: %.3f, num_evs: %.3f, data_int:%.3f\n", i+1, n_pdf, reactor_names[i].c_str(), osc_loss, reactor_osc_pdf[i]->Integral(), data_set_pdf.Integral());

          geos_total += constraint_osc_mean*geo_uth_custom_flux;
	
        }else if (alphans_included){
          Double_t normalisation_unosc = reactor_unosc_pdf[i]->Integral();
          Double_t normalisation_reactor = reactor_osc_pdf[i]->Integral();
          Double_t osc_loss = normalisation_reactor/normalisation_unosc;

          reactor_osc_pdf[i]->Normalise(); // currently custom flux number = total events in livetime
          /////// split alphan pdfs
          /*
          double split_norm_modification = 1.;
          if (split_alpha_n_pdf_1st || split_alpha_n_pdf_2nd){
            TFile *f_split_info = new TFile(split_pdf_params_file.c_str(),"READ");
            file_out->cd(); // switch to output file (for ntuple to use)
            TTree *split_info_tree = (TTree*)f_split_info->Get("alphan_split_pdf_params");
            //std::cout<<reactor_names[i]<<"  split_norm_modification: "<<split_norm_modification<<std::endl;
            split_norm_modification = CalculateSplitPdfTails(split_info_tree, reactor_osc_pdf[i], split_alpha_n_pdf_1st, split_alpha_n_pdf_2nd,e_min,e_max);
          }
          */
          Double_t constraint_osc_mean = constraint_means[i]*osc_loss*mc_scale_factor*bkg_alphan_custom_flux;//*split_norm_modification;
          
          reactor_osc_pdf[i]->Scale(constraint_osc_mean);

          data_set_pdf.Add(*reactor_osc_pdf[i]);
          //for plotting only:
          alpha_n_data_set_pdf.Add(*reactor_osc_pdf[i]);
          if (split_alpha_n_pdf_1st)
            alpha_n_1st_data_set_pdf.Add(*reactor_osc_pdf[i]);
          else if (split_alpha_n_pdf_2nd)
            alpha_n_2nd_data_set_pdf.Add(*reactor_osc_pdf[i]);
          
          printf("  added alpha-n %d/%d: %s, osc_survival: %.3f, num_evs: %.3f, data_int:%.3f\n", i+1, n_pdf, reactor_names[i].c_str(), osc_loss, reactor_osc_pdf[i]->Integral(), data_set_pdf.Integral());

          alphan_13c_total += constraint_osc_mean;
        }
	
    }

    //databincontents = data_set_pdf.GetBinContents();
    printf("End data prep-------------------------------------\n");
    std::cout<<"\nReactor ev total: "<<reactor_total<<std::endl;
    std::cout<<"Geo ev total: "<<geos_total<<std::endl;
    std::cout<<"alpha-n 13C ev total: "<<alphan_13c_total<<std::endl;
    std::cout<<""<<std::endl;
}

int main(int argc, char *argv[]) {

    if (argc != 20){
        std::cout<<"Error: 19 arguments expected."<<std::endl;
        return 1; // return>0 indicates error code
    }
    else{
        const std::string &data_path = argv[1];
        const std::string &info_file = argv[2];
        const std::string &spectrum_phwr_unosc_filepath = argv[3];
        const std::string &spectrum_pwr_unosc_filepath = argv[4];
        const std::string &spectrum_geo_uraniumthorium_unosc_filepath = argv[5];
        const std::string &spectrum_bkg_alphan_unosc_filepath_1 = argv[6];
        const std::string &spectrum_bkg_alphan_unosc_filepath_2 = argv[7];
        const std::string &constraints_info_file = argv[8];
        const double s12 = atof(argv[9]);
        const double d21 = atof(argv[10]);
        const double s13 = atof(argv[11]);
        const double flux_data = atof(argv[12]);
        const double mc_scale_factor = atof(argv[13]);
        const double geo_uth_custom_flux = atof(argv[14]);
        const double bkg_alphan_custom_flux = atof(argv[15]);
        const double e_min = atof(argv[16]);
        const double e_max = atof(argv[17]);
        const size_t n_bins = atoi(argv[18]);
        const std::string &split_pdf_params_file = argv[19];
        printf("Begin--------------------------------------\n");

        // read in reactor information
        std::vector<std::string> reactor_names;
        std::vector<Double_t> distances;
        std::vector<std::string> reactor_types;
        std::vector<ULong64_t> n_cores;
        std::vector<Double_t> powers;
        std::vector<Double_t> power_errs;
        readInfoFile(info_file, reactor_names, distances, reactor_types, n_cores, powers, power_errs);
        
        // read in constraint information
        std::vector<Double_t> constraint_means;
        std::vector<Double_t> constraint_mean_errs;
        std::vector<Double_t> constraint_sigmas;
        std::vector<Double_t> constraint_sigma_errs;
        std::vector<int> constrain_on_vec;

        // read constraint info for each reactor in the info file (one at time to ensure they match correctly)
        for (size_t i=0; i<(size_t)reactor_names.size(); i++){
            double fit_mean, fit_mean_err, fit_sigma, fit_sigma_err;
            int constrain_on;
            readConstraintsInfoFile(constraints_info_file, reactor_names[i].c_str(), fit_mean, fit_mean_err, fit_sigma, fit_sigma_err, constrain_on);
            constraint_means.push_back(fit_mean);
            constraint_mean_errs.push_back(fit_mean_err);
            constraint_sigmas.push_back(fit_sigma);
            constraint_sigma_errs.push_back(fit_sigma_err);
            constrain_on_vec.push_back(constrain_on);
        }

        for (size_t i=0; i<(size_t)reactor_names.size(); i++)
          printf("i:%llu, reactor_name:%s, fit_mean: %.3f, fit_sigma: %.3f, contraining: %i\n", i, reactor_names[i].c_str(), constraint_means[i], constraint_sigmas[i], constrain_on_vec[i]);

        const ULong64_t n_pdf = reactor_names.size();

        AxisCollection axes;
        axes.AddAxis(BinAxis("ev_prompt_fit", e_min, e_max, n_bins));

        ObsSet data_rep(0);
        BinnedED data_set_pdf("data_set_pdf", axes);
        data_set_pdf.SetObservables(data_rep);

        // for plotting separate contributions:
        BinnedED reactors_data_set_pdf("reactors_data_set_pdf", axes);
        reactors_data_set_pdf.SetObservables(data_rep);
        BinnedED geos_data_set_pdf("geos_data_set_pdf", axes);
        geos_data_set_pdf.SetObservables(data_rep);
        BinnedED alpha_n_data_set_pdf("alpha_n_data_set_pdf", axes);
        alpha_n_data_set_pdf.SetObservables(data_rep);
        BinnedED alpha_n_1st_data_set_pdf("alpha_n_1st_data_set_pdf", axes);
        alpha_n_1st_data_set_pdf.SetObservables(data_rep);
        BinnedED alpha_n_2nd_data_set_pdf("alpha_n_2nd_data_set_pdf", axes);
        alpha_n_2nd_data_set_pdf.SetObservables(data_rep);
	
        std::vector<double> databincontents;

        TFile *file_out = new TFile(data_path.c_str(), "RECREATE");
	
        // initialise data
        Make_Fake_Data(data_set_pdf, spectrum_phwr_unosc_filepath,
                       spectrum_pwr_unosc_filepath,
                       spectrum_geo_uraniumthorium_unosc_filepath,
                       spectrum_bkg_alphan_unosc_filepath_1,
                       spectrum_bkg_alphan_unosc_filepath_2,
                       reactor_names, reactor_types,
                       distances,
                       constraint_means, constraint_sigmas,
                       constrain_on_vec,
                       file_out,
                       d21, s12, s13,
                       e_min, e_max, n_bins,
                       flux_data, mc_scale_factor, geo_uth_custom_flux,
                       bkg_alphan_custom_flux, reactors_data_set_pdf,
                       geos_data_set_pdf, alpha_n_data_set_pdf,
                       alpha_n_1st_data_set_pdf, alpha_n_2nd_data_set_pdf,
                       split_pdf_params_file);

        databincontents = data_set_pdf.GetBinContents();

        TNtuple* datacontentNT = new TNtuple("databincontents","databincontents","counts");
        for (int i = 0; i< databincontents.size(); i++)
          datacontentNT->Fill(databincontents[i]);
    
        ////save objects to file
        printf("Save objects to file...\n");

        datacontentNT->Write();

        TH1D DataDist = DistTools::ToTH1D(data_set_pdf);
        DataDist.SetName("h_data");
        DataDist.SetStats(0);
        DataDist.SetLineColor(1);
        DataDist.SetLineWidth(2);
        DataDist.GetXaxis()->SetTitle("Reconstructed Energy");
        DataDist.Write();

        TCanvas* c_all = new TCanvas("c_all","c_all");
        // plotting:
        DataDist.Draw("same");

        TH1D ReactorsDataDist = DistTools::ToTH1D(reactors_data_set_pdf);
        ReactorsDataDist.SetName("h_data_reactors");
        ReactorsDataDist.SetStats(0);
        ReactorsDataDist.SetLineColor(2);
        ReactorsDataDist.SetLineWidth(2);
        ReactorsDataDist.GetXaxis()->SetTitle("Reconstructed Energy");
        ReactorsDataDist.Draw("same");
        
        TH1D GeosDataDist = DistTools::ToTH1D(geos_data_set_pdf);
        GeosDataDist.SetName("h_data_geos");
        GeosDataDist.SetStats(0);
        GeosDataDist.SetLineColor(8);
        GeosDataDist.GetXaxis()->SetTitle("Reconstructed Energy");
        if (GeosDataDist.Integral() > 0)
          GeosDataDist.Draw("same");

        TH1D Alphan1stDataDist = DistTools::ToTH1D(alpha_n_1st_data_set_pdf);
        Alphan1stDataDist.SetName("h_data_alphan1st");
        Alphan1stDataDist.SetStats(0);
        Alphan1stDataDist.SetLineColor(7);
        Alphan1stDataDist.GetXaxis()->SetTitle("Reconstructed Energy");
        if (Alphan1stDataDist.Integral() > 0)
          Alphan1stDataDist.Draw("same");

        TH1D Alphan2ndDataDist = DistTools::ToTH1D(alpha_n_2nd_data_set_pdf);
        Alphan2ndDataDist.SetName("h_data_alphan2nd");
        Alphan2ndDataDist.SetStats(0);
        Alphan2ndDataDist.SetLineColor(8);
        Alphan2ndDataDist.GetXaxis()->SetTitle("Reconstructed Energy");
        if (Alphan2ndDataDist.Integral() > 0)
          Alphan2ndDataDist.Draw("same");

        TH1D AlphanDataDist = DistTools::ToTH1D(alpha_n_data_set_pdf);
        AlphanDataDist.SetName("h_data_alphan");
        AlphanDataDist.SetStats(0);
        AlphanDataDist.SetLineColor(6);
        AlphanDataDist.GetXaxis()->SetTitle("Reconstructed Energy");
        if (Alphan1stDataDist.Integral() == 0 &&\
            Alphan2ndDataDist.Integral() == 0 &&\
            AlphanDataDist.Integral() > 0)
          AlphanDataDist.Draw("same");

        TLegend* leg = new TLegend(0.55,0.55,0.9,0.9);
        leg->AddEntry(&DataDist,"Total Signal+BG","l");
        leg->AddEntry(&ReactorsDataDist,"Reactors (Osc)","l");
        leg->AddEntry(&GeosDataDist,"Geoneutrinos","l");
        leg->AddEntry(&AlphanDataDist,"Alpha-n 13C","l");
        
        leg->Draw("same");

        c_all->Write();

        // close output file
        file_out->Close();

        printf("End--------------------------------------\n");
        return 0; // completed successfully
    }
}
