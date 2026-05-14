/*
MCMC.cpp

This example shows how MCMC is run using OXO, in a realistic analysis situation.
This example is also useful for providing a way of profiling how fast OXO runs
in these sorts of circumstances, which is valuable for us to know.
*/
// C++ Headers
#include <iostream>
// ROOT Headers
#include <TTree.h>
#include <TFile.h>
// OXO Headers
#include <AxisCollection.h>
#include <BinnedNLLH.h>
#include <ParameterDict.h>
#include <FitResult.h>
#include <MetropolisSampler.h>
#include <MCMC.h>
#include <MCMCSamples.h>
#include <Gaussian.h>
#include <DistTools.h>
#include <Scale.h>
#include <VaryingCDF.h>
#include <SquareRootScale.h>
#include <Convolution.h>
#include <GaussianConvolution.h>
#include <GaussianSqrtConvolution.h>
// global consts
constexpr double N_SIGNAL = 100.; // Truth signal rate in dataset
constexpr double N_BACK = 300.; // Truth background rate in dataset
constexpr double SIGMA_BACK = 50.; // Constraint uncertainty on background normalisation
constexpr double SIGMA_ESCALE = 0.01;
constexpr double SIGMA_ESMEAR = 0.05;
constexpr double SIGMA_RSCALE = 0.01;
constexpr double SIGMA_RSMEAR = 0.02;
constexpr size_t N_STEPS = 100;
const std::string outfilename_root = "mcmc_example_output.root";

std::vector<BinnedED> load_mc(const AxisCollection& ax, const std::vector<std::string>& observables) {
    /*
     * Create a set of fake MC PDFs for use in this macro. You'd want to actually load
     * in the real MC at this step.
     */
    // Make a signal PDF: Gaussian in energy, flat in r3
    BinnedED signal("signal", ax);
    signal.SetObservables(observables);
    const Gaussian signal_gaus(3, 1);
    AxisCollection ax_e;
    ax_e.AddAxis(ax.GetAxis("energy"));
    BinnedED signal_energy("signal_energy", DistTools::ToHist(signal_gaus, ax_e));
    std::cout << "\tGenerating fake signal PDF\n";
    for (size_t r_idx = 0; r_idx < ax.GetAxis("r3").GetNBins(); r_idx++) {
        for (size_t e_idx = 0; e_idx < ax.GetAxis("energy").GetNBins(); e_idx++) {
            const double binval = signal_energy.GetBinContent(e_idx);
            const std::vector<size_t> idxs {e_idx, r_idx};
            signal.SetBinContent(signal.FlattenIndices(idxs), binval);
        }
    }
    signal.Normalise();
    // Now make a background PDF: falling in energy, rising in r3
    BinnedED background("background", ax);
    background.SetObservables(observables);
    std::cout << "\tGenerating fake background PDF\n";
    for (size_t r_idx = 0; r_idx < ax.GetAxis("r3").GetNBins(); r_idx++) {
        for (size_t e_idx = 0; e_idx < ax.GetAxis("energy").GetNBins(); e_idx++) {
            const double binval = 1. + 5.*(ax.GetAxis("energy").GetBinCentre(e_idx) - 5.)*(ax.GetAxis("energy").GetBinCentre(e_idx) - 5.) + 10.*(ax.GetAxis("r3").GetBinCentre(r_idx))*(ax.GetAxis("r3").GetBinCentre(r_idx));
            const std::vector<size_t> idxs {e_idx, r_idx};
            background.SetBinContent(background.FlattenIndices(idxs), binval);
        }
    }
    background.Normalise();

    const std::vector<BinnedED> pdfs {signal, background};
    return pdfs;
}

BinnedED load_data(const AxisCollection& ax, const std::vector<std::string>& observables,
                   const std::vector<BinnedED>& mc_pdfs) {
    /*
     * Create a fake dataset for use in this macro. You'd actually want to
     * load your real data in this step.
     */
    std::cout << "\tGenerating (fake) data BinnedED\n";
    // Create empty data
    BinnedED data("test", ax);
    // Add observables list
    data.SetObservables(observables);
    // Now let's fill in the data histogram with a weighted sum of the MC PDFs!
    data.Add(mc_pdfs.at(0), N_SIGNAL);
    data.Add(mc_pdfs.at(1), N_BACK);

    return data;
}

BinnedNLLH create_lh(const BinnedED& data_hist, const std::vector<BinnedED>& mc_pdfs) {
    /*
     * Do the first major set of steps for setting up the likelihood function.
     */
    BinnedNLLH lh_function;
    // Settings for buffer bins which we'll need when we add systematics
    lh_function.SetBufferAsOverflow(false);
    lh_function.SetBuffer("energy", 1, 1);
    lh_function.SetBuffer("r3", 1, 1);
    // Add data distribution
    lh_function.SetDataDist(data_hist);
    // Add MC PDFs
    // It'll use the "INDIRECT" NormFittingStatus for each by default,
    // Which is what we'll want when it comes to systematics here.
    lh_function.AddPdfs(mc_pdfs);

    // Add simple background normalisation constraint (from some sideband, say)
    lh_function.SetConstraint("background", N_BACK, SIGMA_BACK);

    return lh_function;
}

void add_systematics(BinnedNLLH& lh_function) {
    /*
     * Add the systematics to the likelihood, with constraints:
     *  - Energy scale
     *  - Energy smear
     *  - Radial scale
     *  - Radial smear
     */
    // 1: Add energy scale
    Scale* escale = new Scale("energy_scale");
    escale->RenameParameter("scaleFactor", "energy_scale_factor");
    escale->SetScaleFactor(1.0);
    escale->SetAxes(lh_function.GetDataDist().GetAxes());
    escale->SetDistributionObs(ObsSet(std::vector<std::string>({"energy", "r3"})));
    escale->SetTransformationObs(ObsSet("energy"));
    escale->Construct();

    lh_function.AddSystematic(escale);
    lh_function.SetConstraint("energy_scale_factor", 1.0, SIGMA_ESCALE);

    // 2: Add energy smear (less complicated than it used to be!)
    //   A Gaussian kernel...
    // VaryingCDF* smearer = new VaryingCDF("energy_smear");
    // Gaussian* gaus = new Gaussian(0, 0.01, "e_gaus"); // temp sigma value
    // gaus->RenameParameter("means_0", "mean");
    // gaus->RenameParameter("stddevs_0", "sigma");
    // smearer->SetKernel(gaus);
    // //   ...With a width scaling by the square root of the energy...
    // SquareRootScale* smear_sigma_func = new SquareRootScale("e_smear_sigma_func");
    // smear_sigma_func->RenameParameter("grad", "energy_smear_param");
    // smear_sigma_func->SetGradient(0.03); // temp gradient value
    // smearer->SetDependance("sigma", smear_sigma_func);
    // //  ...which smears the PDFs along the energy axis.
    // Convolution* esmear = new Convolution("esmear");
    // esmear->SetConditionalPDF(smearer);

    GaussianSqrtConvolution* esmear = new GaussianSqrtConvolution("esmear");
    esmear->RenameSigma("energy_smear_param");
    esmear->SetAxes(lh_function.GetDataDist().GetAxes());
    esmear->SetDistributionObs(ObsSet(std::vector<std::string>({"energy", "r3"})));
    esmear->SetTransformationObs(ObsSet("energy"));
    esmear->Construct();

    lh_function.AddSystematic(esmear);
    lh_function.SetConstraint("energy_smear_param", 0., SIGMA_ESMEAR);

    // 3: Add radial scale
    Scale* rscale = new Scale("radial_scale");
    rscale->RenameParameter("scaleFactor", "radial_scale_factor");
    rscale->SetScaleFactor(1.0);
    rscale->SetAxes(lh_function.GetDataDist().GetAxes());
    rscale->SetDistributionObs(ObsSet(std::vector<std::string>({"energy", "r3"})));
    rscale->SetTransformationObs(ObsSet("r3"));
    rscale->Construct();

    lh_function.AddSystematic(rscale);
    lh_function.SetConstraint("radial_scale_factor", 1.0, SIGMA_RSCALE);

    // 4: Add radial smear
    GaussianConvolution* rsmear = new GaussianConvolution("radial_smear");
    rsmear->RenameSigma("radial_smear_sigma");
    rsmear->SetAxes(lh_function.GetDataDist().GetAxes());
    rsmear->SetDistributionObs(ObsSet(std::vector<std::string>({"energy", "r3"})));
    rsmear->SetTransformationObs(ObsSet("r3"));
    rsmear->Construct();

    lh_function.AddSystematic(rsmear);
    lh_function.SetConstraint("radial_smear_sigma", 0., SIGMA_RSMEAR);
}

std::pair<double, double> calc_scale_bounds(BinnedNLLH& lh_function, 
    const AxisCollection& axisCollection, const std::string& axis_str) {
    /*
     * The buffer bins on an axis for the pdfs set a bound
     * on the allowable scale factors that can be used before we're
     * in danger of zero-probability bins.
     */
    // Get the number of buffer bins at the bottom and top
    const auto num_buffers = lh_function.GetBuffer(axis_str);
    //
    const auto axis = axisCollection.GetAxis(axis_str);
    const size_t num_bins = axis.GetNBins();
    const double bottom_low = axis.GetBinLowEdge(0);
    const double bottom_high = axis.GetBinLowEdge(num_buffers.first);
    const double top_low = axis.GetBinHighEdge(num_bins - 1 - num_buffers.second);
    const double top_high = axis.GetBinHighEdge(num_bins - 1);
    //
    const double min_escale = top_low/top_high;
    const double max_escale = bottom_high/bottom_low;
    return std::pair<double, double>{min_escale, max_escale};
}

void setup_metaparams(ParameterDict& minima, ParameterDict& maxima,
                      ParameterDict& initial_vals, ParameterDict& initial_err,
                      const std::pair<double, double>& escale_bounds,
                      const std::pair<double, double>& rscale_bounds) {
    /*
     * Set up the minima, maxima, initial values, and step sizes
     * for all parameters in the fit.
     */
    // signal
    minima["signal"] = 0.;
    maxima["signal"] = 5.*N_SIGNAL;
    initial_vals["signal"] = N_SIGNAL;
    initial_err["signal"] = sqrt(N_SIGNAL);
    // background
    minima["background"] = 0.;
    maxima["background"] = 5.*N_BACK;
    initial_vals["background"] = N_BACK;
    initial_err["background"] = SIGMA_BACK;
    // energy_scale_factor
    minima["energy_scale_factor"] = escale_bounds.first;
    maxima["energy_scale_factor"] = escale_bounds.second;
    initial_vals["energy_scale_factor"] = 1.0;
    initial_err["energy_scale_factor"] = SIGMA_ESCALE;
    // energy_smear_param
    minima["energy_smear_param"] = 0.;
    maxima["energy_smear_param"] = 5.*SIGMA_ESMEAR;
    initial_vals["energy_smear_param"] = 0.00001;
    initial_err["energy_smear_param"] = SIGMA_ESMEAR;
    // radial_scale_factor
    minima["radial_scale_factor"] = rscale_bounds.first;
    maxima["radial_scale_factor"] = rscale_bounds.second;
    initial_vals["radial_scale_factor"] = 1.0;
    initial_err["radial_scale_factor"] = SIGMA_RSCALE;
    // radial_smear_sigma
    minima["radial_smear_sigma"] = 0.;
    maxima["radial_smear_sigma"] = 5.*SIGMA_RSMEAR;
    initial_vals["radial_smear_sigma"] = 0.00001;
    initial_err["radial_smear_sigma"] = SIGMA_RSMEAR;
}

TTree* run_mcmc(BinnedNLLH& lh_function, const ParameterDict& minima, const ParameterDict& maxima,
                      const ParameterDict& initial_vals, const ParameterDict& initial_err) {
    /*
     * Given a fully-prepared likelihood function and meta-parameters, 
     * run the MCMC optimiser, using the Metropolis sampling algorithm.
     */
    // Set up Metropolis sampls
    MetropolisSampler sampler;
    sampler.SetSigmas(initial_err);
    // Set up MCMC object, which uses the Metropolis sampler
    MCMC mcmc(sampler);
    mcmc.SetMaxIter(N_STEPS);
    mcmc.SetMinima(minima);
    mcmc.SetMaxima(maxima);
    mcmc.SetInitialTrial(initial_vals);
    mcmc.SetTestStatLogged(true); // We are using the LOG-likelihood test statistic, not the likelihood
    mcmc.SetFlipSign(true); // We are using the NEGATIVE LLH
    mcmc.SetSaveChain(true); //If you want to save the full chain state!

    // Run MCMC!!
    // Unlike a normal optimiser, the result of an MCMC process is not a
    // single best-fit position, but the set of sampled positions themselves.
    std::cout << "...SETUP COMPLETE. RUNNING MCMC NOW!\n";
    auto result = mcmc.Optimise(&lh_function);
    std::cout << "MCMC CHAIN COMPLETE:\n";
    result.SetPrintPrecision(9);
    result.Print();
    // Get sampled positions, in TTree form.
    MCMCSamples mcmcSamples = mcmc.GetSamples();
    TTree* t = mcmcSamples.GetChain();
    return t;
}


int main() {
    std::cout << "SETTING UP MCMC...\n";
    // First - define binned axes for analysis
    // 3D in energy, radius**3, and ITR, 40*10*25 = 10000 bins total
    AxisCollection ax;
    ax.AddAxis(BinAxis("energy", 1, 5, 40));
    ax.AddAxis(BinAxis("r3", 0, 1, 10));
    const std::vector<std::string> observables {"energy", "r3"};
    // "Load" in MC PDFs for fit (just making this up!)
    const std::vector<BinnedED> mc_pdfs = load_mc(ax, observables);
    // "Load" in data for fit (just making this up!)
    const BinnedED data_hist = load_data(ax, observables, mc_pdfs);
    // Create likelihood function
    std::cout << "\tMC and data hist setup.\n";
    BinnedNLLH lh_function = create_lh(data_hist, mc_pdfs);
    add_systematics(lh_function);
    // Set up fit meta-parameters
    const auto escale_bounds = calc_scale_bounds(lh_function, ax, "energy");
    const auto rscale_bounds = calc_scale_bounds(lh_function, ax, "r3");
    ParameterDict minima, maxima, initial_vals, initial_err;
    setup_metaparams(minima, maxima, initial_vals, initial_err, escale_bounds, rscale_bounds);
    // Now build the MCMC object, and run the MCMC!
    auto outfile_samples = new TFile(outfilename_root.c_str(), "RECREATE");
    TTree* chain_output = run_mcmc(lh_function, minima, maxima, initial_vals, initial_err);
    chain_output->Write();
    outfile_samples->Close();

    return 0;
}