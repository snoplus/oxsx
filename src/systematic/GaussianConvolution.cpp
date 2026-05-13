#include <GaussianConvolution.h>
#include <JumpPDF.h>
#include <ConditionalPDF.h>
#include <Gaussian.h>

GaussianConvolution::GaussianConvolution(const std::string &name_) :
Convolution(name_), sigma_name("stddevs_0")
{
    // Set the kernel to be a Gaussian.
    // SetFunction() clones this Gaussian object, packaging it within
    // a JumpPDF object (a kind of ConditionalPDF object).
    Gaussian gauss_tmp(1, "");
    SetFunction(&gauss_tmp);
}

GaussianConvolution::GaussianConvolution(const std::string &name_, double cutoff_) :
Convolution(name_), sigma_name("stddevs_0")
{
    // Set the kernel to be a Gaussian.
    // SetFunction() clones this Gaussian object, packaging it within
    // a JumpPDF object (a kind of ConditionalPDF object).
    Gaussian gauss_tmp(1, "");
    gauss_tmp.SetCdfCutOff(cutoff_); // set a custom CDF cutoff
    SetFunction(&gauss_tmp);
}

void GaussianConvolution::ConstructSubmatrix(std::vector<long long unsigned int> &column_indices, std::vector<long long unsigned int> &row_indices,
                                     std::vector<double> &vals) const
{
    /*
     * Construct the sub-matrix associated with the convolution in the transformed axes.
     * Returns this information, by reference, in the form of column & row indices for non-zero entries of this submatrix,
     * and their associated values. This form is what Armadillo prefers for sparse matrices!
     */
    // variables storing the axes bin information
    std::vector<double> binCentres(fSubMapAxes.GetNDimensions());
    std::vector<double> lowEdges(fSubMapAxes.GetNDimensions());
    std::vector<double> highEdges(fSubMapAxes.GetNDimensions());
    // Pre-allocate memory for the sub-matrix data
    // likely to need only a fraction of this memory!
    column_indices.reserve(fSubMapAxes.GetNBins() * fSubMapAxes.GetNBins());
    row_indices.reserve(fSubMapAxes.GetNBins() * fSubMapAxes.GetNBins());
    vals.reserve(fSubMapAxes.GetNBins() * fSubMapAxes.GetNBins());
    
    // Some tricks we are pulling here to make things go even faster:
    // - In most use cases, the binning along the smearing axis is ~equal,
    //   so the amount of smearing will become translation-invariant!
    //   Cache calculations of the smearing Integral() (which is expensive)
    //   in a map object based on the integration endpoints (the output bin edges)
    //   After the first row, we expect ~all integrals to have already been calculated!
    // - For typical smearing kernels, the contribution necessarily goes monotonically down
    //   for bins further away from the original bin. In other words, we expect this submatrix
    //   to be quite diagonal-heavy, with no random spikes far away.
    //   Also, beyond some distance we expect the contribution to go to zero.
    //   Instead of blithely scanning as usual over the full 2D grid of (origBin, destBin),
    //   we now start ~along the diagonal and work outwards - first forwards, then backwards.
    //   Whenever a zero is reached, we immediately exit the loop! 
    std::map<std::pair<double, double>, double> integral_cache;
    // Loop over all entries of the sub-matrix to determine their values
    for (long long unsigned int origBin = 0; origBin < fSubMapAxes.GetNBins(); origBin++)
    {
        // get the centre of the bin. Need to offset by this for a convolution
        fSubMapAxes.GetBinCentres(origBin, binCentres);

        // loop over the bins it can be smeared into, going forwards
        // from the bin after the diagonal first 
        for (long long unsigned int destBin = origBin+1; destBin < fSubMapAxes.GetNBins(); destBin++)
        {
            fSubMapAxes.GetBinLowEdges(destBin, lowEdges);
            fSubMapAxes.GetBinHighEdges(destBin, highEdges);

            // Calculate destination bin edges along smearing axis, relative
            // to origBin's centre
            const double xlo = lowEdges.at(0) - binCentres.at(0);
            const double xhi = lowEdges.at(0) - binCentres.at(0);
            const std::pair<double, double> edges {xlo, xhi};
            // Has an integral already been calculated for this pair of (relative) bin edges?
            const auto it = std::find_if(integral_cache.begin(), integral_cache.end(), [edges](const std::pair<const std::pair<double, double>, double>& p){ return p.first == edges; });
            double integral = 0.;
            if (it == integral_cache.end())
            {
                // Nope, need to calculate the integral (and add it to the cache)
                integral = fDist->Integral(lowEdges, highEdges, binCentres);
                integral_cache[edges] = integral;
            } else
            {
                // Yep, use that result!
                integral = it->second;
            }

            // Only bother adding to matrix if non-zero!
            // No point looking at any further destination bins in this loop
            // as we know they'll be zero too
            if (integral == 0.)
            {
                break;
            }
            column_indices.push_back(origBin);
            row_indices.push_back(destBin);
            vals.push_back(integral);
        }
        // Now do the same inner loop, but going backwards, starting from
        // the diagonal - this ensures we don't accidentally try negative bins initially
        for (long long unsigned int destBin = origBin; destBin < fSubMapAxes.GetNBins(); destBin--)
        {
            fSubMapAxes.GetBinLowEdges(destBin, lowEdges);
            fSubMapAxes.GetBinHighEdges(destBin, highEdges);

            const double xlo = lowEdges.at(0) - binCentres.at(0);
            const double xhi = lowEdges.at(0) - binCentres.at(0);
            const std::pair<double, double> edges {xlo, xhi};
            const auto it = std::find_if(integral_cache.begin(), integral_cache.end(), [edges](const std::pair<const std::pair<double, double>, double>& p){ return p.first == edges; });
            double integral = 0.;
            if (it == integral_cache.end())
            {
                integral = fDist->Integral(lowEdges, highEdges, binCentres);
                integral_cache[edges] = integral;
            } else
            {
                integral = it->second;
            }

            if (integral == 0.)
            {
                break;
            }
            column_indices.push_back(origBin);
            row_indices.push_back(destBin);
            vals.push_back(integral);
        }
    }
}

double GaussianConvolution::GetSigma() const
{
    return fDist->GetParameter(sigma_name);
}

void GaussianConvolution::SetSigma(double sigma_)
{
    fDist->SetParameter(sigma_name, sigma_);
}

std::string GaussianConvolution::GetSigmaName() const
{
    return sigma_name;
}

void GaussianConvolution::RenameSigma(const std::string& newname_)
{
    fDist->RenameParameter(sigma_name, newname_);
    sigma_name = newname_;
}