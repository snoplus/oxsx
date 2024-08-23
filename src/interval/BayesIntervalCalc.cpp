#include <BayesIntervalCalc.h>
#include <Exceptions.h>
#include <Histogram.h>
#include <gsl/gsl_sf_gamma.h> // imcomplete gamma function
#include <iostream>
#include <cmath>

double
BayesIntervalCalc::UpperBound(Histogram posterior_, double cl_)
{
    // Makes no sense for an empty histogram
    if (!posterior_.Integral())
        throw ValueError("BayesIntervalCalc::Empty histogram passed!");

    if (cl_ <= 0)
        throw ValueError(Formatter() << "BayesIntervalCalc:: cl = " << cl_
                                     << " , must be >0!");

    // only works for 1D histograms currently
    if (posterior_.GetNDims() != 1)
        throw DimensionError("BayesIntervalCal", 1, posterior_.GetNDims(),
                             "Only implemented for 1D, marginalise?");

    posterior_.Normalise();

    // Integrate from the minimum to x until the total probability exceeds cl_
    // find the first bin <b> for which the integral  <0>-><b> > cl
    size_t critBin = 0;
    double sum = 0;
    while (sum < cl_)
    {
        sum = 0;
        for (size_t i = 0; i < critBin; i++)
        {
            sum += posterior_.GetBinContent(i);
        }

        critBin++;
    }

    // now interpolate between the bins
    double upperEdge = posterior_.GetAxes().GetBinHighEdge(critBin, 0);
    double lowerEdge = posterior_.GetAxes().GetBinLowEdge(critBin, 0);
    double content = posterior_.GetBinContent(critBin);

    return upperEdge - (upperEdge - lowerEdge) * (sum - cl_) / content;
}

double
BayesIntervalCalc::OneSidedUpperInterval(double expectedCounts_, int observedCounts_, double upperEdge_)
{
    return (gsl_sf_gamma_inc_P(observedCounts_ + 1, expectedCounts_ + upperEdge_) - gsl_sf_gamma_inc_P(observedCounts_ + 1, expectedCounts_)) / (1 - gsl_sf_gamma_inc_P(observedCounts_ + 1, expectedCounts_));
}

double
BayesIntervalCalc::UpperBound(double expectedCounts_, int observedCounts_, double cl_, double tolerance_)
{
    double currentCL = 0;
    double limit = 0;
    while (std::abs(currentCL - cl_) > tolerance_)
    {
        limit += tolerance_;
        currentCL = OneSidedUpperInterval(expectedCounts_, observedCounts_, limit);
    }
    return limit;
}
