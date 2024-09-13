/***********************************************************************************************/
/* A quadratic constraint on the ratio of two fit parameters.                                  */
/* For log likelihood and Chi-Square tests this is equivalent to a gaussian contraint on the   */
/* ratio.                                                                                      */
/***********************************************************************************************/
#ifndef __OXSX_RATIO_CONSTRAINT__
#define __OXSX_RATIO_CONSTRAINT__
#include <Exceptions.h>

class RatioConstraint
{
public:
    // Default constructor
    RatioConstraint() {}
    // Main constructor
    RatioConstraint(double mean_, double sigma_ ) : fRatioMean(mean_), fRatioSigma(sigma_) { }

    // Evaluate constraint at a particular pair of values
    double Evaluate(double val_1, double val_2) const
    {
        const double ratio = val_1 / val_2;

	return (ratio - fRatioMean) * (ratio - fRatioMean) / (2 * fRatioSigma * fRatioSigma);
    }

private:
    double fRatioMean;
    double fRatioSigma;
};
#endif
