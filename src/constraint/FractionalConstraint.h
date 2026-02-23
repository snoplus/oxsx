/***********************************************************************************************/
/* A quadratic constraint on the fractional difference (a-b)/(a+b) of two fit parameters.      */
/* For log likelihood and Chi-Square tests this is equivalent to a gaussian contraint on       */
/* that quantity.                                                                              */
/***********************************************************************************************/
#ifndef __OXSX_FRACTIONAL_CONSTRAINT__
#define __OXSX_FRACTIONAL_CONSTRAINT__
#include <Exceptions.h>

class FractionalConstraint
{
public:
    // Default constructor
    FractionalConstraint() {}
    // Main constructor
    FractionalConstraint(double mean_, double sigma_) : fFractionalMean(mean_), fFractionalSigma(sigma_) {}

    // Evaluate constraint at a particular pair of values
    double Evaluate(double val_1, double val_2) const
    {

        const double fraction = (val_1 - val_2) / (val_1 + val_2);

        return (fraction - fFractionalMean) * (fraction - fFractionalMean) / (2 * fFractionalSigma * fFractionalSigma);
    }

private:
    double fFractionalMean;
    double fFractionalSigma;
};
#endif
