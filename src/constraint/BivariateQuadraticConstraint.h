/***********************************************************************************************/
/* A bivariate quadratic constraint on two fit parameters.                                     */
/* For log likelihood and Chi-Square tests this is equivalent                                  */
/* to a bivariate gaussian contraint.                                                          */
/* The object is designed for two fit parameters having a correlated constraint.               */
/***********************************************************************************************/
#ifndef __OXSX_BIVARIATE_QUADRATIC_CONSTRAINT__
#define __OXSX_BIVARIATE_QUADRATIC_CONSTRAINT__
#include <Exceptions.h>

class BivariateQuadraticConstraint
{
public:
    // Default constructor
    BivariateQuadraticConstraint() {}
    // Main constructor
    BivariateQuadraticConstraint(double mean_1, double sigma_1, double mean_2, double sigma_2, double correlation) : fMean1(mean_1), fSigma1(sigma_1), fMean2(mean_2), fSigma2(sigma_2), fCorr(correlation)
    {
        if (fCorr == 1.)
        {
            throw ValueError("BivariateQuadraticConstraint: correlation cannot be equal to 1!");
        }
    }
    // Evaluate constraint at a particular pair of values
    double Evaluate(double val_1, double val_2) const
    {
        const double z1 = (val_1 - fMean1) * (val_1 - fMean1) / (2 * fSigma1 * fSigma1);
        const double z2 = (val_2 - fMean2) * (val_2 - fMean2) / (2 * fSigma2 * fSigma2);
        const double z12 = fCorr * (val_1 - fMean1) * (val_2 - fMean2) / (fSigma1 * fSigma2);
        return (z1 - z12 + z2) / (1. - fCorr * fCorr);
    }

private:
    double fMean1;
    double fSigma1;
    double fMean2;
    double fSigma2;
    double fCorr;
};
#endif
