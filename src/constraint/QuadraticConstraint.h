/***********************************************************************************************/
/* A quadratic constraint on a fit parameter, for log likelihood and Chi-Square tests this is  */
/* equivlent to a gaussian contraint.                                                          */
/* If a second width is given, asymmetric errors are used.                                     */
/***********************************************************************************************/
#ifndef __OXSX_QUADRATIC_CONSTRAINT__
#define __OXSX_QUADRATIC_CONSTRAINT__

class QuadraticConstraint{
 public:
    QuadraticConstraint(){}
    // Symmetric constraint constructor
    QuadraticConstraint(double mean_, double width_) : 
        fMean(mean_), fWidth(width_), fWidth_lo(0), is_asym(false) {}
    // Asymmetric constraint constructor
    QuadraticConstraint(double mean_, double width_hi, double width_lo) :
        fMean(mean_), fWidth(width_hi), fWidth_lo(width_lo), is_asym(true) {}

    double Evaluate(double val_) const {
        double sigma;
        if (is_asym) { sigma = (val_ >= fMean) ? fWidth : fWidth_lo; }
        else { sigma = fWidth; }
        return (val_ - fMean) * (val_ - fMean) / (2 * sigma * sigma);
    }
    
 private:
    double fMean;
    double fWidth; //Used as both symmetric width, or asymmetric upper width. 
    double fWidth_lo;
    bool is_asym;
};
#endif
