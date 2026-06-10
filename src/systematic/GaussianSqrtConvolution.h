/****************************************************************/
/* A specialisation of the Convolution class, for the           */
/* specific case of a Gaussian kernel with a width that         */
/* varies proportional to the square root of the smearing axis. */
/* Useful in particular for smearing along in energy, say.      */
/****************************************************************/

#ifndef __OXSX_GAUSSIANSQRTCONVOLUTION__
#define __OXSX_GAUSSIANSQRTCONVOLUTION__
#include <Convolution.h>


class GaussianSqrtConvolution : public Convolution
{
public:
    // Overwrite the allowed constructor
    GaussianSqrtConvolution(const std::string &name);
    // Alternative contructor, if you want a custom CDF cutoff:
    // the number of sigmas away from the mean beyond which the CDF is set to 0 or 1 as appropriate. 5 by default.
    GaussianSqrtConvolution(const std::string &name, double cutoff_);
    // Getters/Setters for underlying SquareRootScale's width parameter -
    // this is the proportionality constant which multiplies the sqrt scaling
    // (FitComponent interface will still work, this is just a quality-of-life feature)
    double GetSigma() const;
    void SetSigma(double sigma_);
    // Get the kernel's width parameter's name
    std::string GetSigmaName() const;
    // Rename the kernel's width parameter's name
    // (default: "grad")
    void RenameSigma(const std::string& newname_);

private:
    std::string fSigmaName; // tracks name of kernel's width parameter: used internally so user doesn't need to remember!

    void ConstructSubmatrix(std::vector<long long unsigned int> &column_indices, std::vector<long long unsigned int> &row_indices,
                            std::vector<double> &vals) const;
};


#endif