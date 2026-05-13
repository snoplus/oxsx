/*******************************************************/
/* A specialisation of the Convolution class, for the  */
/* specific case of a Gaussian kernel.                 */
/*******************************************************/

#ifndef __OXSX_GAUSSIANCONVOLUTION__
#define __OXSX_GAUSSIANCONVOLUTION__
#include <Convolution.h>


class GaussianConvolution : public Convolution
{
public:
    // Overwrite the allowed constructor
    GaussianConvolution(const std::string &name);
    // Alternative contructor, if you want a custom CDF cutoff:
    // the number of sigmas away from the mean beyond which the CDF is set to 0 or 1 as appropriate. 5 by default.
    GaussianConvolution(const std::string &name, double cutoff_);
    // Getters/Setters for underlying Gaussian distribution's width
    // (FitComponent interface will still work, this is just a quality-of-life feature)
    double GetSigma() const;
    void SetSigma(double sigma_);
    // Get the Gaussian kernel's sigma parameter's name
    std::string GetSigmaName() const;
    // Rename the Gaussian kernel's sigma parameter's name
    // (default: stddevs_0)
    void RenameSigma(const std::string& newname_);

private:
    std::string sigma_name; // tracks name of Gaussian's sigma parameter: used internally so user doens't need to remember!

    void ConstructSubmatrix(std::vector<long long unsigned int> &column_indices, std::vector<long long unsigned int> &row_indices,
                            std::vector<double> &vals) const;
};


#endif