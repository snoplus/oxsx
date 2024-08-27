#ifndef __OXSX_GAUSSIAN_FITTER__
#define __OXSX_GAUSSIAN_FITTER__
#include <vector>
#include <string>
#include <ParameterDict.h>
#include <set>
#include <iostream>

class Gaussian;

class GaussianFitter
{
public:
    GaussianFitter(Gaussian *gaus, const size_t &nDims);
    GaussianFitter(Gaussian *gaus, const std::vector<std::string> &, const std::vector<std::string> &);

    void SetParameter(const std::string &name_, double value);
    double GetParameter(const std::string &name_) const;

    void SetParameters(const ParameterDict &);
    ParameterDict GetParameters() const;
    size_t GetParameterCount() const;

    std::set<std::string> GetParameterNames() const;
    void RenameParameter(const std::string &old_, const std::string &new_);

    std::vector<std::string> GetMeanNames() const;
    std::vector<std::string> GetStdDevNames() const;

private:
    Gaussian *fOrignalFunc;
    std::vector<std::string> fMeansNames;
    std::vector<std::string> fStdDevsNames;
};
#endif
