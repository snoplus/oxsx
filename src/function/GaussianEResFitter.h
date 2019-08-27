#ifndef __OXSX_GAUSSIANERES_FITTER__
#define __OXSX_GAUSSIANERES_FITTER__
#include <vector>      
#include <string>
#include <ParameterDict.h>
#include <set>
#include <iostream>

class GaussianERes;

class GaussianEResFitter{
public:

    GaussianEResFitter(GaussianERes* gaus,const size_t& nDims);
    GaussianEResFitter(GaussianERes* gaus, const std::vector<std::string>&);

    void   SetParameter(const std::string& name_, double value);
    double GetParameter(const std::string& name_) const;

    void   SetParameters(const ParameterDict&);
    ParameterDict GetParameters() const;
    size_t GetParameterCount() const;

    std::set<std::string> GetParameterNames() const;
    void   RenameParameter(const std::string& old_, const std::string& new_);

    std::vector<std::string> GetEResNames() const;

private:
    GaussianERes* fOrignalFunc;
    std::vector<std::string> fEResNames;
};
#endif
