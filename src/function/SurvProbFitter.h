#ifndef __OXSX_SURVPROB_FITTER__
#define __OXSX_SURVPROB_FITTER__
#include <vector>      
#include <string>
#include <ParameterDict.h>
#include <set>
#include <iostream>

class SurvProb;

class SurvProbFitter{
public:

    SurvProbFitter(SurvProb* gaus,const size_t& nDims);
    // SurvProbFitter(SurvProb* gaus, const std::vector<std::string>&, const std::vector<std::string>&);
    SurvProbFitter(SurvProb* gaus, const std::vector<std::string>&, const std::vector<std::string>&); //,const std::vector<std::string>&

    void   SetParameter(const std::string& name_, double value);
    double GetParameter(const std::string& name_) const;

    void   SetParameters(const ParameterDict&);
    ParameterDict GetParameters() const;
    size_t GetParameterCount() const;

    std::set<std::string> GetParameterNames() const;
    void   RenameParameter(const std::string& old_, const std::string& new_);

    std::vector<std::string> GetDelmsqr21Names() const;
    std::vector<std::string> GetSinsqrtheta12Names() const;
    //std::vector<std::string> GetSinsqrtheta13Names() const;

private:
    SurvProb* fOrignalFunc;
    std::vector<std::string> fDelmsqr21sNames;
    std::vector<std::string> fSinsqrtheta12sNames;
    //std::vector<std::string> fSinsqrtheta13sNames;
};
#endif
