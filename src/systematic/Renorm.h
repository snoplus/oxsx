/****************************************************************/
/* A systematic that applies a renormalisation                  */
/****************************************************************/
#ifndef __OXSX_RENORM__
#define __OXSX_RENORM__
#include <Systematic.h>

class Renorm : public Systematic
{
public:
Renorm(const std::string &name_) : fParamDict(), fName(name_){}
    void Construct();

    void SetScaleFactor(double);
    double GetScaleFactor() const;

    // FitComponent interface: handling systematic's parameters
    void SetParameter(const std::string &name_, double value);
    double GetParameter(const std::string &name_) const;

    void SetParameters(const ParameterDict &);
    ParameterDict GetParameters() const { return fParamDict; }
    size_t GetParameterCount() const { return fParamDict.size(); }

    std::set<std::string> GetParameterNames() const;
    void RenameParameter(const std::string &old_, const std::string &new_);

    std::string GetName() const { return fName; }
    void SetName(const std::string &name_) { fName = name_; }

private:
    double fScaleFactor;
    ParameterDict fParamDict;
    std::string fName;
};

#endif
