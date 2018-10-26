#ifndef __OXSX_VARYING_PARAMETER_PDF__
#define __OXSX_VARYING_PARAMETER_PDF__
#include <string>
#include <map>
#include <vector>

#include <PDF.h>
#include <ConditionalPDF.h>
#include <ParameterDict.h>
#include <Function.h>
#include <algorithm>

class VaryingCDF : public ConditionalPDF{
 public:
    VaryingCDF(const std::string&);
    ~VaryingCDF();

    VaryingCDF(const VaryingCDF& copy_);

    VaryingCDF& operator=(const VaryingCDF& other_);

    double ConditionalProbability(const std::vector<double>& x_, 
            const std::vector<double>& x2_);
    double Integral(const std::vector<double>& mins_, 
                    const std::vector<double>& maxs_,
                    const std::vector<double>& x2_) const;

    std::vector<double> Sample(const std::vector<double>&) const;

    ConditionalPDF* Clone() const;

    void SetKernel(ConditionalPDF* PDF_);
    void SetKernel(PDF* PDF_);
    void SetDependance(const std::string& paraName_, const Function* func_);

    ParameterDict SetupParameters(const std::vector<double>& x2_) const;

    // FitParameter interface.
    void   SetParameter(const std::string& name_, double value);
    double GetParameter(const std::string& name_) const;

    void   SetParameters(const ParameterDict&);
    ParameterDict GetParameters() const;
    size_t GetParameterCount() const;

    std::set<std::string> GetParameterNames() const;
    void   RenameParameter(const std::string& old_, const std::string& new_);

    std::string GetName() const;
    void SetName(const std::string&);

 private:
    std::map<std::string, Function*> fFunctions;
    ConditionalPDF* fCdf;
    
    std::string fName;
};
#endif
