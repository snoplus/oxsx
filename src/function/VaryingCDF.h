/************************************************************************
 * VaryingCDF is a CDF who's parameters are governed by a general function of
 * inputs.
 *
 * Use case : A spectrum's resolution may have some functional dependence on
 * the spectrum itself. Convolving the spectrum with a gaussian of constant
 * width doesn't encode this functional dependence. The VaryingCDF can take the
 * gaussian (or any function) as it's kernel and change the parameters of this
 * kernel by some general function dependent on the current spectrum point. 
 *
 * Usages:
 *    Gaussian gaus("gaus",0,1);
 *    gaus.RenameParameter("means_0","mean");
 *    gaus.RenameParameter("stddevs_0","std");
 *    VaryingCDF smearer("smearer");
 *    smearer.SetKernal(&gaus)
 *    smearer.SetDependance
 *    smearer.SetDependance("std",ploy(see note));
 *
 * The above creates a VaryingCDF whose gaussian kernel has a variable
 * stddev which goes like the function ploy.

 * Note:
 * Where 'ploy' is a function that takes a vector of doubles that
 * represents the bin in which you are evaluating on e.g. fitting in r,E
 * ploy will take a vector = [r_i,E_j]. The order of the vector is the same
 * as the order in which the axes of the fit were set.
 * 
 ************************************************************************/


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
