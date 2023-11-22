/*****************************************/
/* A simple scale error on an observable */
/*****************************************/
#ifndef __OXSX_SCALE__
#define __OXSX_SCALE__
#include <Systematic.h>
#include <string>
#include <DenseMatrix.h>

class Scale : public Systematic{
 public:
    Scale(const std::string& name_) : fScaleFactor(1), fName(name_), fParamName("scaleFactor"), fCachedBinMapping(false) {}
    void   SetScaleFactor(double);
    double GetScaleFactor() const;
    
    void Construct();

    // Adjustable scale factor
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
    double      fScaleFactor;
    std::string fName;
    std::string fParamName;

    std::vector<size_t> fDistTransBinMapping;
    DenseMatrix fMappingDistAndTrans;
    bool fCachedBinMapping;

    void CacheBinMapping();
    double GetBinContribution(const BinAxis& scaleAxis, double scaledLow,
                                     double scaledHigh, size_t bin_test);
};
#endif
