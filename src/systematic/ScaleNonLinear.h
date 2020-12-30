/*****************************************/
/* A simple scale error on an observable */
/*****************************************/
#ifndef __OXSX_SCALENONLINEAR__
#define __OXSX_SCALENONLINEAR__
#include <Systematic.h>
#include <string>

class ScaleNonLinear : public Systematic{
 public:
   ScaleNonLinear(const std::string& name_) : fScaleFactor1(0.), fScaleFactor2(1), fName(name_), fParamName1("scaleFactor1"), fParamName2("scaleFactor2") {}
    double NonLinearScaleFunction(double);

    void   SetScaleFactor1(const double);
    double GetScaleFactor1() const;

    void   SetScaleFactor2(const double);
    double GetScaleFactor2() const;
    
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
    double    fScaleFactor1;
    double    fScaleFactor2;
    
    std::string fName;
    std::string fParamName1;
    std::string fParamName2;
};
#endif
