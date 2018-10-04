/*****************************************************************************************************/
//
/*****************************************************************************************************/

#ifndef __OXSX_NUOSC__
#define __OXSX_NUOSC__
#include <Systematic.h>

class PDF;
class NuOsc : public Systematic{
 public:
    NuOsc(const std::string& name_): fDist(NULL), fCachedCompatibleBins(false), fName(name_) {}
    ~NuOsc();
    void SetFunction(PDF* function_);
    void Construct();    

    // Make this fittable, by delegating to the underlying pdf
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
    void                     Reset();
    PDF*                     fDist;
    bool                     fCachedCompatibleBins;

    AxisCollection  fSubMapAxes;
    void  CacheCompatibleBins();
    std::vector<std::vector<size_t> > fCompatibleBins;
    // the systematic subMap bin for each global bin of pdf
    std::vector<size_t> fSysBins; 

    std::string fName;
    
};
#endif
