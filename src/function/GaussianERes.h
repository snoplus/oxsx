#ifndef __OXSX_GAUSSIANERES__
#define __OXSX_GAUSSIANERES__
#include <PDF.h>
#include <ParameterManager.h>
#include <string>
#include <GaussianEResFitter.h>

class GaussianERes : public PDF{
 public:
    // Constructory things
    GaussianERes();
    GaussianERes(size_t nDims_, const std::string& name_ = "");// means = 0, stdDevs = 1
    GaussianERes(double ERes_, const std::string& name_ = "");
    GaussianERes(const std::vector<double>& ERes_, const std::string& name_ = "");

    GaussianERes(const GaussianERes& copy_);

    GaussianERes& operator=(const GaussianERes& other_);

    virtual   Function* Clone() const;

    // Probability
    virtual double operator()(const std::vector<double>& vals_) const;
    double  Cdf(size_t dim_, double val_) const;
    double  Cdf(size_t dim_, double val_, const double bincentre) const;
    double Integral(const std::vector<double>& mins_, 
                    const std::vector<double>& maxs_) const;
    double Integral(const std::vector<double>& mins_, 
                    const std::vector<double>& maxs_,
                    const double& bincentre) const;
    double Integral() const {return 1;} // normalised by definition
    std::vector<double> Sample() const;

    // Getters/Setters
    double GetERes(size_t dimension_) const;    

    void SetERes(const size_t& dim_ , const double& value_);

    std::vector<double> GetERess() const;
    double GetCdfCutOff() const;
    void   SetCdfCutOff(double);
    int    GetNDims() const;
    
    // Make this object fittable
    void   SetParameter(const std::string& name_, double value);
    double GetParameter(const std::string& name_) const;
    
    void   SetParameters(const ParameterDict&);
    ParameterDict GetParameters() const;
    size_t GetParameterCount() const;
    
    std::set<std::string> GetParameterNames() const;
    void   RenameParameter(const std::string& old_, const std::string& new_);
    
    std::string GetName() const;
    void SetName(const std::string&);
    void SetERess(const std::vector<double>& ERes_);
    std::vector<std::string> GetEResNames() const;
 private:
    GaussianEResFitter fFitter;
    std::vector<double> fERess;
    
    double fCdfCutOff; // number of stDevs away from the mean
                       // assumed to be zero or 1 for speed integration
    int fNDims;
    std::string fName;
    
    void   Initialise(const std::vector<double>& ERess_,
                      const std::string& name_);
    
};
#endif
