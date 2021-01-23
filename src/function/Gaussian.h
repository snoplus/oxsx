#ifndef __OXSX_GAUSSIAN__
#define __OXSX_GAUSSIAN__
#include <PDF.h>
#include <ParameterManager.h>
#include <string>
#include <GaussianFitter.h>

class Gaussian : public PDF{
 public:
    // Constructory things
    Gaussian();
    Gaussian(size_t nDims_, const std::string& name_ = "");// means = 0, stdDevs = 1
    Gaussian(double mean_, double stdDev_, const std::string& name_ = "");
    Gaussian(const std::vector<double>& mean_, 
            const std::vector<double>& stdDev_, const std::string& name_ = "");

    Gaussian(const Gaussian& copy_);

    Gaussian& operator=(const Gaussian& other_);

    virtual   Function* Clone() const;

    // Probability
    virtual double operator()(const std::vector<double>& vals_) const;
    double  Cdf(size_t dim_, double val_) const;
    double Integral(const std::vector<double>& mins_, 
                    const std::vector<double>& maxs_) const;
    double Integral() const {return 1;} // normalised by definition
    std::vector<double> Sample() const;

    // Getters/Setters
    double GetMean(size_t dimension_) const;
    double GetStDev(size_t dimension_) const;    

    void SetMean(const size_t& dim_ , const double& value_);
    void SetStDev(const size_t& dim_ , const double& value_);

    std::vector<double> GetMeans() const;
    std::vector<double> GetStdDevs() const;
    double GetCdfCutOff() const;
    void   SetCdfCutOff(double);
    size_t GetNDims() const;
    
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
    void SetMeans(const std::vector<double>& means_);
    void SetStDevs(const std::vector<double>& stddev_);
    std::vector<std::string> GetMeanNames() const;
    std::vector<std::string> GetStDevNames() const;
 private:
    GaussianFitter fFitter;
    std::vector<double> fMeans;
    std::vector<double> fStdDevs;
    
    double fCdfCutOff; // number of stDevs away from the mean
                       // assumed to be zero or 1 for speed integration
    size_t fNDims;
    std::string fName;
    
    void   Initialise(const std::vector<double>& means_, 
                      const std::vector<double>& stdDevs_,
                      const std::string& name_);
    
    // this is private, we want the dimensionality to be fixed at creation
    void SetMeansStdDevs(const std::vector<double>& means_, 
                         const std::vector<double>& stdDevs_);
};
#endif
