#ifndef __OXSX_SURVPROB__
#define __OXSX_SURVPROB__
#include <PDF.h>
#include <ParameterManager.h>
#include <string>
#include <SurvProbFitter.h>
#include <TF1.h>

class SurvProb : public PDF{
 public:
    // Constructory things
    SurvProb();
    SurvProb(size_t nDims_, const std::string& name_ = "");
    SurvProb(double delmsqr21_, double sinsqrtheta12_, double baseline_, const std::string& name_ = ""); //double sinsqrtheta13_,
    SurvProb(const std::vector<double>& delmsqr21_, const std::vector<double>& sinsqrtheta12_, const std::string& name_ = ""); //const std::vector<double>& sinsqrtheta13_,

    SurvProb(const SurvProb& copy_);

    SurvProb& operator=(const SurvProb& other_);

    virtual   Function* Clone() const;

    // Probability
    virtual double operator()(const std::vector<double>& vals_) const;
    double  Cdf(size_t dim_, double val_) const;
    double Integral(const std::vector<double>& mins_,
                    const std::vector<double>& maxs_) const;
    double Integral() const {return 1;} // normalised by definition
    std::vector<double> Sample() const;

    // Getters/Setters
    double Getdelmsqr21(size_t dimension_) const;
    double Getsinsqrtheta12(size_t dimension_) const;
    //double Getsinsqrtheta13(size_t dimension_) const;

    void Setdelmsqr21(const size_t& dim_ , const double& value_);
    void Setsinsqrtheta12(const size_t& dim_ , const double& value_);
    void Setsinsqrtheta13(const size_t& dim_ , const double& value_);

    std::vector<double> Getdelmsqr21s() const;
    std::vector<double> Getsinsqrtheta12s() const;
    //std::vector<double> Getsinsqrtheta13s() const;
    double Getsinsqrtheta13s();
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
    void Setdelmsqr21s(const std::vector<double>& delmsqr21s_);
    void Setsinsqrtheta12s(const std::vector<double>& sinsqrtheta12s_);
    //void Setsinsqrtheta13s(const std::vector<double>& sinsqrtheta13s_);
    void Setsinsqrtheta13s(double sinsqrtheta13s_);
    std::vector<std::string> Getdelmsqr21Names() const;
    std::vector<std::string> Getsinsqrtheta12Names() const;
    //std::vector<std::string> Getsinsqrtheta13Names() const;

 private:
    SurvProbFitter fFitter;
    std::vector<double> fdelmsqr21s;
    std::vector<double> fsinsqrtheta12s;
    //std::vector<double> fsinsqrtheta13s;
    double fsinsqrtheta13s;
    double fBaseline;

    double fCdfCutOff; // number of stDevs away from the mean
                       // assumed to be zero or 1 for speed integration
    int fNDims;
    std::string fName;

    TF1 *fsincfunction;
    double Si(double val_) const;

    void   Initialise(const std::vector<double>& delmsqr21s_,
                      const std::vector<double>& sinsqrtheta12s_,
                      const std::string& name_); //const std::vector<double>& sinsqrtheta13s_,

    // this is private, we want the dimensionality to be fixed at creation
    void Setdelmsqrssinsqrtheta12ssinsqrtheta13s(const std::vector<double>& delmsqr21s_,
						 const std::vector<double>& sinsqrtheta12s_); //, const std::vector<double>& sinsqrtheta13s_
};
#endif
