#include <Gaussian.h>
#include <GaussianFitter.h>
#include <Exceptions.h>
#include <ContainerParameter.h>
#include <ContainerTools.hpp>
#include <Formatter.hpp>
#include <gsl/gsl_cdf.h>
#include <Rand.h>
#include <sstream>
#include <math.h>

/////////////////////////
// Constructory Things //
/////////////////////////

// Constructory things
Gaussian::Gaussian() : fFitter(this,1){
    Initialise(std::vector<double>(1, 0), std::vector<double>(1, 1), "");
}

Gaussian::Gaussian(size_t nDims_, const std::string& name_): fFitter(this,nDims_){
    Initialise(std::vector<double>(nDims_, 0), std::vector<double>(nDims_, 1), name_);
}// means = 0, stdDevs = 1

Gaussian::Gaussian(double mean_, double stdDev_, const std::string& name_ ): fFitter(this,1){
    Initialise(std::vector<double>(1, mean_), std::vector<double>(1, stdDev_), name_);
}

Gaussian::Gaussian(const std::vector<double>& mean_, 
                   const std::vector<double>& stdDev_,
                   const std::string& name_ ): fFitter(this,mean_.size()){
    Initialise(mean_, stdDev_, name_);
}

Gaussian::Gaussian(const Gaussian& copy_): fFitter(this,copy_.GetMeanNames(),copy_.GetStDevNames()){
    fMeans = copy_.fMeans;
    fStdDevs = copy_.fStdDevs;
    fCdfCutOff = copy_.fCdfCutOff;
    fNDims = copy_.fNDims;
    fName = std::string(copy_.fName+"_copy");
}

Gaussian&
Gaussian::operator=(const Gaussian& copy_){
    fMeans = copy_.fMeans;
    fStdDevs = copy_.fStdDevs;
    fCdfCutOff = copy_.fCdfCutOff;
    fNDims = copy_.fNDims;
    fName = std::string(copy_.fName+"_copy");
    fFitter = GaussianFitter(this,copy_.GetMeanNames(),copy_.GetStDevNames());
    return *this;
}

void
Gaussian::Initialise(const std::vector<double>& means_, const std::vector<double>& stdDevs_, 
                     const std::string& name_){
    if (name_ == "")
        fName = "gaussian";
    else
        fName = name_;
    SetMeansStdDevs(means_, stdDevs_);
    fNDims   = means_.size() ;
    fMeans   = means_;
    fStdDevs = stdDevs_;
    fCdfCutOff = 6; // default val
}

Function* 
Gaussian::Clone() const{
    return static_cast<Function*> (new Gaussian(*this));
}

/////////////////////
// Getters/Setters //
/////////////////////

double 
Gaussian::GetMean(size_t dimension_) const{
    try{
        return fMeans.at(dimension_);
    }
    catch(const std::out_of_range& e_){
        throw NotFoundError("Gaussian::Requested Gaussian mean beyond function dimensionality!");
    }
}

double 
Gaussian::GetStDev(size_t dimension_) const{
    try{
        return fStdDevs.at(dimension_);
    }
    catch(const std::out_of_range& e_){
        throw NotFoundError("Gaussian::Requested Gaussian stdDev beyond function dimensionality!");
    }
}

void
Gaussian::SetMean(const size_t& dim_ , const double& value_) {
    fMeans[dim_]= value_;
}

void
Gaussian::SetStDev(const size_t& dim_ , const double& value_) {
    fStdDevs[dim_]= value_;
}

void
Gaussian::SetMeans(const std::vector<double>& means_) {
    fMeans = means_;
}

void
Gaussian::SetStDevs(const std::vector<double>& stddevs_) {
    fStdDevs = stddevs_;
}

std::vector<double>
Gaussian::GetMeans() const {
    return fMeans;
}

void
Gaussian::SetMeansStdDevs(const std::vector<double>& means_, 
                          const std::vector<double>& stdDevs_){
    if (means_.size() != stdDevs_.size())
        throw DimensionError("Gaussian::Tried to set Gaussian function with #means != #stdDevs!");

    for(size_t i = 0; i < stdDevs_.size(); i++)
        if(stdDevs_.at(i) <= 0)
            throw ValueError("Gaussian::Gaussian standard deviation must be greater than 0!");

    fMeans = means_;
    fStdDevs = stdDevs_;
    fNDims = means_.size();
}

std::vector<double>
Gaussian::GetStdDevs() const {
    return fStdDevs;
}

double
Gaussian::GetCdfCutOff() const{
    return fCdfCutOff;
}

void 
Gaussian::SetCdfCutOff(double cutOff_){
    fCdfCutOff = cutOff_;
}

int 
Gaussian::GetNDims() const{
    return fNDims;
}

/////////////////
// Probability //
/////////////////

double 
Gaussian::operator() (const std::vector<double>& vals_) const{
    if (vals_.size() != GetNDims())
        throw DimensionError("Gaussian::Gaussian dimensionality does not match the input vector to evaluate!");

    double exponent = 0;
    double mean;
    double stdDev;
    double nDevs;
    for(size_t i = 0; i < GetNDims(); i++){
        mean  = fMeans.at(i);
        stdDev = fStdDevs.at(i);

        nDevs = (vals_.at(i) - mean)/stdDev;
        exponent += nDevs * nDevs;
    }

    double norm = 1;
    for(size_t i = 0; i < GetNDims(); i++)
        norm *= sqrt(2*M_PI) * fStdDevs.at(i);

    return exp(- 0.5 * exponent) / norm; 
}

double 
Gaussian::Cdf(size_t dim_, double val_) const{
    double nDevs = (val_ - GetMean(dim_))/GetStDev(dim_);
    if (nDevs > fCdfCutOff)
        return 1;
    if(nDevs < -1 * fCdfCutOff)
        return 0;

    return gsl_cdf_gaussian_P(val_ - GetMean(dim_), GetStDev(dim_));
}

double 
Gaussian::Integral(const std::vector<double>& mins_, const std::vector<double>& maxs_) const{
    if(mins_.size() != GetNDims() || maxs_.size() != GetNDims())
        throw DimensionError("Gaussian::Gaussian, tried to integrate over interval of wrong dimensionality");

    double integral = 1;
    for(size_t i = 0; i < mins_.size(); i++)
        integral *= ( Cdf(i, maxs_[i]) - Cdf(i, mins_[i]));
  
    return integral;  
}

double
Gaussian::Integral(const std::vector<double>& mins_, const std::vector<double>& maxs_, const double& bincentre) const{
    std::cout<<"Are you using energy depency for SurvProb? This function is currently unusable"<<std::endl;
    exit(1);
    if(mins_.size() != GetNDims() || maxs_.size() != GetNDims())
        throw DimensionError("Gaussian::Gaussian, tried to integrate over interval of wrong dimensionality");

    double integral = 1;
    for(size_t i = 0; i < mins_.size(); i++)
        integral *= ( Cdf(i, maxs_[i]) - Cdf(i, mins_[i]));
  
    return integral;  
}


std::vector<double>
Gaussian::Sample() const{
  std::vector<double> sample(GetNDims(), 0);
  for(size_t i = 0; i < GetNDims(); i++){
    sample[i] = Rand::Gaus(fMeans.at(i), fStdDevs.at(i));
  }
  return sample;
}

////////////////////////
// Make it fittable   //
////////////////////////
void
Gaussian::RenameParameter(const std::string& old_, const std::string& new_){
    fFitter.RenameParameter(old_, new_);
}

void
Gaussian::SetParameter(const std::string& name_, double value_){
    fFitter.SetParameter(name_, value_);
}

double
Gaussian::GetParameter(const std::string& name_) const{
    return fFitter.GetParameter(name_);
}

void
Gaussian::SetParameters(const ParameterDict& ps_){
     fFitter.SetParameters(ps_);
}

ParameterDict
Gaussian::GetParameters() const{
    return fFitter.GetParameters();
}

size_t
Gaussian::GetParameterCount() const{
    return fFitter.GetParameterCount();
}

std::set<std::string>
Gaussian::GetParameterNames() const{
    return fFitter.GetParameterNames();
}

std::vector<std::string>
Gaussian::GetMeanNames() const{
    return fFitter.GetMeanNames();
}

std::vector<std::string>
Gaussian::GetStDevNames() const{
    return fFitter.GetStdDevNames();
}

std::string
Gaussian::GetName() const{
    return fName;
}

void
Gaussian::SetName(const std::string& name_){
    fName = name_;
}
