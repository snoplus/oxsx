#include <GaussianERes.h>
#include <GaussianEResFitter.h>
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
GaussianERes::GaussianERes() : fFitter(this,1){
    Initialise(std::vector<double>(1, 1), "");
}

GaussianERes::GaussianERes(size_t nDims_, const std::string& name_): fFitter(this,nDims_){
    Initialise(std::vector<double>(nDims_, 1), name_);
}// means = 0, stdDevs = 1

GaussianERes::GaussianERes(double ERes_, const std::string& name_ ): fFitter(this,1){
    Initialise(std::vector<double>(1, ERes_), name_);
}

GaussianERes::GaussianERes(const std::vector<double>& ERes_,
                   const std::string& name_ ): fFitter(this,ERes_.size()){
    Initialise(ERes_, name_);
}

GaussianERes::GaussianERes(const GaussianERes& copy_): fFitter(this,copy_.GetEResNames()){
    fERess = copy_.fERess;
    fCdfCutOff = copy_.fCdfCutOff;
    fNDims = copy_.fNDims;
    fName = std::string(copy_.fName+"_copy");
}

GaussianERes&
GaussianERes::operator=(const GaussianERes& copy_){
    fERess = copy_.fERess;
    fCdfCutOff = copy_.fCdfCutOff;
    fNDims = copy_.fNDims;
    fName = std::string(copy_.fName+"_copy");
    fFitter = GaussianEResFitter(this,copy_.GetEResNames());
    return *this;
}

void
GaussianERes::Initialise(const std::vector<double>& ERess_, 
                     const std::string& name_){
    if (name_ == "")
        fName = "gaussianeres";
    else
        fName = name_;
    SetERess(ERess_);
    fNDims   = ERess_.size() ;
    fERess = ERess_;
    fCdfCutOff = 6; // default val
}

Function* 
GaussianERes::Clone() const{
    return static_cast<Function*> (new GaussianERes(*this));
}

/////////////////////
// Getters/Setters //
/////////////////////

double 
GaussianERes::GetERes(size_t dimension_) const{
    try{
        return fERess.at(dimension_);
    }
    catch(const std::out_of_range& e_){
        throw NotFoundError("GaussianERes::Requested GaussianERes ERes beyond function dimensionality!");
    }
}

void
GaussianERes::SetERes(const size_t& dim_ , const double& value_) {
    fERess[dim_]= value_;
}

void
GaussianERes::SetERess(const std::vector<double>& ERess_) {
    fERess = ERess_;
}

std::vector<double>
GaussianERes::GetERess() const {
    return fERess;
}

double
GaussianERes::GetCdfCutOff() const{
    return fCdfCutOff;
}

void 
GaussianERes::SetCdfCutOff(double cutOff_){
    fCdfCutOff = cutOff_;
}

int 
GaussianERes::GetNDims() const{
    return fNDims;
}

/////////////////
// Probability //
/////////////////

double 
GaussianERes::operator() (const std::vector<double>& vals_) const{
    if (vals_.size() != GetNDims())
        throw DimensionError("GaussianERes::GaussianERes dimensionality does not match the input vector to evaluate!");

    double exponent = 0;
    double stdDev;
    double nDevs;
    for(size_t i = 0; i < GetNDims(); i++){
        stdDev = sqrt(abs(vals_.at(i)))*sqrt(pow(1+fERess.at(i),2) - 1); //fERess.at(i)*sqrt(vals_.at(i));

        nDevs = (vals_.at(i))/stdDev;
        exponent += nDevs * nDevs;
    }

    double norm = 1;
    for(size_t i = 0; i < GetNDims(); i++)
        norm *= sqrt(2*M_PI) * (sqrt(abs(vals_.at(i)))*sqrt(pow(1+fERess.at(i),2) - 1)); //fERess.at(i)*sqrt(vals_.at(i));
    
    return exp(- 0.5 * exponent) / norm; 
}

double 
GaussianERes::Cdf(size_t dim_, double val_) const{
    double nDevs = (val_)/GetERes(dim_);
    if (nDevs > fCdfCutOff)
        return 1;
    if(nDevs < -1 * fCdfCutOff)
        return 0;

    return gsl_cdf_gaussian_P(val_, GetERes(dim_));
}

double 
GaussianERes::Cdf(size_t dim_, double val_, const double bincentre) const{
    double stdDev = sqrt(bincentre)*sqrt(pow(1+GetERes(dim_),2) - 1);
    double nDevs = (val_ - bincentre)/stdDev;
    if (nDevs > fCdfCutOff)
        return 1;
    if(nDevs < -1 * fCdfCutOff)
        return 0;

    //std::cout<<"val "<<val_<<" centre "<<bincentre<<" eres "<<GetERes(dim_)<<" stddev "<<stdDev<<" ndevs "<<nDevs<<" = "<<gsl_cdf_gaussian_P(val_ - bincentre, stdDev)<<std::endl;
    return gsl_cdf_gaussian_P(val_ - bincentre, stdDev);
}

double 
GaussianERes::Integral(const std::vector<double>& mins_, const std::vector<double>& maxs_) const{
    if(mins_.size() != GetNDims() || maxs_.size() != GetNDims())
        throw DimensionError("GaussianERes::GaussianERes, tried to integrate over interval of wrong dimensionality");

    double integral = 1;
    for(size_t i = 0; i < mins_.size(); i++)
        integral *= ( Cdf(i, maxs_[i]) - Cdf(i, mins_[i]));

    return integral;  
}

double 
GaussianERes::Integral(const std::vector<double>& mins_, const std::vector<double>& maxs_, const double& bincentre) const{
    if(mins_.size() != GetNDims() || maxs_.size() != GetNDims())
        throw DimensionError("GaussianERes::GaussianERes, tried to integrate over interval of wrong dimensionality");

    double integral = 1;
    for(size_t i = 0; i < mins_.size(); i++)
      integral *= ( Cdf(i, maxs_[i], bincentre) - Cdf(i, mins_[i], bincentre));
    
    return integral;  
}


std::vector<double>
GaussianERes::Sample() const{
  std::cout<<"this sampling function is not been finished for EResolution"<<std::endl;
  std::vector<double> sample(GetNDims(), 0);
  for(size_t i = 0; i < GetNDims(); i++){
    sample[i] = Rand::Gaus(0, fERess.at(i));
  }
  return sample;
}

////////////////////////
// Make it fittable   //
////////////////////////
void
GaussianERes::RenameParameter(const std::string& old_, const std::string& new_){
    fFitter.RenameParameter(old_, new_);
}

void
GaussianERes::SetParameter(const std::string& name_, double value_){
    fFitter.SetParameter(name_, value_);
}

double
GaussianERes::GetParameter(const std::string& name_) const{
    return fFitter.GetParameter(name_);
}

void
GaussianERes::SetParameters(const ParameterDict& ps_){
     fFitter.SetParameters(ps_);
}

ParameterDict
GaussianERes::GetParameters() const{
    return fFitter.GetParameters();
}

size_t
GaussianERes::GetParameterCount() const{
    return fFitter.GetParameterCount();
}

std::set<std::string>
GaussianERes::GetParameterNames() const{
    return fFitter.GetParameterNames();
}

std::vector<std::string>
GaussianERes::GetEResNames() const{
    return fFitter.GetEResNames();
}

std::string
GaussianERes::GetName() const{
    return fName;
}

void
GaussianERes::SetName(const std::string& name_){
    fName = name_;
}
