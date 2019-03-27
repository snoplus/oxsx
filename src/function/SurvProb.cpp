#include <SurvProb.h>
#include <SurvProbFitter.h>
#include <Exceptions.h>
#include <ContainerParameter.h>
#include <ContainerTools.hpp>
#include <Formatter.hpp>
#include <gsl/gsl_cdf.h>
#include <Rand.h>
#include <sstream>
#include <math.h>
#include <TF1.h>

/////////////////////////
// Constructory Things //
/////////////////////////

// Constructory things
SurvProb::SurvProb() : fFitter(this,1),fBaseline(0){
  Initialise(std::vector<double>(1, 0),std::vector<double>(1, 0), ""); //, std::vector<double>(1, 0)
}

SurvProb::SurvProb(size_t nDims_, const std::string& name_): fFitter(this,nDims_),fBaseline(0){
  Initialise(std::vector<double>(nDims_, 0),std::vector<double>(nDims_, 0), name_); //, std::vector<double>(nDims_, 0)
}// means = 0, stdDevs = 1

SurvProb::SurvProb(double delmsqr21_, double sinsqrtheta12_, double baseline_, const std::string& name_ ): fFitter(this,1), fBaseline(baseline_){
  Initialise(std::vector<double>(1, delmsqr21_), std::vector<double>(1, sinsqrtheta12_), name_); // , std::vector<double>(1, sinsqrtheta13_) double sinsqrtheta13_,
}

SurvProb::SurvProb(const std::vector<double>& delmsqr21_, 
                   const std::vector<double>& sinsqrtheta12_, 
                   const std::string& name_ ): fFitter(this,delmsqr21_.size()),fBaseline(0){ //		   const std::vector<double>& sinsqrtheta13_,
  Initialise(delmsqr21_, sinsqrtheta12_, name_); //sinsqrtheta13_,
}

SurvProb::SurvProb(const SurvProb& copy_): fFitter(this,copy_.Getdelmsqr21Names(),copy_.Getsinsqrtheta12Names()){ //,copy_.Getsinsqrtheta13Names()
    fdelmsqr21s = copy_.fdelmsqr21s;
    fsinsqrtheta12s = copy_.fsinsqrtheta12s; //    fsinsqrtheta13s = copy_.fsinsqrtheta13s;
    fCdfCutOff = copy_.fCdfCutOff;
    fNDims = copy_.fNDims;
    fBaseline = copy_.fBaseline;
    fName = std::string(copy_.fName+"_copy");
}

SurvProb&
SurvProb::operator=(const SurvProb& copy_){
    fdelmsqr21s = copy_.fdelmsqr21s;
    fsinsqrtheta12s = copy_.fsinsqrtheta12s; //    fsinsqrtheta13s = copy_.fsinsqrtheta13s;
    fCdfCutOff = copy_.fCdfCutOff;
    fNDims = copy_.fNDims;
    fBaseline = copy_.fBaseline;
    fName = std::string(copy_.fName+"_copy");
    fFitter = SurvProbFitter(this,copy_.Getdelmsqr21Names(),copy_.Getsinsqrtheta12Names()); //,copy_.Getsinsqrtheta13Names()
    return *this;
}

void
SurvProb::Initialise(const std::vector<double>& delmsqr21s_, const std::vector<double>& sinsqrtheta12s_,const std::string& name_){ //const std::vector<double>& sinsqrtheta13s_
    if (name_ == "")
        fName = "survprob";
    else
        fName = name_;
    Setdelmsqrssinsqrtheta12ssinsqrtheta13s(delmsqr21s_, sinsqrtheta12s_); //, sinsqrtheta13s_
    fNDims   = delmsqr21s_.size() ;
    fdelmsqr21s  = delmsqr21s_;
    fsinsqrtheta12s = sinsqrtheta12s_;
    //fsinsqrtheta13s = sinsqrtheta13s_;
    fCdfCutOff = 6; // default val
    fsincfunction = new TF1("sinc_function","sin(x)/x",0,1000);
}

Function* 
SurvProb::Clone() const{
    return static_cast<Function*> (new SurvProb(*this));
}

/////////////////////
// Getters/Setters //
/////////////////////

double 
SurvProb::Getdelmsqr21(size_t dimension_) const{
    try{
        return fdelmsqr21s.at(dimension_);
    }
    catch(const std::out_of_range& e_){
        throw NotFoundError("SurvProb::Requested SurvProb mean beyond function dimensionality!");
    }
}

double 
SurvProb::Getsinsqrtheta12(size_t dimension_) const{
    try{
        return fsinsqrtheta12s.at(dimension_);
    }
    catch(const std::out_of_range& e_){
        throw NotFoundError("SurvProb::Requested SurvProb stdDev beyond function dimensionality!");
    }
}

//double 
//SurvProb::Getsinsqrtheta13(size_t dimension_) const{
    //try{
    //    return fsinsqrtheta13s.at(dimension_);
    //}
    //catch(const std::out_of_range& e_){
    //    throw NotFoundError("SurvProb::Requested SurvProb stdDev beyond function dimensionality!");
    //}
//}
double 
SurvProb::Getsinsqrtheta13s() {
    return fsinsqrtheta13s;
}

void
SurvProb::Setdelmsqr21(const size_t& dim_ , const double& value_) {
    fdelmsqr21s[dim_]= value_;
}

void
SurvProb::Setsinsqrtheta12(const size_t& dim_ , const double& value_) {
    fsinsqrtheta12s[dim_]= value_;
}

//void
//SurvProb::Setsinsqrtheta13(const size_t& dim_ , const double& value_) {
//    fsinsqrtheta13s[dim_]= value_;
//}

void
SurvProb::Setdelmsqr21s(const std::vector<double>& delmsqr21s_) {
    fdelmsqr21s = delmsqr21s_;
}

void
SurvProb::Setsinsqrtheta12s(const std::vector<double>& sinsqrtheta12s_) {
    fsinsqrtheta12s = sinsqrtheta12s_;
}

//void
//SurvProb::Setsinsqrtheta13s(const std::vector<double>& sinsqrtheta13s_) {
//    fsinsqrtheta13s = sinsqrtheta13s_;
//}
void
SurvProb::Setsinsqrtheta13s(const double sinsqrtheta13s_) {
    fsinsqrtheta13s = sinsqrtheta13s_;
}

std::vector<double>
SurvProb::Getdelmsqr21s() const {
    return fdelmsqr21s;
}

void
SurvProb::Setdelmsqrssinsqrtheta12ssinsqrtheta13s(const std::vector<double>& delmsqr21s_, 
						  const std::vector<double>& sinsqrtheta12s_){ //const std::vector<double>& sinsqrtheta13s_
    if (delmsqr21s_.size() != sinsqrtheta12s_.size()) //!= sinsqrtheta13s_.size()
        throw DimensionError("SurvProb::Tried to set SurvProb function with #delmqrs != #sinsqrtheta12s or #sinsqrtheta13s !");

    for(size_t i = 0; i < delmsqr21s_.size(); i++)
        if(delmsqr21s_.at(i) <= 0)
            throw ValueError("SurvProb::SurvProb delmsqr21s must be greater than 0!");
    for(size_t i = 0; i < sinsqrtheta12s_.size(); i++)
        if(sinsqrtheta12s_.at(i) <= 0)
            throw ValueError("SurvProb::SurvProb sinsqrtheta12 must be greater than 0!");

    //for(size_t i = 0; i < sinsqrtheta13s_.size(); i++)
    //    if(sinsqrtheta13s_.at(i) <= 0)
    //        throw ValueError("SurvProb::SurvProb sinsqrtheta13 must be greater than 0!");

    fdelmsqr21s = delmsqr21s_;
    fsinsqrtheta12s = sinsqrtheta12s_;
    //fsinsqrtheta13s = sinsqrtheta13s_;
    fNDims = delmsqr21s_.size();
}

std::vector<double>
SurvProb::Getsinsqrtheta12s() const {
    return fsinsqrtheta12s;
}

//std::vector<double>
//SurvProb::Getsinsqrtheta13s() const {
//    return fsinsqrtheta13s;
//}

double
SurvProb::GetCdfCutOff() const{
    return fCdfCutOff;
}

void 
SurvProb::SetCdfCutOff(double cutOff_){
    fCdfCutOff = cutOff_;
}

int 
SurvProb::GetNDims() const{
    return fNDims;
}

/////////////////
// Probability //
/////////////////

//need to have Energy as the variable, distance is a fixed parameter depending on Reactor, PDF name
double 
SurvProb::operator() (const std::vector<double>& vals_) const{
    if (vals_.size() != GetNDims())
        throw DimensionError("SurvProb::SurvProb dimensionality does not match the input vector to evaluate!");

    // double baseline = ~ GetDistance, another arugment in this function "double Distance"
    double survprob = 0;
    double dmsqr21;
    double ssqr12;
    double ssqr13;
    for(size_t i = 0; i < GetNDims(); i++){
        dmsqr21  = fdelmsqr21s.at(i);
        ssqr12 = fsinsqrtheta12s.at(i);
	ssqr13 = fsinsqrtheta13s;//.at(i);
	double fSSqr2Theta12 = pow(sin(2.0 * asin(sqrt(ssqr12))), 2.0);
	double fS4 = pow(ssqr13, 2.0);
	double fC4 = pow(1.0-ssqr13, 2.0);
	double scale = 1.267e3; // for nuE in [MeV] and baseline in [km]
	double sSqrDmBE = pow(sin(scale * dmsqr21 * fBaseline / vals_.at(i)), 2.0);

        survprob +=  (fC4 * (1.0 - fSSqr2Theta12 * sSqrDmBE) + fS4); 
    }
    return survprob;
}

double
SurvProb::Si(double val_) const{
    return (double)fsincfunction->Integral(0.,val_);
}

double SurvProb::Cdf(size_t dim_, double val_) const{

    double dmsqr21  = fdelmsqr21s.at(dim_);
    double ssqr12 = fsinsqrtheta12s.at(dim_);
    double ssqr13 = fsinsqrtheta13s;//.at(dim_);
    double scale = 1.267e3; // for nuE in [MeV] and baseline in [km]
    double a = (1.0-ssqr13)*(1.0-ssqr13);
    double b = pow(sin(2.0 * asin(sqrt(ssqr12))), 2.0);
    double c =scale * dmsqr21 * fBaseline; 
    double d = ssqr13*ssqr13;
    double e = 2.*c/val_;
    double rtnVal = 0;

    rtnVal = (a*b/2)*val_*cos(e);
    rtnVal += val_*(a+d-(a*b)/2.);
    rtnVal += a*b*c*Si(e);
    rtnVal -= a*b*c*3.1415926535898/2.;
    return rtnVal;
}

double 
SurvProb::Integral(const std::vector<double>& mins_, const std::vector<double>& maxs_) const{
    if(mins_.size() != GetNDims() || maxs_.size() != GetNDims())
        throw DimensionError("SurvProb::SurvProb, tried to integrate over interval of wrong dimensionality");

    double integral = 1;
    for(size_t i = 0; i < mins_.size(); i++){
      integral *= ( Cdf(i, maxs_[i]) - Cdf(i, mins_[i]))/( maxs_[i] - mins_[i] );
    }
    return integral;  
}

//BL - implement this!!
std::vector<double>
SurvProb::Sample() const{
  std::vector<double> sample(GetNDims(), 0);
  for(size_t i = 0; i < GetNDims(); i++){
    // sample[i] = Rand::Gaus(fdelmsqr21s.at(i), fsinsqrtheta12s.at(i), fsinsqrtheta12s.at(i));
  }
  return sample;
}

////////////////////////
// Make it fittable   //
////////////////////////
void
SurvProb::RenameParameter(const std::string& old_, const std::string& new_){
    fFitter.RenameParameter(old_, new_);
}

void
SurvProb::SetParameter(const std::string& name_, double value_){
    fFitter.SetParameter(name_, value_);
}

double
SurvProb::GetParameter(const std::string& name_) const{
    return fFitter.GetParameter(name_);
}

void
SurvProb::SetParameters(const ParameterDict& ps_){
     fFitter.SetParameters(ps_);
}

ParameterDict
SurvProb::GetParameters() const{
    return fFitter.GetParameters();
}

size_t
SurvProb::GetParameterCount() const{
    return fFitter.GetParameterCount();
}

std::set<std::string>
SurvProb::GetParameterNames() const{
    return fFitter.GetParameterNames();
}

std::vector<std::string>
SurvProb::Getdelmsqr21Names() const{
    return fFitter.GetDelmsqr21Names();
}

std::vector<std::string>
SurvProb::Getsinsqrtheta12Names() const{
    return fFitter.GetSinsqrtheta12Names();
}

//std::vector<std::string>
//SurvProb::Getsinsqrtheta13Names() const{
//    return fFitter.GetSinsqrtheta13Names();
//}

std::string
SurvProb::GetName() const{
    return fName;
}

void
SurvProb::SetName(const std::string& name_){
    fName = name_;
}
