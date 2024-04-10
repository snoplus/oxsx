#include <JumpPDF.h>
#include <PDF.h>
#include <Exceptions.h>
#include <algorithm>
#include <functional>

//////////////////////////////
// CONSTRUCTORS/DESTRUCTORS //
//////////////////////////////

JumpPDF::JumpPDF(const std::string& name_, PDF* f_){
  fName = name_;
  if(f_)
    fPDF = static_cast<PDF*>(f_->Clone());
  
  else
    fPDF = NULL;  
}

JumpPDF::~JumpPDF(){
  delete fPDF;
}

JumpPDF::JumpPDF(const JumpPDF& other_){
  if(other_.fPDF)
    fPDF = static_cast<PDF*>(other_.fPDF->Clone());
  else
    fPDF = NULL;

  fName = other_.fName;
}

JumpPDF
JumpPDF::operator=(const JumpPDF& other_){
  if(other_.fPDF)
    fPDF = static_cast<PDF*>(other_.fPDF->Clone());
  else
    fPDF = NULL;

  fName = other_.fName;
  return *this;
}

std::vector<double>
JumpPDF::Diff(const std::vector<double>& x_, 
                                const std::vector<double>& x2_) const{
  std::vector<double> diff = x_;
  std::transform(diff.begin(), diff.end(), x2_.begin(), 
                 diff.begin(), std::minus<double>());
  return diff;
}

std::vector<double>
JumpPDF::Sum(const std::vector<double>& x_, 
                               const std::vector<double>& x2_) const{
  std::vector<double> sum = x_;
  std::transform(sum.begin(), sum.end(), x2_.begin(), 
                 sum.begin(), std::plus<double>());
  return sum;
}

ConditionalPDF*
JumpPDF::Clone() const{
  return static_cast<ConditionalPDF*>(new JumpPDF(*this));
}

////////////////////////////////////////
// CONDITIONAL DISTRIBUTION INTERFACE //
////////////////////////////////////////

double 
JumpPDF::ConditionalProbability(const std::vector<double>& x_,
                                 const std::vector<double>& x2_){
  if(!fPDF)
    throw NULLPointerAccessError("JumpPDF::Probability",
                                 "Have you set the function?");
     
  return fPDF->operator()(Diff(x_, x2_));
}

std::vector<double> 
JumpPDF::Sample(const std::vector<double>& x2_) const{
  if(!fPDF)
    throw NULLPointerAccessError("JumpPDF::Integral",
                                 "Have you set the function?");
  std::vector<double> samp = fPDF->Sample();
  return Sum(samp, x2_);
}
double
JumpPDF::Integral(const std::vector<double>& mins_,
                                    const std::vector<double>& maxs_,
                                    const std::vector<double>& x2_) const{
  if(!fPDF)
    throw NULLPointerAccessError("JumpPDF::Integral",
                                 "Have you set the function?");
  
  return fPDF->Integral( Diff(mins_, x2_) , Diff(maxs_, x2_) );
}
                                      

/////////////////////////////
// FIT COMPONENT INTERFACE //
/////////////////////////////

std::string
JumpPDF::GetName() const{
    return fName;
}

void 
JumpPDF::SetName(const std::string& name){
    fName=name;
}
void
JumpPDF::RenameParameter(const std::string& old_, const std::string& new_){
    if(!fPDF)
         throw NULLPointerAccessError("JumpPDF::Probability",
                                      "Have you set the function?");
    fPDF->RenameParameter(old_, new_);
}

void
JumpPDF::SetParameter(const std::string& name_, double value_){
    if(!fPDF)
        throw NULLPointerAccessError("JumpPDF::Probability",
                                     "Have you set the function?");
    fPDF->SetParameter(name_, value_);
}

double
JumpPDF::GetParameter(const std::string& name_) const{
    if(!fPDF)
         throw NULLPointerAccessError("JumpPDF::Probability",
                                      "Have you set the function?");
    return fPDF->GetParameter(name_);
}

void
JumpPDF::SetParameters(const ParameterDict& ps_){
    if(!fPDF)
         throw NULLPointerAccessError("JumpPDF::Probability",
                                      "Have you set the function?");
    try{
        fPDF->SetParameters(ps_);
    }
    catch(const ParameterError& e_){
        throw ParameterError("JumpPDF internal function: " + std::string(e_.what()));
    }
}

ParameterDict
JumpPDF::GetParameters() const{
    if(!fPDF)
         throw NULLPointerAccessError("JumpPDF::Probability",
                                      "Have you set the function?");
    return fPDF->GetParameters();
}

size_t
JumpPDF::GetParameterCount() const{
    if(!fPDF)
         throw NULLPointerAccessError("JumpPDF::Probability",
                                      "Have you set the function?");
    return fPDF->GetParameterCount();
}

std::set<std::string>
JumpPDF::GetParameterNames() const{
    if(!fPDF)
         throw NULLPointerAccessError("JumpPDF::Probability",
                                      "Have you set the function?");
    return fPDF->GetParameterNames();
}
