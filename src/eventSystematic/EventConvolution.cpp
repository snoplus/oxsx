#include <EventConvolution.h>
#include <PDF.h>
#include <JumpPDF.h>
#include <Exceptions.h>
#include <Rand.h>
#include <Event.h>
#include <iostream>

/////////////////////////////
// CONSTRUCTORS/DESTRUCTOR //
/////////////////////////////

EventConvolution::~EventConvolution(){
  delete fDist;
}

EventConvolution::EventConvolution(const EventConvolution& other_){
  if(other_.fDist)
    fDist = other_.fDist -> Clone();
  else
    fDist = NULL;
  fInObservables = other_.fInObservables;
  fOutObservables = other_.fOutObservables;
}

EventConvolution
EventConvolution::operator=(const EventConvolution& other_){ 
  if(other_.fDist)
    fDist = other_.fDist -> Clone();
  else                                 
    fDist = NULL;

  fInObservables = other_.fInObservables;
  fOutObservables = other_.fOutObservables;
  return *this;
}

/////////////
// Get/Set //
/////////////

void
EventConvolution::SetPDF(PDF* f_){
  if(!f_)
    fDist = NULL;

  else{
    if(f_->GetNDims() != 1)
      throw DimensionError("EventConvolution::SetFunction", 1, 
                           f_->GetNDims(), "Only implemented for 1D functions");
    fDist = static_cast<ConditionalPDF*>(new JumpPDF("kernel", f_));
  } 
}

void
EventConvolution::SetConditionalPDF(ConditionalPDF* c_){
  if(!c_)
    fDist = NULL;
  else
    fDist  = c_->Clone();
}

/////////////////////////////////////////////////////////////////////////////
// FIT COMPONENT INTERFACE : just forward the call to the underlying funcn //
/////////////////////////////////////////////////////////////////////////////

void
EventConvolution::RenameParameter(const std::string& old_, const std::string& new_){
    if(!fDist)
        throw NULLPointerAccessError("EventConvolution::GetParameterNames", 
                                     "Have you set the sampling function?");
    fDist->RenameParameter(old_, new_);
}

void
EventConvolution::SetParameter(const std::string& name_, double value_){
    if(!fDist)
        throw NULLPointerAccessError("EventConvolution::GetParameterNames", 
                                     "Have you set the sampling function?");
    fDist->SetParameter(name_, value_);
}

double
EventConvolution::GetParameter(const std::string& name_) const{
    if(!fDist)
        throw NULLPointerAccessError("EventConvolution::GetParameterNames", 
                                     "Have you set the sampling function?");
    return fDist->GetParameter(name_);
}

void
EventConvolution::SetParameters(const ParameterDict& ps_){
    if(!fDist)
        throw NULLPointerAccessError("EventConvolution::GetParameterNames", 
                                     "Have you set the sampling function?");
    try{
        fDist->SetParameters(ps_);
    }
    catch(const ParameterError& e_){
        throw ParameterError("EventConvolution internal function: " + std::string(e_.what()));
    }
}

ParameterDict
EventConvolution::GetParameters() const{
    if(!fDist)
        throw NULLPointerAccessError("EventConvolution::GetParameterNames", 
                                     "Have you set the sampling function?");
    return fDist->GetParameters();
}

size_t
EventConvolution::GetParameterCount() const{
    if(!fDist)
        throw NULLPointerAccessError("EventConvolution::GetParameterNames", 
                                     "Have you set the sampling function?");
    return fDist->GetParameterCount();
}

std::set<std::string>
EventConvolution::GetParameterNames() const{
    if(!fDist)
        throw NULLPointerAccessError("EventConvolution::GetParameterNames", 
                                     "Have you set the sampling function?");
    return fDist->GetParameterNames();
}

std::string
EventConvolution::GetName() const{
    return fName;
}
void
EventConvolution::SetName(const std::string& n_){
    fName = n_;
}

// Event Systematic Interface
Event
EventConvolution::operator()(const Event& event_){
  if(!fDist)
    throw NULLPointerAccessError("EventConvolution::operator()", 
                                 "Have you set the sampling function?");

  std::string obsToActOn = fOutObservables.GetNames().at(0);
  // change below once we support more than 1D
  double newVal = fDist->Sample(event_.ToObsSet(fOutObservables)).at(0);

  Event newEv(event_);
  newEv.SetDatum(obsToActOn, newVal);

  return newEv;
}
