#include <EventReconvolution.h>
#include <Event.h>
#include <Formatter.hpp>
#include <Exceptions.h>
#include <iostream>
#include <ContainerTools.hpp>

using ContainerTools::ToString;
using ContainerTools::GetKeys;

Event
EventReconvolution::operator()(const Event& event_){
    std::string truthObs   = fInObservables.GetNames().at(0);
    std::string smearedObs = fOutObservables.GetNames().at(0);


    double truthVal  = event_.GetDatum(truthObs);;  //  eg smeared energy
    double smearedVal  = event_.GetDatum(smearedObs); // eg mc energy
    
    Event newEvent(event_);
    newEvent.SetDatum(smearedObs,   truthVal + fCorrection * (smearedVal - truthVal) );

    return newEvent;
}

void
EventReconvolution::SetCorrection(double correction_){
    if(correction_ < 0)
        throw ValueError(Formatter() << "Tried to set Reconvolution systematic correction factor to " << correction_ << " - must be >= 0");
    fCorrection = correction_;
}
    

// Fit Component Interface

void
EventReconvolution::SetParameter(const std::string& name_, double value){
    if(name_ != fParamName)
        throw ParameterError("EventReconvolution: can't set " + name_ + ", " + fParamName + " is the only parameter" );
    fCorrection = value;
}

double
EventReconvolution::GetParameter(const std::string& name_) const{
   if(name_ != fParamName)
        throw ParameterError("EventReconvolution: can't get " + name_ + ", " + fParamName + " is the only parameter" );
   return fCorrection;
}

void
EventReconvolution::SetParameters(const ParameterDict& pd_){
    try{
        fCorrection = pd_.at(fParamName);
    }
    catch(const std::out_of_range& e_){
        throw ParameterError("Set dictionary is missing " + fParamName + ". I did contain: \n" + ToString(GetKeys(pd_)));
    }
}

ParameterDict
EventReconvolution::GetParameters() const{
    ParameterDict d;
    d[fParamName] = fCorrection;
    return d;
}

size_t
EventReconvolution::GetParameterCount() const{
    return 1;
}

std::set<std::string>
EventReconvolution::GetParameterNames() const{
    std::set<std::string> ret;
    ret.insert(fParamName);
    return ret;
}

void
EventReconvolution::RenameParameter(const std::string& old_, const std::string& new_){
    if(old_ != fParamName)
        throw ParameterError("EventReconvolution: can't rename " + old_ + ", " + fParamName + " is the only parameter" );
    fParamName = new_;
}

std::string
EventReconvolution::GetName() const{
    return fName;
}

void
EventReconvolution::SetName(const std::string& name_){
    fName = name_;
}
