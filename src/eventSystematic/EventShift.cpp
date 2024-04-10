#include <EventShift.h>
#include <Event.h>
#include <Exceptions.h>
#include <ContainerTools.hpp>

using ContainerTools::ToString;
using ContainerTools::GetKeys;

Event
EventShift::operator()(const Event& inEvent_){
    Event newEvent = inEvent_;
    
    std::string obsToChange = fOutObservables.GetNames().at(0);
    
    // pull out the relevant data point
    newEvent.SetDatum(obsToChange, inEvent_.GetDatum(obsToChange) + GetShift());
    return newEvent;
}

void
EventShift::SetShift(double val_){
    fShift = val_;
}

double 
EventShift::GetShift() const{
    return fShift;
}

// Fit Component Interface
void
EventShift::SetParameter(const std::string& name_, double value){
    if(name_ != fParamName)
        throw ParameterError("EventShift: can't set " + name_ + ", " + fParamName + " is the only parameter" );
    fShift = value;
}

double
EventShift::GetParameter(const std::string& name_) const{
   if(name_ != fParamName)
        throw ParameterError("EventShift: can't get " + name_ + ", " + fParamName + " is the only parameter" );
   return fShift;
}

void
EventShift::SetParameters(const ParameterDict& pd_){
    try{
        fShift = pd_.at(fParamName);
    }
    catch(const std::out_of_range& e_){
        throw ParameterError("Set dictionary is missing " + fParamName + ". I did contain: \n" + ToString(GetKeys(pd_)));
    }
}

ParameterDict
EventShift::GetParameters() const{
    ParameterDict d;
    d[fParamName] = fShift;
    return d;
}

size_t
EventShift::GetParameterCount() const{
    return 1;
}

std::set<std::string>
EventShift::GetParameterNames() const{
    std::set<std::string> set;
    set.insert(fParamName);
    return set;
}

void
EventShift::RenameParameter(const std::string& old_, const std::string& new_){
    if(old_ != fParamName)
        throw ParameterError("EventShift: can't rename " + old_ + ", " + fParamName + " is the only parameter" );
    fParamName = new_;
}

std::string
EventShift::GetName() const{
    return fName;
}

void
EventShift::SetName(const std::string& name_){
    fName = name_;
}

