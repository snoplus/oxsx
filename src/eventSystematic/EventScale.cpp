#include <EventScale.h>
#include <Exceptions.h>
#include <Event.h>
#include <string>
#include <Formatter.hpp>
#include <ContainerTools.hpp>

using ContainerTools::ToString;
using ContainerTools::GetKeys;

Event
EventScale::operator()(const Event& inEvent_){
    Event newEvent = inEvent_;

    std::string obsToChange = fOutObservables.GetNames().at(0);

    // pull out the relevant data point
    newEvent.SetDatum(obsToChange, inEvent_.GetDatum(obsToChange) * GetScale());
    return newEvent;
}

void 
EventScale::SetScale(double scale_){
    if(scale_ <= 0)
        throw ValueError(Formatter() << "Tried to set EventScale systematic scale to " << scale_ << " - must be >= 0");
    fScale = scale_;
}

double
EventScale::GetScale() const{
    return fScale;
}


// Fit Component Interface
void
EventScale::SetParameter(const std::string& name_, double value){
    if(name_ != fParamName)
        throw ParameterError("EventScale: can't set " + name_ + ", " + fParamName + " is the only parameter" );
    fScale = value;
}

double
EventScale::GetParameter(const std::string& name_) const{
   if(name_ != fParamName)
        throw ParameterError("EventScale: can't get " + name_ + ", " + fParamName + " is the only parameter" );
   return fScale;
}

void
EventScale::SetParameters(const ParameterDict& pd_){
    try{
        fScale = pd_.at(fParamName);
    }
    catch(const std::out_of_range& e_){
        throw ParameterError("Set dictionary is missing " + fParamName + ". It did contain: \n" + ToString(GetKeys(pd_)));
    }
}

ParameterDict
EventScale::GetParameters() const{
    ParameterDict d;
    d[fParamName] = fScale;
    return d;
}

size_t
EventScale::GetParameterCount() const{
    return 1;
}

std::set<std::string>
EventScale::GetParameterNames() const{
    std::set<std::string> ret;
    ret.insert(fParamName);
    return ret;
}

void
EventScale::RenameParameter(const std::string& old_, const std::string& new_){
    if(old_ != fParamName)
        throw ParameterError("EventScale: can't rename " + old_ + ", " + fParamName + " is the only parameter" );
    fParamName = new_;
}

std::string
EventScale::GetName() const{
    return fName;
}

void
EventScale::SetName(const std::string& name_){
    fName = name_;
}



