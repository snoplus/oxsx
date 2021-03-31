#include <VaryingCDF.h>
#include <Exceptions.h>
#include <ContainerTools.hpp>
#include <ContainerParameter.h>
#include <ContainerTools.hpp>
#include <Formatter.hpp>
#include <Rand.h>
#include <PDF.h>
#include <JumpPDF.h>

using ContainerTools::ToString;
using ContainerTools::GetKeys;
using ContainerTools::GetValues;

/////////////////////////
// Constructory Things //
/////////////////////////

VaryingCDF::VaryingCDF(const std::string& name_ = ""){
    if(name_=="")
        fName=std::string("VaryingCDF");
    else
        fName=name_;
}
VaryingCDF::~VaryingCDF(){
    delete fCdf;
    for(std::map<std::string,Function*>::const_iterator parameter_ = fFunctions.begin(); parameter_ !=fFunctions.end(); ++parameter_) {
        delete parameter_->second;
    }
}

VaryingCDF::VaryingCDF(const VaryingCDF& copy_){
    fFunctions = copy_.fFunctions;
    fCdf = copy_.fCdf->Clone();
    fName = copy_.fName;
}

VaryingCDF&
VaryingCDF::operator=(const VaryingCDF& copy_){
    fFunctions = copy_.fFunctions;
    fCdf = copy_.fCdf->Clone();
    fName = copy_.fName;
    return *this;
}

ConditionalPDF*
VaryingCDF::Clone() const{
    return static_cast<ConditionalPDF*> (new VaryingCDF(*this));
}

ParameterDict
VaryingCDF::SetupParameters(const std::vector<double>& x2_) const {
    ParameterDict parameters;
    for(std::map<std::string,Function*>::const_iterator parameter_ = fFunctions.begin(); parameter_ !=fFunctions.end(); ++parameter_) {
        parameters[parameter_->first] = (parameter_->second)->operator()(x2_);
    }
    return parameters;
}

double
VaryingCDF::ConditionalProbability(const std::vector<double>& x_, 
                                   const std::vector<double>& x2_){
    if (!fCdf)
        throw NULLPointerAccessError(Formatter()<<"VaryingCDF:: fCdf is not pointing to a PDF.");
    fCdf->SetParameters(SetupParameters(x2_));
    return fCdf->ConditionalProbability(x_,x2_);
}

void
VaryingCDF::SetDependance(const std::string& paraName_, const Function* func_){
    fFunctions[paraName_] = func_->Clone();
}

double
VaryingCDF::Integral(const std::vector<double>& mins_, 
                              const std::vector<double>& maxs_,
                              const std::vector<double>& x2_) const{
    if (!fCdf)
        throw NULLPointerAccessError(Formatter()<<"VaryingCDF:: fCdf is not pointing to a ConditionalPDF.");
    fCdf->SetParameters(SetupParameters(x2_));
    return fCdf->Integral(mins_,maxs_,x2_);
}

void 
VaryingCDF::SetKernel(ConditionalPDF* PDF_){
    fCdf=PDF_->Clone();
}

void 
VaryingCDF::SetKernel(PDF* PDF_){
    // BL : If we are passed a PDF we transform to a JumpPDF in order to preserve
    // the structure of shifting a function to a bin center and integrating over
    // over bins relative to that shift.
    fCdf=static_cast<ConditionalPDF*>(new JumpPDF("kernel",PDF_));
}

std::vector<double>
VaryingCDF::Sample(const std::vector<double>& x2_) const{
    if (!fCdf)
        throw NULLPointerAccessError(Formatter()<<"VaryingCDF:: fCdf is not pointing to a ConditionalPDF.");
    fCdf->SetParameters(SetupParameters(x2_));
    return fCdf->Sample(x2_);
}

void 
VaryingCDF::SetParameter(const std::string& name_, double value_){
    for (std::map<std::string,Function*>::iterator function = fFunctions.begin(); function != fFunctions.end(); ++function) {
        function->second->SetParameter(name_,value_);
    }
}

double 
VaryingCDF::GetParameter(const std::string& name_) const{
    //This returns parameters of the parameter functions that have been defined.
    size_t changed = 0;
    double holder=0;
    for (std::map<std::string,Function*>::const_iterator function = fFunctions.begin(); function != fFunctions.end(); ++function) {
        try {
            holder = function->second->GetParameter(name_);
            changed = 1;
        }catch(const NotFoundError& e_) {
            ;
        }
        return holder;
    }
    if (changed !=1){
            std::vector<std::string> functionNames;
            for (std::map<std::string,Function*>::const_iterator function = fFunctions.begin(); function != fFunctions.end(); ++function) {
                functionNames.push_back(function->second->GetName());
            }
            throw NotFoundError(Formatter() << "VaryingCDF:: Parameter : " << name_ <<" not found in fFunctions available. Functions available: "
                    << ToString(functionNames)<<" Parameters available : "<<ToString(GetParameterNames()) ); 
    }
}

void
VaryingCDF::SetParameters(const ParameterDict& paraDict_){
    for (std::map<std::string,Function*>::const_iterator function = fFunctions.begin(); function != fFunctions.end(); ++function) {
        function->second->SetParameters(paraDict_);
    }
}

ParameterDict
VaryingCDF::GetParameters() const{
    ParameterDict paraDict_;
    for (std::map<std::string,Function*>::const_iterator function = fFunctions.begin(); function != fFunctions.end(); ++function) {
        ParameterDict holder = function->second->GetParameters();
        paraDict_.insert(holder.begin(), holder.end());
    }
    return paraDict_;
}

size_t
VaryingCDF::GetParameterCount() const{
    size_t number=0;
    for (std::map<std::string,Function*>::const_iterator function = fFunctions.begin(); function != fFunctions.end(); ++function) {
        number += function->second->GetParameterCount();
    }
    return number;
}

std::set<std::string>
VaryingCDF::GetParameterNames() const {
    std::set<std::string> names_;
    for (std::map<std::string,Function*>::const_iterator function = fFunctions.begin(); function != fFunctions.end(); ++function) {
        std::set<std::string> holder = function->second->GetParameterNames();
        names_.insert(holder.begin(), holder.end());
    }
    return names_;
}

void
VaryingCDF::RenameParameter(const std::string& old_, const std::string& new_){
    for (std::map<std::string,Function*>::const_iterator function = fFunctions.begin(); function != fFunctions.end(); ++function) {
        function->second->RenameParameter(old_,new_);
    }
}

std::string
VaryingCDF::GetName() const{
    return fName;   
}

void 
VaryingCDF::SetName(const std::string& name_){
    fName= name_;
}
