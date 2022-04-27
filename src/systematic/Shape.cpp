#include <Shape.h>
#include <Exceptions.h>

void Shape::SetShapeFunction(ShapeFunction& shape_func,
                             std::vector<std::string>& param_names) {
    /*
    * Setter for the shape function object; also the only way to add parameters
    * into the parameter dict: the FitComponent interface allows one to set 
    * values or rename params, but not add params.
    */
   fShapeFunc = shape_func;
   // Clear any existing param stuff, then add param names with default values
   fParamDict.clear();
   for (auto &&param : param_names) {
       fParamDict[param] = 0;
   }
}

void Shape::Construct() {
    /*
    * 
    */
}

// FitComponent interface: handling systematic's parameters
// Note that once the parameter names have been set in SetShapeFunction(),
// The number of parameters cannot be changed with this interface.
// You can change the values and names of these parameters

void
Shape::SetParameter(const std::string& name_, double value){
    if(fParamDict.count(name_) == 0)
        throw ParameterError("Shift: can't set " + name_);
    fParamDict[name_] = value;
}

double
Shape::GetParameter(const std::string& name_) const{
   if(fParamDict.count(name_) == 0)
        throw ParameterError("Shift: can't set " + name_);
    return fParamDict.at(name_);
}

void
Shape::SetParameters(const ParameterDict& pd_){
    try{
        fScaleFactor = pd_.at(fParamName);
    }
    catch(const std::out_of_range& e_){
        throw ParameterError("Set dictionary is missing " + fParamName + ". I did contain: \n" + ContainerTools::ToString(ContainerTools::GetKeys(pd_)));
    }
}

std::set<std::string>
Shape::GetParameterNames() const {
    std::set<std::string> set;
    for (auto &&pair : fParamDict) {
        set.insert(pair.first);
    }
    return set;
}

void
Shape::RenameParameter(const std::string& old_, const std::string& new_){
    if(old_ != fParamName)
        throw ParameterError("Scale: can't rename " + old_ + ", " + fParamName + " is the only parameter" );
    fParamName = new_;
}