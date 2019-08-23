#include <GaussianEResFitter.h>
#include <GaussianERes.h>
#include <iostream>
#include <sstream>
#include <Formatter.hpp>
#include <ContainerTools.hpp>
#include <Exceptions.h>
#include <Formatter.hpp>
#include <algorithm>

using ContainerTools::ToString;

GaussianEResFitter::GaussianEResFitter(GaussianERes* gaus, const size_t& nDims_){
    fOrignalFunc = gaus; 
    std::stringstream ss;
    for (int i = 0; i <nDims_; ++i) {
        ss.str("");
        ss << "eres" << "_" << i;
        fEResNames.push_back(ss.str());
        ss.str("");
    }
}

GaussianEResFitter::GaussianEResFitter(GaussianERes* gaus, const std::vector<std::string>& EResNames_){
    fOrignalFunc = gaus;
    std::stringstream ss;
    for (int i = 0; i <EResNames_.size(); ++i)
        fEResNames.push_back(EResNames_.at(i));
}

void
GaussianEResFitter::RenameParameter(const std::string& old_, const std::string& new_){
    std::vector<std::string>::iterator it;
    if(find(fEResNames.begin(),fEResNames.end(),old_)==fEResNames.end())
            throw NotFoundError(Formatter()<<"GaussianEResFitter:: When attempting to renaming the parameter "<< old_<<", it wasn't found. Available names: "<<
                    ToString(GetParameterNames()) );
 
    it=find(fEResNames.begin(),fEResNames.end(), old_);
    while(it!=fEResNames.end()){
        *it=new_;
        it=find(it++,fEResNames.end(),old_);
    }
}

std::vector<std::string>
GaussianEResFitter::GetEResNames() const{
    return fEResNames; 
}

void
GaussianEResFitter::SetParameter(const std::string& name_, double value_){
    std::vector<std::string>::iterator it;
    it=find(fEResNames.begin(),fEResNames.end(), name_);
    while(it!=fEResNames.end()){
        fOrignalFunc->SetERes(it-fEResNames.begin(),value_);
        it=find(++it,fEResNames.end(), name_);
    }
}

double
GaussianEResFitter::GetParameter(const std::string& name_) const{
    // BL: If n parameters have the same name (either across both means
    // and stddevs or not) the value will be the same. So the value of the
    // first instance is sufficient.

    std::vector<std::string>::const_iterator it;
    it=find(fEResNames.begin(),fEResNames.end(), name_);
    if(it==fEResNames.end())
        throw NotFoundError(Formatter()<<"GaussianEResFitter:: Parameter : "<<
                                name_<<
                                " was not known to the GaussianEResFitter. Available names: "<<
                                ToString(GetParameterNames()) );
    return fOrignalFunc->GetERes(it-fEResNames.end());
	
}

void
GaussianEResFitter::SetParameters(const ParameterDict& ps_){
    for (ParameterDict::const_iterator i = ps_.begin(); i != ps_.end(); ++i) {
        SetParameter(i->first,i->second);
    }
}

ParameterDict
GaussianEResFitter::GetParameters() const{
    std::vector<double> stddevs= fOrignalFunc->GetERess();
    std::vector<double> values;

    values.reserve( stddevs.size() ); // preallocate memory
    values.insert( values.end(), stddevs.begin(), stddevs.end() );

    std::vector<std::string> names;
    names.reserve( fEResNames.size() ); // preallocate memory
    names.insert( names.end(), fEResNames.begin(), fEResNames.end() );

    return ContainerTools::CreateMap(names,values);
}

size_t
GaussianEResFitter::GetParameterCount() const{
    return fEResNames.size();
}

std::set<std::string>
GaussianEResFitter::GetParameterNames() const{
    std::set<std::string> names;
    for (int i = 0; i < fEResNames.size(); ++i)
        names.insert(fEResNames.at(i));
    
    return names;
}

