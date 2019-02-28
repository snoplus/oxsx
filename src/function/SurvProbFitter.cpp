#include <SurvProbFitter.h>
#include <SurvProb.h>
#include <iostream>
#include <sstream>
#include <Formatter.hpp>
#include <ContainerTools.hpp>
#include <Exceptions.h>
#include <Formatter.hpp>
#include <algorithm>

using ContainerTools::ToString;

SurvProbFitter::SurvProbFitter(SurvProb* gaus, const size_t& nDims_){
    fOrignalFunc = gaus; 
    std::stringstream ss;
    for (int i = 0; i <nDims_; ++i) {
        ss << "delmsqr21" << "_" << i;
        fDelmsqr21sNames.push_back(ss.str());
        ss.str("");
        ss << "sinsqrtheta12" << "_" << i;
        fSinsqrtheta12sNames.push_back(ss.str());
        ss.str("");
        //ss << "sinsqrtheta13" << "_" << i;
        //fSinsqrtheta13sNames.push_back(ss.str());
        //ss.str("");
    }
}

SurvProbFitter::SurvProbFitter(SurvProb* gaus, const std::vector<std::string>& delmsqr21Names_, const std::vector<std::string>& Sinsqrtheta12Names_){ //,const std::vector<std::string>& Sinsqrtheta13Names_
    if (delmsqr21Names_.size()!= Sinsqrtheta12Names_.size())
        throw OXSXException(Formatter()<<"SurvProbFitter:: #delmsqr21Name != #Sinsqrtheta12Names");
    //if (Sinsqrtheta13Names_.size()!= Sinsqrtheta12Names_.size())
    //    throw OXSXException(Formatter()<<"SurvProbFitter:: #delmsqr21Name != #Sinsqrtheta12Names");
    fOrignalFunc = gaus;
    std::stringstream ss;
    for (int i = 0; i < delmsqr21Names_.size(); ++i)
        fDelmsqr21sNames.push_back(delmsqr21Names_.at(i));
    for (int i = 0; i <Sinsqrtheta12Names_.size(); ++i)
        fSinsqrtheta12sNames.push_back(Sinsqrtheta12Names_.at(i));
    //for (int i = 0; i <Sinsqrtheta13Names_.size(); ++i)
    //    fSinsqrtheta13sNames.push_back(Sinsqrtheta13Names_.at(i));
}

void
SurvProbFitter::RenameParameter(const std::string& old_, const std::string& new_){
    std::vector<std::string>::iterator it;
    if(find(fDelmsqr21sNames.begin(),fDelmsqr21sNames.end(),old_)==fDelmsqr21sNames.end() && find(fSinsqrtheta12sNames.begin(),fSinsqrtheta12sNames.end(),old_)==fSinsqrtheta12sNames.end() )
	//&& find(fSinsqrtheta13sNames.begin(),fSinsqrtheta13sNames.end(),old_)==fSinsqrtheta13sNames.end()
            throw NotFoundError(Formatter()<<"SurvProbFitter:: When attempting to renaming the parameter "<< old_<<", it wasn't found. Available names: "<<
                    ToString(GetParameterNames()) );
 
    it=find(fDelmsqr21sNames.begin(),fDelmsqr21sNames.end(),old_);
    while(it!=fDelmsqr21sNames.end()){
        *it=new_;
        it=find(it++,fDelmsqr21sNames.end(),old_);
    }

    it=find(fSinsqrtheta12sNames.begin(),fSinsqrtheta12sNames.end(), old_);
    while(it!=fSinsqrtheta12sNames.end()){
        *it=new_;
        it=find(it++,fSinsqrtheta12sNames.end(),old_);
    }

    //it=find(fSinsqrtheta13sNames.begin(),fSinsqrtheta13sNames.end(), old_);
    //while(it!=fSinsqrtheta13sNames.end()){
    //    *it=new_;
    //    it=find(it++,fSinsqrtheta13sNames.end(),old_);
    //}
}

std::vector<std::string>
SurvProbFitter::GetDelmsqr21Names() const{
    return fDelmsqr21sNames; 
}

std::vector<std::string>
SurvProbFitter::GetSinsqrtheta12Names() const{
    return fSinsqrtheta12sNames; 
}

//std::vector<std::string>
//SurvProbFitter::GetSinsqrtheta13Names() const{
//    return fSinsqrtheta13sNames; 
//}

void
SurvProbFitter::SetParameter(const std::string& name_, double value_){
    std::vector<std::string>::iterator it;
    it=find(fDelmsqr21sNames.begin(),fDelmsqr21sNames.end(), name_);
    while(it!=fDelmsqr21sNames.end()){
        fOrignalFunc->Setdelmsqr21(it-fDelmsqr21sNames.begin(),value_);
        it=find(++it,fDelmsqr21sNames.end(), name_);
    }
    it=find(fSinsqrtheta12sNames.begin(),fSinsqrtheta12sNames.end(), name_);
    while(it!=fSinsqrtheta12sNames.end()){
        fOrignalFunc->Setsinsqrtheta12(it-fSinsqrtheta12sNames.begin(),value_);
        it=find(++it,fSinsqrtheta12sNames.end(), name_);
    }
    //it=find(fSinsqrtheta13sNames.begin(),fSinsqrtheta13sNames.end(), name_);
    //while(it!=fSinsqrtheta13sNames.end()){
    //    fOrignalFunc->Setsinsqrtheta13(it-fSinsqrtheta13sNames.begin(),value_);
    //    it=find(++it,fSinsqrtheta13sNames.end(), name_);
    //}
}

double
SurvProbFitter::GetParameter(const std::string& name_) const{
    // BL: If n parameters have the same name (either across both delmsqr21s
    // and stddevs or not) the value will be the same. So the value of the
    // first instance is sufficient.

    std::vector<std::string>::const_iterator it;
    it=find(fDelmsqr21sNames.begin(),fDelmsqr21sNames.end(), name_);
    if(it==fDelmsqr21sNames.end()){
      
        it=find(fSinsqrtheta12sNames.begin(),fSinsqrtheta12sNames.end(), name_);

        if(it==fDelmsqr21sNames.end()){
            //it=find(fSinsqrtheta13sNames.begin(),fSinsqrtheta13sNames.end(), name_);

            //if(it==fSinsqrtheta13sNames.end())
            //    throw NotFoundError(Formatter()<<"SurvProbFitter:: Parameter : "<<
            //                        name_<<
            //                        " was not known to the SurvProbFitter. Available names: "<<
            //                        ToString(GetParameterNames()) );
            //return fOrignalFunc->Getsinsqrtheta13(it-fSinsqrtheta13sNames.end());
        }
      return fOrignalFunc->Getsinsqrtheta12(it-fSinsqrtheta12sNames.end());
    }
    return fOrignalFunc->Getdelmsqr21(it-fDelmsqr21sNames.end());
}

void
SurvProbFitter::SetParameters(const ParameterDict& ps_){
    for (ParameterDict::const_iterator i = ps_.begin(); i != ps_.end(); ++i) {
        SetParameter(i->first,i->second);
    }
}

ParameterDict
SurvProbFitter::GetParameters() const{
    std::vector<double> delmsqr21s = fOrignalFunc->Getdelmsqr21s();
    std::vector<double> sinsqrtheta12s= fOrignalFunc->Getsinsqrtheta12s();
    //std::vector<double> sinsqrtheta13s= fOrignalFunc->Getsinsqrtheta13s();
    std::vector<double> values;

    //values.reserve( delmsqr21s.size() + sinsqrtheta12s.size()+ sinsqrtheta13s.size() ); // preallocate memory
    values.reserve( delmsqr21s.size() + sinsqrtheta12s.size() ); // preallocate memory
    values.insert( values.end(), delmsqr21s.begin(), delmsqr21s.end() );
    values.insert( values.end(), sinsqrtheta12s.begin(), sinsqrtheta12s.end() );
    //values.insert( values.end(), sinsqrtheta13s.begin(), sinsqrtheta13s.end() );

    std::vector<std::string> names;
    //names.reserve( fDelmsqr21sNames.size() + fSinsqrtheta12sNames.size() + fSinsqrtheta13sNames.size() ); // preallocate memory
    names.reserve( fDelmsqr21sNames.size() + fSinsqrtheta12sNames.size() ); // preallocate memory
    names.insert( names.end(), fDelmsqr21sNames.begin(), fDelmsqr21sNames.end() );
    names.insert( names.end(), fSinsqrtheta12sNames.begin(), fSinsqrtheta12sNames.end() );
    //names.insert( names.end(), fSinsqrtheta13sNames.begin(), fSinsqrtheta13sNames.end() );

    return ContainerTools::CreateMap(names,values);
}

size_t
SurvProbFitter::GetParameterCount() const{
    return fDelmsqr21sNames.size()+fSinsqrtheta12sNames.size(); //+fSinsqrtheta13sNames.size()
}

std::set<std::string>
SurvProbFitter::GetParameterNames() const{
    std::set<std::string> names;
    for (int i = 0; i < fDelmsqr21sNames.size(); ++i){
        names.insert(fDelmsqr21sNames.at(i));
        names.insert(fSinsqrtheta12sNames.at(i));
        //names.insert(fSinsqrtheta13sNames.at(i));
    }
    return names;
}
