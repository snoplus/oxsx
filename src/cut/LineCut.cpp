#include <LineCut.h>
#include <Event.h>
#include <Exceptions.h>

LineCut::LineCut(const std::string& name_, const std::string& obs_, double val_, 
                 const std::string& side_){
    fName = name_;
    fObs = obs_;
    fVal = val_;
    SetSidedness(side_);
}

Cut* 
LineCut::Clone() const{
    return new LineCut(*this);
}

bool 
LineCut::PassesCut(const Event& ev_) const{
    double val = 0;
    try{
        val = ev_.GetDatum(fObs);
    }
    catch(const NotFoundError&){
        throw NotFoundError("Cut::Cut to non-existent data observable requested!");
    }
    
    if (fSidedness == "lower")
        return val > fVal;
    else
        return val < fVal;
}
    
void   
LineCut::SetCutValue(double val_){
    fVal = val_;
}

double 
LineCut::GetCutValue() const{
    return fVal;
}

void 
LineCut::SetSidedness(const std::string& side_){
    if(side_ == "upper" || side_ == "lower")
        fSidedness = side_;
    else
        throw ValueError("LineCut::Invalid sidedness! should be 'upper' or 'lower'");
} 

std::string 
LineCut::GetSidedness() const{
    return fSidedness;
}


std::string
LineCut::GetName() const{
  return fName;
}
