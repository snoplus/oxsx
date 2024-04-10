#include <CutCollection.h>
#include <Event.h>
#include <Cut.h>
#include <CutLog.h>
#include <Formatter.hpp>
#include <iostream>

bool 
CutCollection::PassesCuts(const Event& event_) const{
  for(size_t i = 0; i < fCuts.size(); i++)
    if (!fCuts.at(i)->PassesCut(event_))
        return false;
  return true;
}

bool CutCollection::PassesCuts(const Event& event_, CutLog& log_) const{
    for(size_t i = 0; i < fCuts.size(); i++)
        if(!fCuts.at(i) -> PassesCut(event_)){
            log_.LogCut(i);
            return false;
        }
    log_.LogPass();
    return true;
}

void
CutCollection::AddCut(const Cut& cut_){
  fCuts.push_back(cut_.Clone());
}

CutCollection::~CutCollection(){
  for(size_t i = 0; i < fCuts.size(); i++){
    delete fCuts.at(i);
  }
}

CutCollection::CutCollection(const CutCollection& other_){
  fCuts.reserve(other_.fCuts.size());
  for(size_t i = 0; i < other_.fCuts.size(); i++)
    fCuts.push_back(other_.fCuts.at(i)->Clone());
}

CutCollection&
CutCollection::operator= (const CutCollection& other_){
  for(size_t i = 0; i < fCuts.size(); i++)
    delete fCuts.at(i);

  fCuts.resize(0);
  fCuts.reserve(other_.fCuts.size());
  for(size_t i = 0; i < other_.fCuts.size(); i++)
    fCuts.push_back(other_.fCuts.at(i)->Clone());
  return *this;
}

std::vector<std::string> 
CutCollection::GetCutNames() const{
  std::vector<std::string> names;
  for(size_t i = 0; i < fCuts.size(); i++)
    names.push_back(fCuts.at(i)->GetName());
  return names;
}
