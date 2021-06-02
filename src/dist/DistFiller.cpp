#include <DistFiller.h>
#include <CutCollection.h>
#include <CutLog.h>
#include <BinnedED.h>
#include <DataSet.h>
#include <Event.h>
#include <iostream>

void
DistFiller::FillDist(BinnedED& pdf_, const DataSet& data_, 
                     const CutCollection& cuts_,
                     EventSystematicManager sysMan_,                      
                     int nEv_){
  if(nEv_ < 0)
    nEv_ = data_.GetNEntries();
  for(int i = 0; i < nEv_; i++){
    if(!(i% 10000000))
      std::cout << i << "/" << nEv_ << std::endl;
    Event ev = sysMan_.ApplySystematics(data_.GetEntry(i)); 
    if(cuts_.PassesCuts(ev))
      pdf_.Fill(ev);
  }
}

void
DistFiller::FillDist(BinnedED& pdf_, const DataSet& data_, 
                     const CutCollection& cuts_,
                     CutLog& log_, 
                     EventSystematicManager sysMan_,                      
                     int nEv_){
  if(nEv_ < 0)
    nEv_ = data_.GetNEntries();
  for(int i = 0; i < nEv_; i++){
    if(!(i% 10000000))
      std::cout << i << "/" << nEv_ << std::endl;
    Event ev = sysMan_.ApplySystematics(data_.GetEntry(i)); 
    if(cuts_.PassesCuts(ev, log_))
      pdf_.Fill(ev);
  }
}
