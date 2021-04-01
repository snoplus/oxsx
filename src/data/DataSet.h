/******************************************************************/
/* Abstract base for a class that hands out data from some source */
/******************************************************************/
#ifndef __OXSX_DATA_SET__
#define __OXSX_DATA_SET__
#include <stddef.h>
#include <vector>
#include <string>

class Event;
class DataSet{
 public:
    virtual Event GetEntry(size_t iEvent_) const = 0;
    virtual ~DataSet() {}

    virtual unsigned GetNEntries() const = 0; 
    virtual unsigned GetNObservables() const = 0; 

    virtual std::vector<std::string> GetObservableNames() const = 0;

};
#endif
