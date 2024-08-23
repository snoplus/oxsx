/**********************************************************************************************/
/* OXSX own data object, can be written and read from HTF5 compressed files                   */
/* Simply a wrapper on a vector of Event objects, dynamically allocated to stop the stack */
/* filling up.                                                                                */
/**********************************************************************************************/

#ifndef __OXSX_OXSX_DATA_SET__
#define __OXSX_OXSX_DATA_SET__
#include <vector>
#include <DataSet.h>
#include <Event.h>
#include <string>

class OXSXDataSet : public DataSet
{
public:
    OXSXDataSet() : fNObservables(0) {}

    void AddEntry(const Event &);
    Event GetEntry(size_t eventIndex_) const;
    unsigned GetNEntries() const;
    unsigned GetNObservables() const;

    void SetObservableNames(const std::vector<std::string> &names_);
    std::vector<std::string> GetObservableNames() const;

    // concatenate
    OXSXDataSet operator+(const OXSXDataSet &other_);

    void Reserve(int);

private:
    std::vector<Event> fData;
    std::vector<std::string> fObservableNames;
    unsigned fNObservables;
};
#endif
