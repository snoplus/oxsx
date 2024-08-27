/**********************************************************************************************/
/* OXSX own data object, can be written and read from HTF5 compressed files                   */
/* Simply a wrapper on a vector of Event objects, dynamically allocated to stop the stack */
/* filling up.                                                                                */
/**********************************************************************************************/

#ifndef __OXSX_LAZY_OXSX_DATA_SET__
#define __OXSX_LAZY_OXSX_DATA_SET__
#include <OXSXDataSet.h>
#include <Event.h>
#include <string>

class LazyOXSXDataSet : public DataSet
{
public:
    LazyOXSXDataSet(const std::string &filename_) : fFileName(filename_), fDataSet(NULL) { Register(); }
    ~LazyOXSXDataSet(); // frees fData
    LazyOXSXDataSet(const LazyOXSXDataSet &);
    LazyOXSXDataSet operator=(const LazyOXSXDataSet &);

    // Memory management
    void Load() const;
    void Close() const;

    void Register();
    void Deregister();

    // DataSet interface
    unsigned GetNEntries() const;
    unsigned GetNObservables() const;
    Event GetEntry(size_t) const;
    std::vector<std::string> GetObservableNames() const;

private:
    static std::vector<LazyOXSXDataSet *> fFiles;
    std::string fFileName;
    mutable OXSXDataSet *fDataSet;
};
#endif
