#ifndef __OXSX_CUT_LOG__
#define __OXSX_CUT_LOG__
#include <vector>
#include <string>

class CutLog{
 public:
    CutLog() : fNEvents(0) {}
    CutLog(const std::vector<std::string> cutNames_) : fCutCounts(std::vector<int>(cutNames_.size(), 0)), fCutNames(cutNames_), fNEvents(0) {}    
       
    void LogCut(size_t cutIndex_);
    void LogPass();

    std::vector<int>         GetCutCounts() const;
    std::vector<int>         GetRemainderCounts() const;
    std::vector<double>      GetCutPercentages() const;
    std::vector<double>      GetRemainderPercentages() const;
    std::vector<std::string> GetCutNames() const;
    
    std::string AsString() const;
    void CalculateMeta();
    void Print() const;

 private:
    std::vector<int>         fCutCounts;
    std::vector<int>         fRemainderCounts;
    std::vector<double>      fCutPercentages;
    std::vector<double>      fRemainderPercentages;

    std::vector<std::string> fCutNames;    
    int fNEvents;
};
#endif


