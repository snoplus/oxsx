#ifndef __OXSX_HIST_TOOLS__
#define __OXSX_HIST_TOOLS__
#include <vector>
#include <string>
#include <map>
#include <set>


class AxisCollection;
class Histogram;

class HistTools{
 public:
    static std::map<std::string, Histogram> MakeAllHists(const AxisCollection&, 
                                                         const std::set<std::pair<std::string, std::string> >& combinations_);

    template<typename T1>
    static std::map<std::string, Histogram> MakeAllHists(const AxisCollection&, 
                                                         const T1& names_);

    static void FillAllHists(std::map<std::string, Histogram>&, const std::map<std::string, double>& fillvals);


    // this works for histograms and binnedED -> put it somewhere else
    template<typename T>
    static std::vector<T> GetVisualisableProjections(const T&);

    template<typename T>
    static std::vector<T> Get1DSlices(const T&, const std::string&);
};

#include <HistTools.hpp>
#endif
