#include <AxisCollection.h>
#include <Histogram.h>
#include <Combinations.hpp>

template<typename T1>
std::map<std::string, Histogram> 
HistTools::MakeAllHists(const AxisCollection& axes_,
                        const T1& names_){
    std::map<std::string, Histogram> hists;

    typedef typename T1::const_iterator It;
    for(It it = names_.begin(); it != names_.end(); ++it){
        AxisCollection axes;
        axes.AddAxis(axes_.GetAxis(*it));
        hists[*it] = Histogram(axes);
    }
    return hists;
}

template<typename T>
std::vector<T>
HistTools::GetVisualisableProjections(const T& hist_){
  typedef std::pair<std::string, std::string> AxisPair;
  typedef std::set<AxisPair> AxisCombs;

  std::vector<T> hists;
  size_t dim  = hist_.GetNDims();
  hists.reserve(dim + dim * (dim -1)/2);
  
  std::vector<std::string> axisNames = hist_.GetAxes().GetAxisNames();
  // 1D
  for(size_t i = 0; i < axisNames.size(); i++){
    hists.push_back(hist_.Marginalise(axisNames.at(i)));
  }

  if(hist_.GetNDims() > 2){  
      // 2D
      AxisCombs combs = Combinations::AllCombsNoDiag(axisNames);
      for(AxisCombs::iterator it = combs.begin(); it != combs.end(); ++it){
          // sadly this method takes a vector not a pair so conv
          std::vector<std::string> vec;
          vec.push_back(it->first);
          vec.push_back(it->second);
          hists.push_back(hist_.Marginalise(vec));
      }
  }
  return hists;
}

template<typename T>
std::vector<T>
HistTools::Get1DSlices(const T& hist_, const std::string& axisName_){

  std::vector<std::string> axisNames = hist_.GetAxes().GetAxisNames();
  std::vector<std::string> fixedAxisNames;
  for(std::vector<std::string>::iterator it = axisNames.begin(); it!=axisNames.end(); it++){
    if ((*it)!=axisName_) fixedAxisNames.push_back(*it);
  }

  // get no of bins for the corresponding axes
  std::vector<size_t> maxes;
  for(std::vector<std::string>::iterator it = fixedAxisNames.begin(); it!=fixedAxisNames.end(); it++){
    maxes.push_back(hist_.GetAxes().GetAxis(*it).GetNBins());
  }
  
  // reserve space for slice histos
  std::vector<T> hists;
  size_t totalNBins = hist_.GetNBins();
  size_t axisNBins = hist_.GetAxes().GetAxis(axisName_).GetNBins();
  hists.reserve(totalNBins/axisNBins);
  
  // get all combinations of bin indices
  size_t dim  = hist_.GetNDims();
  std::vector< std::vector<size_t> > combs;
  combs.reserve(totalNBins/axisNBins);
  std::vector<size_t> vec;
  vec.resize(fixedAxisNames.size());
  
  Combinations::RecursiveCombinations(fixedAxisNames.size(), combs, vec, maxes);
  
  for(std::vector< std::vector<size_t> >::iterator it = combs.begin(); it!=combs.end(); it++){

    // associate axis names with bin indices for the given combination
    std::map<std::string, size_t> sliceMap;
    for(size_t i=0; i<dim-1; i++){
      sliceMap[fixedAxisNames.at(i)]= (*it).at(i);
    }

    // get the slice histo
    T slice = hist_.GetSlice(sliceMap);
    hists.push_back(slice); 
  }
  return hists;
}
