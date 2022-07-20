#include <BinnedEDShrinker.h>
#include <Exceptions.h>
#include <iostream>

BinnedEDShrinker::BinnedEDShrinker(){
    fUsingOverflows = false;
}

void
BinnedEDShrinker::SetBuffer(const std::string& obs_, unsigned lowerBuff_, unsigned upperBuff_){
    fBuffers[obs_] = std::pair<unsigned, unsigned> (lowerBuff_, upperBuff_);
}

std::pair<unsigned, unsigned>
BinnedEDShrinker::GetBuffer(const std::string& obs_) const{
    try{
        return fBuffers.at(obs_);
    }
    catch(const std::out_of_range&){
        throw NotFoundError(Formatter() << "BinnedEDShrinker::Requested buffer boundaries on non-existent observable " << obs_ << "!");
    }
}

std::map<std::string, std::pair<unsigned, unsigned> > 
BinnedEDShrinker::GetBuffers() const{
    return fBuffers;
}

BinAxis
BinnedEDShrinker::ShrinkAxis(const BinAxis& axis_, const unsigned lowerBuff_, 
                              const unsigned upperBuff_) const{
    // no buffer no problem
    if (!lowerBuff_ && !upperBuff_)
        return axis_;

    if ((lowerBuff_  + upperBuff_) >= axis_.GetNBins())
        throw ValueError(Formatter() << "BinnedEDShrinker::Buffersize ("
                                     << lowerBuff_ << ", "  << upperBuff_ << ")" 
                                     << " exceeds number of bins (" << axis_.GetNBins() << ")");
    
    const size_t oldBinCount = axis_.GetNBins();
    const size_t newBinCount = oldBinCount - upperBuff_ - lowerBuff_;

    std::vector<double> lowEdges(newBinCount, 0);
    std::vector<double> highEdges(newBinCount, 0);

    // fill the vectors of the new edges using a simple offset
    for(size_t i = 0; i < newBinCount; i++){
        size_t equivilentOldBin = lowerBuff_ + i;
        lowEdges[i]   = axis_.GetBinLowEdge(equivilentOldBin);
        highEdges[i]  = axis_.GetBinHighEdge(equivilentOldBin);
    }
    
    return BinAxis(axis_.GetName(), lowEdges, highEdges, axis_.GetLatexName());
}

void
BinnedEDShrinker::SetBinMap(const BinnedED& dist_ ) {

  // No buffer no problem. FIXME: what about if all the values are zero?
  if (!fBuffers.size())
    return ;

  size_t nDims = dist_.GetNDims();

  // FIXME Add a check to see if the non zero entries of fBuffers are in the pdf and give warning

  // 1. Build new axes. ShrinkPdf method just makes a copy if buffer size is zero
  
  for(size_t i = 0; i < nDims; i++){
    const std::string& axisName = dist_.GetAxes().GetAxis(i).GetName();
    if (!fBuffers.count(axisName))
      fNewAxes.AddAxis(dist_.GetAxes().GetAxis(i));
    else
      fNewAxes.AddAxis(ShrinkAxis(dist_.GetAxes().GetAxis(i),
				 fBuffers.at(axisName).first,
				 fBuffers.at(axisName).second));
  }

  // 2. Find corresponding bins
  std::vector<size_t> newIndices(dist_.GetNDims());  // same as old, just corrected for overflow
  int   offsetIndex = 0; // note taking difference of two unsigneds
  size_t newBin = 0;     //  will loop over dims and use this to assign bin # corrected for overflow

  const AxisCollection& axes = dist_.GetAxes();

  // bin by bin of old pdf
  for(size_t i = 0; i < dist_.GetNBins(); i++){
    fUseContent.push_back(true);
    // work out the index of this bin in the new shrunk pdf.
    for(size_t j = 0; j < nDims; j++){
      std::string axisName = axes.GetAxis(j).GetName();
      offsetIndex = axes.UnflattenIndex(i, j);            // the index in old pdf
      if (fBuffers.count(axisName))          // offset by lower buffer if nonzero
	offsetIndex -= fBuffers.at(axisName).first;

      // Correct the ones that fall in the buffer regions
      // bins in the lower buffer have negative index. Put in first bin in fit region or ignore
      if (offsetIndex < 0){
	offsetIndex = 0;
	// If not using overflows, flag this bin as not having it's contents included
	if(!fUsingOverflows)
	  fUseContent.at(i) = false;
      }
      // bins in the upper buffer have i > number of bins in axis j. Do the same
      if (offsetIndex >= fNewAxes.GetAxis(j).GetNBins()){
        offsetIndex = fNewAxes.GetAxis(j).GetNBins() - 1;
	// If not using overflows, flag this bin as not having it's contents included
	if(!fUsingOverflows)
	  fUseContent.at(i) = false;
      }

      newIndices[j] = offsetIndex;
    }
    // Fill
    newBin = fNewAxes.FlattenIndices(newIndices);
    fBinVec.push_back(newBin);
  }
} 

BinnedED
BinnedEDShrinker::ShrinkDist(const BinnedED& dist_) const{

    // No buffer no problem. FIXME: what about if all the values are zero?
    if (!fBuffers.size())
        return dist_;

    // Initialise the new pdf with same observables

    BinnedED newDist(dist_.GetName() + "_shrunk", fNewAxes);
    newDist.SetObservables(dist_.GetObservables());
    

    // Fill the axes
    size_t newBin = 0;     //  will loop over dims and use this to assign bin # corrected for overflow
    double content = 0;
    
    // bin by bin of old pdf    
    for(size_t i = 0; i < dist_.GetNBins(); i++){

      content = dist_.GetBinContent(i);
      if(!content) // no content no problem
	continue;

      newBin = fBinVec.at(i);
      // Check if binn has been flagged to not use it's content
      if(fUseContent.at(i))
	newDist.AddBinContent(newBin, content);
      }
    return newDist;
}

void
BinnedEDShrinker::SetUsingOverflows(bool b_){
    fUsingOverflows = b_;
}

bool
BinnedEDShrinker::GetUsingOverflows() const{
    return fUsingOverflows;
}

