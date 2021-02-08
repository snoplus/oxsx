#include <BinnedEDShrinker.h>
#include <Exceptions.h>
#include <iostream>

BinnedEDShrinker::BinnedEDShrinker(){
    fUsingOverflows = false;
}

void
BinnedEDShrinker::SetBuffer(size_t dim_, unsigned lowerBuff_, unsigned upperBuff_){
    fBuffers[dim_] = std::pair<unsigned, unsigned> (lowerBuff_, upperBuff_);
}

std::pair<unsigned, unsigned>
BinnedEDShrinker::GetBuffer(size_t dim_) const{
    try{
        return fBuffers.at(dim_);
    }
    catch(const std::out_of_range&){
        throw NotFoundError(Formatter() << "BinnedEDShrinker::Requested buffer boundaries on non-existent dim " << dim_ << "!");
    }
}

std::map<size_t, std::pair<unsigned, unsigned> > 
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

BinnedED
BinnedEDShrinker::ShrinkDist(const BinnedED& dist_) const{

    // No buffer no problem. FIXME: what about if all the values are zero?
    if (!fBuffers.size())
        return dist_;
    
    size_t nDims = dist_.GetNDims();

    // FIXME Add a check to see if the non zero entries of fBuffers are in the pdf and give warning

    // 1. Build new axes. ShrinkPdf method just makes a copy if buffer size is zero
    AxisCollection newAxes;
    const std::vector<size_t> distDataIndices = dist_.GetObservables().GetIndices();
    size_t dataIndex = 0;
    
    for(size_t i = 0; i < nDims; i++){
        dataIndex = distDataIndices.at(i);
        if (!fBuffers.count(dataIndex))
            newAxes.AddAxis(dist_.GetAxes().GetAxis(i));
        else
            newAxes.AddAxis(ShrinkAxis(dist_.GetAxes().GetAxis(i), 
                                       fBuffers.at(dataIndex).first,  
                                       fBuffers.at(dataIndex).second));
    }

    // 2. Initialise the new pdf with same data rep
    BinnedED newDist(dist_.GetName() + "_shrunk", newAxes);
    newDist.SetObservables(dist_.GetObservables());

    // 3. Fill the axes
    std::vector<size_t> newIndices(dist_.GetNDims());  // same as old, just corrected for overflow
    int   offsetIndex = 0; // note taking difference of two unsigneds
    size_t newBin = 0;     //  will loop over dims and use this to assign bin # corrected for overflow

    const AxisCollection& axes = dist_.GetAxes();
    double content = 0;
    // bin by bin of old pdf
    for(size_t i = 0; i < dist_.GetNBins(); i++){
        content = dist_.GetBinContent(i);
        if(!content) // no content no problem
            continue;

        // work out the index of this bin in the new shrunk pdf. 
        for(size_t j = 0; j < nDims; j++){
            offsetIndex = axes.UnflattenIndex(i, j);            // the index in old pdf
            if (fBuffers.count(distDataIndices.at(j)))          // offset by lower buffer if nonzero
                offsetIndex -= fBuffers.at(distDataIndices.at(j)).first;

            // Correct the ones that fall in the buffer regions
            // bins in the lower buffer have negative index. Put in first bin in fit region or ignore
            if (offsetIndex < 0){
                offsetIndex = 0;
                if(!fUsingOverflows)
                    content = 0;

            }

            // bins in the upper buffer have i > number of bins in axis j. Do the same
            if (offsetIndex >= static_cast<int>(newAxes.GetAxis(j).GetNBins())){
                offsetIndex = newAxes.GetAxis(j).GetNBins() - 1;

                if (!fUsingOverflows)
                    content = 0;
            }

            newIndices[j] = offsetIndex;
        }
        // Fill 
        newBin = newAxes.FlattenIndices(newIndices);
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

