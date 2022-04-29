#include <Systematic.h>
#include <Exceptions.h>
#include <algorithm>
#include <iostream>

BinnedED 
Systematic::operator() (const BinnedED& pdf_, double* norm) const{
    try{
        BinnedED afterSmear = pdf_;
        afterSmear.SetBinContents(fResponse(pdf_.GetBinContents()));
        if (norm != nullptr) {
            *norm = afterSmear.Integral();
            afterSmear.Scale(1./(*norm));
        } else {
            afterSmear.Normalise();
        }
        return afterSmear;
    }
    catch(const DimensionError& e_){
        throw DimensionError(std::string("Failed to apply systematic :: ") + e_.what());
    }
}

void 
Systematic::SetResponse(const SparseMatrix& responseMatrix_){
    fResponse = responseMatrix_;
}

const SparseMatrix& 
Systematic::GetResponse() const{
    return fResponse;
}

void
Systematic::SetTransformationObs(const ObsSet& rep_) {
    fTransObs = rep_;
}

ObsSet
Systematic::GetTransformationObs() const {
    return fTransObs;
}

void
Systematic::SetDistributionObs(const ObsSet& rep_) {
    fDistObs = rep_;
}

ObsSet
Systematic::GetDistributionObs() const {
    return fDistObs;
}

bool
Systematic::BinsCompatible(size_t bin1_, size_t bin2_) const{
    std::vector<size_t> bin1Indices = fAxes.UnpackIndices(bin1_);
    std::vector<size_t> bin2Indices = fAxes.UnpackIndices(bin2_);

    // Where are the indices the systematic cares about in the pdfs index scheme
    std::vector<size_t> relativeIndices = fTransObs.GetRelativeIndices(fDistObs);

    // Do the two global bin numbers have the same indices except for in the dimisensions 
    // this systematic affects?
    for(size_t i = 0; i < bin1Indices.size(); i++){
        if (VectorContains(relativeIndices, i))
            continue;
        if (bin1Indices.at(i) != bin2Indices.at(i))
            return false;
    }
        
    return true;
}


bool
Systematic::VectorContains(const std::vector<size_t>& vec_,  size_t val_) const{
    return std::find(vec_.begin(), vec_.end(), val_) != vec_.end();
}


const AxisCollection&
Systematic::GetAxes() const{
    return fAxes;
}

void
Systematic::SetAxes(const AxisCollection& axes_){
    fAxes = axes_;
    SetResponse(SparseMatrix(axes_.GetNBins(),axes_.GetNBins()));

}
