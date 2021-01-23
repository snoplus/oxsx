#include <AxisCollection.h>
#include <iostream>
#include <Exceptions.h>
#include <Formatter.hpp>
#include <ContainerTools.hpp>
#include <algorithm>

size_t 
AxisCollection::GetNBins() const {
    return fNBins;
}
size_t AxisCollection::GetNDimensions() const {
    return fNDimensions;
}

const BinAxis& 
AxisCollection::GetAxis(size_t axisIndex_) const{
    try {
        return fAxes.at(axisIndex_);
    }

    catch (const std::out_of_range&){
        throw NotFoundError(Formatter() << "AxisCollection: asked for axis "
                            << axisIndex_
                            << ", there are only 0 - "
                            << fNDimensions -1
                            );
    }
}

const BinAxis& 
AxisCollection::GetAxis(const std::string& axisName_) const{
    return fAxes.at(GetAxisIndex(axisName_));
}

size_t
AxisCollection::GetAxisIndex(const std::string& name_) const{
    std::vector<std::string>::const_iterator it = std::find(fAxisNames.begin(), fAxisNames.end(), name_);
    if(it == fAxisNames.end())
        throw NotFoundError("No axis by the name of " + name_ + " found!" + "\n axes are: " + ContainerTools::ToString(fAxisNames));
    
    return it - fAxisNames.begin();
}

std::vector<std::string> 
AxisCollection::GetAxisNames() const{
    return fAxisNames;
}

void 
AxisCollection::CountBins(){
    if (!fNDimensions)
        fNBins = 0;

    else{
        fNBins = 1;
        for(size_t i = 0; i < fNDimensions; i++)
            fNBins *= fAxisNbins[i];
    }
}

void 
AxisCollection::AddAxis(const BinAxis& axis_){
    if (HasAxis(axis_.GetName())){
        throw ValueError(Formatter() << "Can't add axis " 
                         << axis_.GetName() <<  " already exists!"); 
    }
    fAxes.push_back(axis_);
    fAxisNames.push_back(axis_.GetName());
    fAxisNbins.push_back(axis_.GetNBins());
    fNDimensions++;
    CountBins();
}

void 
AxisCollection::AddAxes(const std::vector<BinAxis>& axes_){
    for(size_t i = 0; i < axes_.size(); i++)
        AddAxis(axes_[i]);
}

size_t 
AxisCollection::FlattenIndices(const std::vector<size_t>& indices_) const{
    if(indices_.size() != fNDimensions)
        throw DimensionError("AxisCollection::FlattenIndices", 
                             fNDimensions , indices_.size());
                             

    if(fNDimensions == 1)
        return indices_.at(0);

    size_t index = indices_.back();
    for (size_t i = 0; i < fNDimensions-1; i++){
        size_t x = indices_[i];
        for(size_t j = i+1; j < fNDimensions; j++)
            x *= fAxisNbins[j];
        index += x;
    }
    return index;
}


size_t 
AxisCollection::FindBin(const std::vector<double>& vals_) const{
    if (vals_.size() != fNDimensions)
        throw DimensionError("AxisCollection::FindBin", fNDimensions, 
                             vals_.size());

    std::vector<size_t> indices(fNDimensions, 0);
    for (size_t i = 0; i < fNDimensions; i++)
        indices[i] = fAxes.at(i).FindBin(vals_[i]);

    return FlattenIndices(indices);
}

size_t 
AxisCollection::UnflattenIndex(size_t index_, size_t dim_) const{
    if(fNDimensions == 1)
        return index_;

    if (index_ >= fNBins)
        throw NotFoundError(Formatter() 
                            << "Can't unflatten index in non-existent"
                            << "dimension  "
                            << dim_);

    size_t x =1;
    for (size_t i = dim_+1 ; i < fNDimensions; i++ )
        x *= fAxisNbins[i];
    return (index_/x) % fAxisNbins[dim_];
}

std::vector<size_t>
AxisCollection::UnpackIndices(size_t index_) const{
    std::vector<size_t> indices(fNDimensions, 0);
    for(size_t i = 0; i < fNDimensions; i++)
        indices[i] = UnflattenIndex(index_, i);
    return indices;
}

bool 
AxisCollection::HasAxis(const std::string& name_){
    for(size_t i = 0; i < fNDimensions; i++)
        if(name_ == fAxes[i].GetName())
            return true;
    return false;
}

void
AxisCollection::GetBinLowEdges(size_t bin_, std::vector<double>& output_) const{

    if(bin_ > fNBins)
        throw NotFoundError(Formatter() << "GetBinLowEdges on bin "
                            << bin_
                            << " nBins = " << fNBins);

    if(output_.size() != fNDimensions)
        throw DimensionError("GetBinEdges: Input return vector.size() != fNDims ");

    for(size_t i = 0; i <fNDimensions; i++){
        output_[i] = fAxes[i].GetBinLowEdge(UnflattenIndex(bin_, i));
    }
}

void
AxisCollection::GetBinHighEdges(size_t bin_, std::vector<double>& output_) const{

    if(bin_ > fNBins)
        throw NotFoundError(Formatter() << "GetBinHighEdges on bin "
                            << bin_
                            << " nBins = " << fNBins);
    
    if(output_.size() != fNDimensions)
        throw DimensionError("GetBinEdges: Input return vector.size() != fNDims ");

    for(size_t i = 0; i <fNDimensions; i++){
        output_[i] = fAxes[i].GetBinHighEdge(UnflattenIndex(bin_, i));
    }
}

void
AxisCollection::GetBinCentres(size_t bin_, std::vector<double>& output_) const{
    if(bin_ > fNBins)
        throw NotFoundError(Formatter() << "GetBinHighEdges on bin "
                            << bin_
                            << " nBins = " << fNBins);
    
    if(output_.size() != fNDimensions)
        throw DimensionError("GetBinEdges: Input return vector.size() != fNDims ");
    

    for(size_t i = 0; i <fNDimensions; i++){
        output_[i] = fAxes[i].GetBinCentre(UnflattenIndex(bin_, i));
    }
}

double
AxisCollection::GetBinLowEdge(size_t bin_, size_t dim_) const{
    if(dim_ >= fNDimensions)
        throw NotFoundError(Formatter() << "Axis Collection::Bin edge on non existent dimension " << dim_  << "!" );

    if(bin_ >= fNBins)
        throw NotFoundError(Formatter() << "AxisCollection::Called bin edge on out of bounds bin " << bin_ << "!");
    
    return fAxes.at(dim_).GetBinLowEdge(UnflattenIndex(bin_, dim_));
}


double
AxisCollection::GetBinHighEdge(size_t bin_, size_t dim_) const{
    if(dim_ >= fNDimensions)
        throw NotFoundError(Formatter() << "Axis Collection::Bin edge on non existent dimension" << dim_ <<  "!");
   
    if(bin_ >= fNBins)
        throw NotFoundError(Formatter() << "AxisCollection::Called bin edge on out of bounds bin " << bin_ << "!");
    
    return fAxes.at(dim_).GetBinHighEdge(UnflattenIndex(bin_, dim_));
}

double
AxisCollection::GetBinCentre(size_t bin_, size_t dim_) const{
    if(dim_ >= fNDimensions)
        throw NotFoundError(Formatter() << "Axis Collection::Bin center on non existent dimension" << dim_ <<  "!");
   
    if(bin_ >= fNBins)
        throw NotFoundError(Formatter() << "AxisCollection::Called bin center on out of bounds bin " << bin_ << "!");
    
    return fAxes.at(dim_).GetBinCentre(UnflattenIndex(bin_, dim_));
}

double
AxisCollection::GetBinWidth(size_t bin_, size_t dim_) const{
    if(dim_ >= fNDimensions)
        throw NotFoundError(Formatter() << "Axis Collection::Bin width on non existent dimension" << dim_ <<  "!");
   
    if(bin_ >= fNBins)
        throw NotFoundError(Formatter() 
                            << "AxisCollection::Called bin width on out of bounds bin " << bin_ << "!");
    
    return fAxes.at(dim_).GetBinWidth(UnflattenIndex(bin_, dim_));
}

bool
AxisCollection::operator==(const AxisCollection& other_) const{
    return fAxes == other_.fAxes;
}

bool
AxisCollection::operator!=(const AxisCollection& other_) const{
    return fAxes != other_.fAxes;
}
