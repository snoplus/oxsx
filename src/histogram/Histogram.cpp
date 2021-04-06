#include <Histogram.h>
#include <Exceptions.h>
#include <Formatter.hpp>
#include <Combinations.hpp>
#include <ContainerTools.hpp>
#include <iostream>
#include <set>
#include <algorithm>

Histogram::Histogram(const AxisCollection& axes_){
    SetAxes(axes_);
}

void 
Histogram::SetAxes(const AxisCollection& axes_){
    fAxes  = axes_;
    fNBins = fAxes.GetNBins();
    fNDims = fAxes.GetNDimensions();
    fBinContents.resize(fNBins, 0);
    
}

const AxisCollection& 
Histogram::GetAxes() const{
    return fAxes;
}

double 
Histogram::Integral() const{
    double sum = 0;
    for(size_t i = 0; i < fNBins; i++)
        sum += fBinContents[i];
    return sum;
}

void 
Histogram::Normalise(){
    double sum = Integral();
    for(size_t i = 0; i < fNBins; i++)
        fBinContents[i] /= sum;
}

void
Histogram::Scale(double s_){
    for(size_t i = 0; i < fNBins; i++)
        fBinContents[i] *= s_;
}

void 
Histogram::Fill(const std::vector<double>& vals_, double weight_){
    if(vals_.size() != fNDims)                             
        throw DimensionError("Histogram::Fill", fNDims, vals_.size());

    fBinContents[FindBin(vals_)] += weight_;
}

void 
Histogram::Fill(const std::map<std::string, double>& vals_, double weight_){
    try{
        Fill(ContainerTools::GetValues(vals_, GetAxisNames()), weight_);
    }
    catch(const std::out_of_range& e_){
        throw NotFoundError("Tried to fill a histogram with incomplete dictionary!");
    }
}


void 
Histogram::Fill(double vals_, double weight_){
    Fill(std::vector<double>(1, vals_), weight_);
}

size_t 
Histogram::FindBin(const std::vector<double>& vals_) const{
    return fAxes.FindBin(vals_);
    
}

double 
Histogram::GetBinContent(size_t bin_) const{
    if(bin_ > fNBins)
        throw NotFoundError(Formatter() 
                             << "Out of bounds bin access attempted on bin "
                             << bin_ <<  " !");
    return fBinContents[bin_];
}

void 
Histogram::AddBinContent(size_t bin_, double content_){
    if(bin_ > fNBins)
        throw NotFoundError(Formatter() 
                             << "Out of bounds bin access attempted on bin "
                             << bin_ <<  " !");
    fBinContents[bin_] += content_;
}

void 
Histogram::SetBinContent(size_t bin_, double content_){
    if(bin_ > fNBins)
        throw NotFoundError(Formatter()  << "Out of bounds bin access attempted on bin " << bin_ <<  " !");
    fBinContents[bin_] = content_;
}

size_t 
Histogram::GetNBins() const{
    return fNBins;
}

size_t 
Histogram::GetNDims() const{
  return fNDims;
}

void 
Histogram::Empty(){
    for(size_t i = 0; i < fNBins; i++)
        fBinContents[i] = 0;
}

size_t 
Histogram::FlattenIndices(const std::vector<size_t>& indices_) const{
    return fAxes.FlattenIndices(indices_);
}

std::vector<size_t> 
Histogram::UnpackIndices(size_t bin_) const{
    return fAxes.UnpackIndices(bin_);
}

std::vector<double> 
Histogram::GetBinContents() const{
    return fBinContents;
}
void 
Histogram::SetBinContents(const std::vector<double>& data_){
    if (data_.size() != fNBins)
        throw DimensionError("Histogram::SetBinContents", fNBins, 
                             data_.size());

    fBinContents = data_;
    fNBins = fBinContents.size();
}

std::vector<double>
Histogram::Means() const{
    std::vector<double> means(fNDims, 0);    
    for(size_t i = 0; i < fNBins; i++)
        for(size_t j = 0; j < fNDims; j++)
            means[j] += fBinContents.at(i) * fAxes.GetBinCentre(i, j);
    return means;
}

std::vector<double>
Histogram::Variances() const{
    std::vector<double> variances(fNDims, 0);

    for(size_t i = 0; i < fNBins; i++)
        for(size_t j = 0; j < fNDims; j++){
            double binCent = fAxes.GetBinCentre(i, j);
            variances[j] += binCent * binCent *  fBinContents.at(i);
        }
    
    

    std::vector<double> means = Means();
    for(size_t i = 0; i < fNDims; i++)
        variances[i] -= means.at(i) * means.at(i);

    return variances;
}

Histogram
Histogram::Marginalise(const std::vector<std::string>& axes_) const{
    std::vector<std::string> allAxisNames = fAxes.GetAxisNames();

    // check the pdf does contain the axes asked for
    for(size_t i = 0; i < axes_.size(); i++){

        if (std::find(allAxisNames.begin(), allAxisNames.end(), axes_.at(i)) == allAxisNames.end())
            throw NotFoundError(Formatter()
                                << "Histogram::Marginalise::Tried "
                                << "to project out non existent axis "
                                << axes_.at(i) << "!");
    }


    // work out which axis number corresponds to each name
    std::vector<size_t> indices;
    for(size_t i = 0; i < axes_.size(); i++)
        indices.push_back(fAxes.GetAxisIndex(axes_.at(i)));
        
    // Get the axes you are interested in, in the order requested
    AxisCollection newAxes;
    for(size_t i = 0;  i < axes_.size(); i++)
        newAxes.AddAxis(fAxes.GetAxis(axes_.at(i)));

    // New histogram
    Histogram marginalised(newAxes);

    std::vector<size_t> oldIndices(fNDims);
    std::vector<size_t> newIndices(axes_.size());
    size_t newBin = -1;

    // Now loop over the bins in old and fill new pdfs 
    for(size_t bin = 0; bin < fNBins; bin++){
        for(size_t i = 0; i < fNDims; i++)
            oldIndices[i] = fAxes.UnflattenIndex(bin, i);

        for(size_t i = 0; i < indices.size(); i++)
            newIndices[i] = oldIndices.at(indices.at(i));

        newBin = marginalised.FlattenIndices(newIndices);
        marginalised.AddBinContent(newBin, fBinContents.at(bin));
    }
    return marginalised;
}

Histogram
Histogram::Marginalise(const std::string& index_) const {
    return Marginalise(std::vector<std::string>(1, index_));
}

void 
Histogram::Recurse(size_t numFreeIdx, std::vector<size_t> binsEachAxis, std::vector<size_t> coords, std::vector<std::vector<size_t> > &localIdx, std::vector<size_t> free_ax) const{
    // function used in GetSlice method to generate every possible combination 
    // of local indicies using on numFreeIdx recursively nested for-loops.

    // check if at deepest level 
    std::cout << "Inside recursive function! numFreeIdx is " << numFreeIdx << std::endl;
    if(numFreeIdx < 1){
        // assign the current combination of indexes to vector
        std::vector<size_t> current_idx; 
        for(size_t idx = 0; idx < coords.size(); idx++){
            //std::cout << "Assigning current coords to a vector" << std::endl;
            current_idx.push_back(coords[idx]); 
        }

        // push this current coordinate vector to be saved with all the 
        // others 
        //std::cout << "Pushing current_idx to the giga list localIDx" << std::endl;
        localIdx.push_back(current_idx); 
        std::cout << "LocalIdx size (inside) " << localIdx.size() << std::endl;
        //std::cout << current_idx[0] << current_idx[1] << current_idx[2] << current_idx[3] << std::endl;
    } 
    else{
        // need to go deeper! Begin this level's loop
        //std::cout << "Going deeper!" << std::endl;
        for(size_t i = 0; i < binsEachAxis.at(numFreeIdx-1); i++){
            // fill relevant idx 
            //std::cout << "Inside a loop, filling coords element " << free_ax[numFreeIdx-1] << " with " << i << std::endl;
            coords[free_ax[numFreeIdx-1]] = i; 
            //std::cout << "FILLED" << std::endl; 
            // go DEEPER 
            Recurse(numFreeIdx-1, binsEachAxis, coords, localIdx, free_ax);
        }
    }
}

Histogram
Histogram::GetSlice(const std::map<std::string,size_t>& fixedBins_) const{
    
    // all the axis names in the initial (pre sliced) histogram 
    const std::vector<std::string> allAxisNames = fAxes.GetAxisNames();
    
    // check the pdf does contain the axes and bins asked for
    for(std::map<std::string, size_t>::const_iterator it = fixedBins_.begin(); it!=fixedBins_.end(); it++){
        if (std::find(allAxisNames.begin(), allAxisNames.end(), it->first) == allAxisNames.end())
            throw NotFoundError(Formatter()
			      << "Histogram::GetSlice::Tried "
			      << "to get slice of non existent axis "
			      << it->first << "!");
    }

    // want to support multidimensional slices 
    if (fixedBins_.size() >= fNDims)
        throw DimensionError("Histogram::GetSlice", fNDims -1, fixedBins_.size());

    // create vectors to hold the slicing idxs 
    std::vector<std::string> newAxisNames; 
    std::vector<size_t> newAxisIndexes;
    std::vector<size_t> sliceIndices; 
    std::vector<bool> all_bins; // is this axis a fixed bin?
    for(std::vector<std::string>::const_iterator it = allAxisNames.begin(); it!=allAxisNames.end(); it++){
        
        // goes through entire fixedBins_ map and reaches end without finding axis in fixedBins 
        if(fixedBins_.find(*it) == fixedBins_.end()){
	        // this is the axis we want to slice out! 
            newAxisNames.push_back(*it); 
            newAxisIndexes.push_back(fAxes.GetAxisIndex(*it));
            std::cout << "Free Index: " << *it << std::endl; 
            all_bins.push_back(false);
	    } else{
	    // record the indices for the fixed bins in other dimensions
	    //sliceIndices[fAxes.GetAxisIndex(*it)] = fixedBins_.at(*it);
        sliceIndices.push_back(fAxes.GetAxisIndex(*it)); 
        all_bins.push_back(true); 
        std::cout << "Fixed Indices recorded in Axis " << *it << " with idx " << fixedBins_.at(*it) << std::endl; 
	    }
    }

    // new histogram loop 
    AxisCollection newAxes; 
    std::vector<size_t> binsEachAxis;  

    // vector of vectors to store each local idx
    std::vector<std::vector<size_t> > localIdx; 

    // vector to store current local idx 
    std::vector<size_t> coords(fNDims); 

    // creating the slice  
    for(size_t i = 0; i < newAxisNames.size(); i++){
        newAxes.AddAxis(fAxes.GetAxis(newAxisIndexes[i]));
        
        // get num bins in each axis
        BinAxis individualAxis = newAxes.GetAxis(newAxisNames[i]);
        size_t bins = individualAxis.GetNBins();
        binsEachAxis.push_back(bins);
    }
    Histogram slice(newAxes);
    
    // fill in the fixed bins 
    std::cout << "all bins " << all_bins[0] << " " << all_bins[1] << " " << all_bins[2] << " " << all_bins[3] << std::endl;
    for(size_t i = 0; i < fNDims; i++){
        if(all_bins.at(i) == true){
            std::cout << " Filling fixed bin idx " << i << " with " << fixedBins_.at(allAxisNames[i]) << std::endl;
            coords[i] = fixedBins_.at(allAxisNames[i]); 
        }
    }
    std::cout << "coords before: " << coords[0] << " " << coords[1] << " " << coords[2] << coords[3] << std::endl;  
    // recursively fill localIdxs vector with N nested for loops, where N is number of free axis 
    std::cout << "About to call recursive function!" << std::endl;
    Recurse(newAxisNames.size(), binsEachAxis, coords, localIdx, newAxisIndexes);
    std::cout << "LocalIdx size (outside) " << localIdx.size() << std::endl;
    for(size_t i = 0; i < localIdx.size(); i++){
        std::vector<size_t> val = localIdx[i];
        std::cout << "saved value: " << val[0] << " " << val[1] << " " << val[2] << " " << val[3] << std::endl;
          
    }
    // loop over the (now filled?) localIdxs
    size_t bin = 0;   
    for(size_t idx = 0; idx < localIdx.size(); idx++){
        // convert a given localIdx to globalIdx
        size_t oldGlobalBin = fAxes.FlattenIndices(localIdx.at(idx)); 
        slice.AddBinContent(bin, fBinContents.at(oldGlobalBin)); 
        bin += 1; 
    } 
    // // array holds the fNDims dimensional coordinates of every bin in slice in terms 
    // // of original histogram   
    // //size_t localIdx[slice.GetNBins()][allAxisNames.size()];
    // std::vector<size_t> localIdx(fNDims); 
    // // loop over each bin in slice
    // //for(size_t bin = 0; bin < slice.GetNBins(); bin++){ 
    //     // fill in the fixed indexes 
        
    //     for(size_t axFixIdx = 0; axFixIdx < sliceIndices.size(); axFixIdx++){
    //         localIdx[sliceIndices[axFixIdx]] = fixedBins_.at(allAxisNames[sliceIndices[axFixIdx]]);   
    //     }
    //     // loop over the free idx in each dimension 
    //     for(size_t freeBin = 0; freeBin < binsEachAxis[i]; freeBin++){
    //         std::cout << "New Slice Bin! " << std::endl; 
                
    //         // set the free idx 
    //         localIdx[newAxisIndexes[i]] = freeBin; 
            
    //         std::cout<< "localIdx = "; 
    //         for(size_t j = 0; j < localIdx.size(); j++){
    //             std::cout << localIdx[j] << " "; 
    //         }
    //         std::cout << "\n"; 
                
    //         // compute the global idx coodinates based on these local ones in 
    //         // original histogram  
    //         size_t oldGlobalBin = fAxes.FlattenIndices(localIdx);
	//         slice.AddBinContent(bin, fBinContents.at(oldGlobalBin));
    //         std::cout << "Added " << fBinContents.at(oldGlobalBin) << " to bin " << bin << std::endl; 
            
    //     }
    //}

    return slice;
}

double
Histogram::GetBinLowEdge(size_t bin_, size_t index_) const{
    return fAxes.GetBinLowEdge(bin_, index_);
}

double
Histogram::GetBinHighEdge(size_t bin_, size_t index_) const{
    return fAxes.GetBinHighEdge(bin_, index_);
}

double
Histogram::GetBinCentre(size_t bin_, size_t index_) const{
    return fAxes.GetBinCentre(bin_, index_);
}

void
Histogram::Add(const Histogram& other_, double weight_){
    if(other_.GetAxes() != GetAxes())
        throw ValueError(Formatter() << "Histogram::Add can't add histograms with different binning definitions!");
    
    for(size_t i = 0; i < GetNBins(); i++)
        AddBinContent(i, other_.GetBinContent(i) * weight_);
}


void
Histogram::Multiply(const Histogram& other_){
    if(other_.GetAxes() != GetAxes())
        throw ValueError(Formatter() << "Histogram::Add can't add histograms with different binning definitions!");
    
    for(size_t i = 0; i < GetNBins(); i++)
        SetBinContent(i, GetBinContent(i) * other_.GetBinContent(i));
}


void
Histogram::Divide(const Histogram& other_){
    if(other_.GetAxes() != GetAxes())
        throw ValueError(Formatter() << "Histogram::Add can't add histograms with different binning definitions!");
    
    for(size_t i = 0; i < GetNBins(); i++)
        SetBinContent(i, GetBinContent(i) / other_.GetBinContent(i));
}

std::vector<std::string>
Histogram::GetAxisNames() const{
    return fAxes.GetAxisNames();
}


size_t
Histogram::GetAxisIndex(const std::string& name_) const{
    return fAxes.GetAxisIndex(name_);
}
void
Histogram::AddPadding(double padding_){
  std::vector<double> newBinContents;
  for(size_t i =0; i<fBinContents.size();i++){
    if(fBinContents.at(i)==0)
      newBinContents.push_back(padding_);
    else
      newBinContents.push_back(fBinContents[i]);
  }  
  fBinContents = newBinContents;
}
