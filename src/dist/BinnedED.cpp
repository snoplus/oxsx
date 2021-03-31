#include <BinnedED.h>
#include <Exceptions.h>
#include <Event.h>
#include <Combinations.hpp>
#include <iostream>


BinnedED::BinnedED(const std::string& name_, const AxisCollection& axes_){
    fName = name_;
    fHistogram.SetAxes(axes_);
}

BinnedED::BinnedED(const std::string& name_, const Histogram& histo_){
    fName = name_;
    fHistogram = histo_;
}

void
BinnedED::SetObservables(const std::vector<std::string>& rep_){
    fObservables = ObsSet(rep_);
}

const std::vector<std::string>&
BinnedED::GetObservables() const {
    return fObservables.GetNames();
}

const Histogram&
BinnedED::GetHistogram() const{
    return fHistogram;
}

void
BinnedED::SetHistogram(const Histogram& hist_){
    fHistogram = hist_;
}

void 
BinnedED::Fill(const Event& data_, double weight_){
    try{
        fHistogram.Fill(data_.ToObsSet(fObservables), weight_);
    }
    catch (const DimensionError& e_){
        throw RepresentationError(std::string("Representation in compatible with pdf ") + e_.what());
    }
}

size_t 
BinnedED::FindBin(const Event& data_) const{
    try{
        return fHistogram.FindBin(data_.ToObsSet(fObservables));    
    }
    catch (const DimensionError& e_){
        throw RepresentationError(std::string("Representation in compatible with pdf ") + e_.what());
    }    
}

EventDistribution* 
BinnedED::Clone() const{
    return static_cast<EventDistribution*>(new BinnedED(*this));
}


std::string 
BinnedED::GetName() const{
    return fName;
}

void 
BinnedED::SetName(const std::string& name_){
    fName = name_;
}

//////////////////////////////////////////////////////////////////////////////////////////
// All methods below this line just forward the call to the underlying histogram object //
//////////////////////////////////////////////////////////////////////////////////////////

void 
BinnedED::SetAxes(const AxisCollection& axes_){
    fHistogram.SetAxes(axes_);

}

const AxisCollection& 
BinnedED::GetAxes() const{
    return fHistogram.GetAxes();
}

double 
BinnedED::Probability(const std::vector<double>& vals_) const{
    return fHistogram.GetBinContent(fHistogram.FindBin(vals_))/fHistogram.Integral();
}

double
BinnedED::Probability(const Event& oberservations_) const{
    try{
        return Probability(oberservations_.ToObsSet(fObservables));
    }

    catch(const RepresentationError& e_){
        throw RepresentationError("BinnedED::Probability() failed with  " 
                                  + std::string(e_.what()) 
                                  + " is the rep set correctly?");
    }
}

double 
BinnedED::Integral() const{
    return fHistogram.Integral();
}

void 
BinnedED::Normalise(){
    fHistogram.Normalise();
}

void
BinnedED::Scale(double s_){
    fHistogram.Scale(s_);
}

void 
BinnedED::Fill(const std::vector<double>& vals_, double weight_){
    fHistogram.Fill(vals_, weight_);
}


void 
BinnedED::Fill(double vals_, double weight_){
    fHistogram.Fill(vals_, weight_);
}

size_t 
BinnedED::FindBin(const std::vector<double>& vals_) const{
    return fHistogram.FindBin(vals_);
    
}

double 
BinnedED::GetBinContent(size_t bin_) const{
    return fHistogram.GetBinContent(bin_);
}

void 
BinnedED::AddBinContent(size_t bin_, double content_){
    fHistogram.AddBinContent(bin_, content_);
}

void 
BinnedED::SetBinContent(size_t bin_, double content_){
    fHistogram.SetBinContent(bin_, content_);
}

size_t 
BinnedED::GetNBins() const{
    return fHistogram.GetNBins();
}

unsigned
BinnedED::GetNDims() const{
    return fHistogram.GetNDims();
}

void 
BinnedED::Empty(){
    fHistogram.Empty();
}

size_t 
BinnedED::FlattenIndices(const std::vector<size_t>& indices_) const{
    return fHistogram.FlattenIndices(indices_);
}

std::vector<size_t> 
BinnedED::UnpackIndices(size_t bin_) const{
    return fHistogram.UnpackIndices(bin_);
}

std::vector<double> 
BinnedED::GetBinContents() const{
    return fHistogram.GetBinContents();
}
void 
BinnedED::SetBinContents(const std::vector<double>& data_){
    return fHistogram.SetBinContents(data_);
}

std::vector<double>
BinnedED::Means() const{
    return fHistogram.Means();
}

std::vector<double>
BinnedED::Variances() const{
    return fHistogram.Variances();
}

BinnedED 
BinnedED::Marginalise(const std::vector<std::string>& names_) const{
    // Find the relative indicies indicies in 
    ObsSet newRep = ObsSet(names_);

    // create a name for the projection
    Formatter f;
    f << fName;
    for(size_t i = 0; i < names_.size(); i++)
        f << "_" << names_.at(i);
    f << "_proj";

    // Marginalise the histogram
    BinnedED proj(std::string(f), fHistogram.Marginalise(names_));
    proj.SetObservables(names_);
    return proj;
}

BinnedED
BinnedED::Marginalise(const std::string& name_) const{
    return Marginalise(std::vector<std::string>(1, name_));
}

BinnedED 
BinnedED::GetSlice(const std::map<std::string, size_t>& namesAndBins_) const{
    
    const std::vector<std::string> allAxisNames = fHistogram.GetAxes().GetAxisNames();
    std::vector<std::string> newObsNames;

    // create a name for the slice and get the observable name
    Formatter f;
    f << fName;
    for(std::vector<std::string>::const_iterator it = allAxisNames.begin(); it!=allAxisNames.end(); it++){
        if(namesAndBins_.find(*it) == namesAndBins_.end()){
	    //this is the axis we want the slice in
	    newObsNames.push_back(*it);
	}else{
	    f << "_" << *it << namesAndBins_.at(*it);
	}
    }
    f << "_slice";

    // slice the histogram
    BinnedED slice(std::string(f), fHistogram.GetSlice(namesAndBins_));
    slice.SetObservables(newObsNames);
    return slice;
}

void
BinnedED::Add(const BinnedED& other_, double weight){
    if(other_.fObservables != fObservables)
        std::cout << "Warning::Adding distributions with different observables" << std::endl;
    try{
        fHistogram.Add(other_.fHistogram, weight);
    }
    catch(const ValueError& e_){
        throw ValueError("BinnedED::Add can't add distributions with different binning definitions!");
    }
}


void
BinnedED::Multiply(const BinnedED& other_){
    if(other_.fObservables != fObservables)
        std::cout << "Warning::Multiplying distributions with different observables" << std::endl;
    try{
        fHistogram.Multiply(other_.fHistogram);
    }
    catch(const ValueError& e_){
        throw ValueError("BinnedED::Add can't add distributions with different binning definitions!");
    }
}


void
BinnedED::Divide(const BinnedED& other_){
    if(other_.fObservables != fObservables)
        std::cout << "Warning::Dividing distributions with different observables" << std::endl;
    try{
        fHistogram.Divide(other_.fHistogram);
    }
    catch(const ValueError& e_){
        throw ValueError("BinnedED::Add can't add distributions with different binning definitions!");
    }
}


void
BinnedED::AddPadding(double padding_){
  fHistogram.AddPadding(padding_);
}
