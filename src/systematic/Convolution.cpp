#include <Convolution.h>
#include <PDF.h>
#include <JumpPDF.h>
#include <DenseMatrix.h>
#include <Exceptions.h>
#include <string>
void 
Convolution::SetFunction(PDF* function_){
  // wrap this up if position independent kernel of the form P(x | x2) = P(x - x2)
  delete fDist;
  fDist = static_cast<ConditionalPDF*>(new JumpPDF("kernel", function_));
}

void
Convolution::SetConditionalPDF(ConditionalPDF* c_){
  delete fDist;
  fDist = c_->Clone();
}

Convolution::~Convolution(){
  delete fDist;
}

void 
Convolution::Construct(){
  if (!fDist || !fAxes.GetNBins())
    throw LogicError("Convolution::Construct() : Tried to construct convolution without axes or function/distribution, or both!!");
    
  if(!fCachedCompatibleBins)
    CacheCompatibleBins();

  size_t nBins = fAxes.GetNBins();
  size_t nDims = fAxes.GetNDimensions();
  const AxisCollection& axes = fAxes;
  std::vector<size_t> relativeIndices = fTransObs.GetRelativeIndices(fDistObs);

  // Work out the transition probabilitites within this sub set of the bins
  std::vector<double> binCentres(fSubMapAxes.GetNDimensions());
  std::vector<double> lowEdges(fSubMapAxes.GetNDimensions());
  std::vector<double> highEdges(fSubMapAxes.GetNDimensions());

  DenseMatrix subMap(fSubMapAxes.GetNBins(), fSubMapAxes.GetNBins());

  for (size_t origBin = 0; origBin < fSubMapAxes.GetNBins(); origBin++){
  // get the centre of the bin. Need to offset by this for a convolution
  fSubMapAxes.GetBinCentres(origBin, binCentres);

  // loop over the bins it can be smeared into 
  for(size_t destBin = 0; destBin < fSubMapAxes.GetNBins(); destBin++){
  fSubMapAxes.GetBinLowEdges(destBin, lowEdges);
  fSubMapAxes.GetBinHighEdges(destBin, highEdges);
            
  subMap.SetComponent(destBin, origBin, fDist -> Integral(lowEdges, highEdges, binCentres));
}        
  }

// Now expand to the full size matrix. Elements are zero by default
// compatible bins are cached, values must match the smaller matrix above
size_t destBin = -1;
std::vector<unsigned> nonZeroRowIndices;
std::vector<unsigned> nonZeroColIndices;
std::vector<double> values;
nonZeroRowIndices.reserve(fCompatibleBins.at(0).size());
nonZeroColIndices.reserve(fCompatibleBins.at(0).size());
    

for(size_t origBin = 0; origBin < axes.GetNBins(); origBin++){
for(size_t i = 0; i < fCompatibleBins.at(origBin).size(); i++){
destBin = fCompatibleBins.at(origBin).at(i);
nonZeroRowIndices.push_back(origBin);
nonZeroColIndices.push_back(destBin);
values.push_back( subMap.GetComponent(fSysBins.at(origBin),
                                        fSysBins.at(destBin)));
}
}
        
fResponse.SetComponents(nonZeroRowIndices, nonZeroColIndices, values);
}


void
Convolution::CacheCompatibleBins(){
  fCompatibleBins.resize(fAxes.GetNBins());
  // only need to look at one side of the matrix, its symmetric
  for(size_t i = 0; i < fAxes.GetNBins(); i++){
    fCompatibleBins.at(i).push_back(i); // always true
    for(size_t j = i+1;  j < fAxes.GetNBins(); j++){
      if(BinsCompatible(i , j)){
        fCompatibleBins.at(i).push_back(j);
        fCompatibleBins.at(j).push_back(i);
      }
    }
  }

  std::vector<size_t> relativeIndices = fTransObs.GetRelativeIndices(fDistObs);
  const AxisCollection& axes = fAxes;

  //  get the axes that this systematic will act on
  fSubMapAxes = AxisCollection();
  for(size_t i = 0; i < relativeIndices.size(); i++)
    fSubMapAxes.AddAxis(axes.GetAxis(relativeIndices.at(i)));
    
  // cache the equivilent index in the binning system of the systematic
  fSysBins.resize(fAxes.GetNBins());
  std::vector<size_t> sysIndices(relativeIndices.size(), 0);
  for(size_t i = 0; i < axes.GetNBins(); i++){
  for(size_t dim = 0; dim < relativeIndices.size(); dim++)
    sysIndices[dim] = axes.UnflattenIndex(i, relativeIndices.at(dim));

  fSysBins[i] = fSubMapAxes.FlattenIndices(sysIndices);
}
  fCachedCompatibleBins = true;
}

  ///////////////////////////////
  // Make this object fittable //
  ///////////////////////////////

  // Fitting this dist to data means adjusting the underlying function

void
  Convolution::RenameParameter(const std::string& old_, const std::string& new_){
  fDist->RenameParameter(old_, new_);
}

void
  Convolution::SetParameter(const std::string& name_, double value_){
  fDist->SetParameter(name_, value_);
}

double
  Convolution::GetParameter(const std::string& name_) const{
  return fDist->GetParameter(name_);
}

void
  Convolution::SetParameters(const ParameterDict& ps_){
  try{
  fDist->SetParameters(ps_);
}
  catch(const ParameterError& e_){
  throw ParameterError("Convolution internal function: " + std::string(e_.what()));
}
}

ParameterDict
  Convolution::GetParameters() const{
  return fDist->GetParameters();
}

size_t
  Convolution::GetParameterCount() const{
  return fDist->GetParameterCount();
}

 std::set<std::string>
   Convolution::GetParameterNames() const{
  return fDist->GetParameterNames();
}

 std::string
   Convolution::GetName() const{
  return fName;
}
void
  Convolution::SetName(const std::string& n_){
  fName = n_;
}
