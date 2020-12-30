#include <ScaleNonLinear.h>
#include <sstream>
#include <DoubleParameter.h>
#include <Exceptions.h>
#include <ContainerTools.hpp>

double 
ScaleNonLinear::NonLinearScaleFunction(double val_){
    return fScaleFactor1*val_ + fScaleFactor2;
}

void 
ScaleNonLinear::Construct(){
    if (fScaleFactor2 <= 0)
        throw ValueError("Scale factor must be >0 !");
    
    if(fTransObs.GetNObservables() != 1)
        throw RepresentationError("Scale systematic must have a 1D representation!");

    const AxisCollection& axes       = fAxes;
    // the axis to scale
    const size_t  scaleAxisDataIndex = fTransObs.GetIndex(0);
    const BinAxis& scaleAxis         = axes.GetAxis(fDistObs.GetDataIndexPos(scaleAxisDataIndex));


    const size_t nBins               = axes.GetNBins(); 
    const size_t scaleAxisNBins      = scaleAxis.GetNBins(); 
 
    for(size_t i = 0; i < nBins; i++){
        // For each old bin, work out the contributions into all of the new bins
        // indices in other components should be unaffected
        std::vector<size_t> oldIndices = axes.UnpackIndices(i);
        size_t scaleBin                = oldIndices.at(fDistObs.GetDataIndexPos(scaleAxisDataIndex));
        
        double binCentre   = scaleAxis.GetBinLowEdge(scaleBin);
        double scaleFactor = NonLinearScaleFunction(binCentre);

        double scaledLow   = scaleAxis.GetBinLowEdge(scaleBin)  * scaleFactor;
        double scaledHigh  = scaleAxis.GetBinHighEdge(scaleBin) * scaleFactor;
        double scaledWidth = scaledHigh - scaledLow;

        // new bin to map into, mapping only happens if the indices are the same except the one to 
        // scale so, loop over the bins in the scale axes and leave other indices the same
        // the others are zero from initialisation 
        
        std::vector<size_t> newIndices = oldIndices;
        for(size_t j = 0; j < scaleAxisNBins; j++){
            newIndices[fDistObs.GetDataIndexPos(scaleAxisDataIndex)] = j;
            size_t newScaleBin = j;
                        
            double newLow  = scaleAxis.GetBinLowEdge(newScaleBin);
            double newHigh = scaleAxis.GetBinHighEdge(newScaleBin);

            double contribution;
            // Is it in the scale region at all?
            if (newLow > scaledHigh || newHigh < scaledLow) 
                contribution = 0;

            else{
                // Is it fully in the region?
                bool includedFromBelow = newLow > scaledLow;
                bool includedFromAbove = newHigh < scaledHigh;

                if (includedFromBelow && includedFromAbove) 
                    // fully inside
                    contribution = 1/scaleFactor;//GetScaleFactor();

                else if (includedFromBelow) 
                    // spills partly over the top
                    contribution = (scaledHigh - newLow)/scaledWidth;

                else 
                    // spills partly over the bottom
                    contribution = (newHigh - scaledLow)/scaledWidth;
                
            }
            
            fResponse.SetComponent(axes.FlattenIndices(newIndices), i, contribution);
        }
               
    }
    return;
}

void
ScaleNonLinear::SetScaleFactor1(const double scaleFactor1_){
    fScaleFactor1 = scaleFactor1_;
}

double
ScaleNonLinear::GetScaleFactor1() const{
    return fScaleFactor1;
}

void
ScaleNonLinear::SetScaleFactor2(const double scaleFactor2_){
    fScaleFactor2 = scaleFactor2_;
}

double
ScaleNonLinear::GetScaleFactor2() const{
    return fScaleFactor2;
}

//////////////////////////////////////////////////////////////////
// Make this object fittable, so the scale factor is adjustable //
//////////////////////////////////////////////////////////////////

void
ScaleNonLinear::SetParameter(const std::string& name_, double value){
    if(name_ == fParamName1)
        fScaleFactor1 = value;
    else if(name_ == fParamName2)
        fScaleFactor2 = value;
    else
        throw ParameterError("Scale: can't set " + name_ + ", " + fParamName1 + " and " + fParamName2  + " are the only parameters" );
}

double
ScaleNonLinear::GetParameter(const std::string& name_) const{
    if(name_ == fParamName1)
        return fScaleFactor1;
    else if(name_ == fParamName2)
        return fScaleFactor2;
    else
        throw ParameterError("Scale: can't set " + name_ + ", " + fParamName1 + " and " + fParamName2  + " are the only parameters" );
}

void
ScaleNonLinear::SetParameters(const ParameterDict& pd_){
    for (ParameterDict::const_iterator i = pd_.begin(); i != pd_.end(); ++i)
      SetParameter(i->first,i->second);
}

ParameterDict
ScaleNonLinear::GetParameters() const{
    ParameterDict d;
    d[fParamName1] = fScaleFactor1;
    d[fParamName2] = fScaleFactor2;
    return d;
}

size_t
ScaleNonLinear::GetParameterCount() const{
    return 2;
}

std::set<std::string>
ScaleNonLinear::GetParameterNames() const{
    std::set<std::string> set;
    set.insert(fParamName1);
    set.insert(fParamName2);
    return set;
}

void
ScaleNonLinear::RenameParameter(const std::string& old_, const std::string& new_){
    if(old_ == fParamName1)
        fParamName1 = new_;
    else if(old_ == fParamName2)
        fParamName2 = new_;
    else
        throw ParameterError("Scale: can't set " + new_ + ", " + fParamName1 + " and " + fParamName2  + " are the only parameters" );
}

std::string
ScaleNonLinear::GetName() const{
    return fName;
}

void
ScaleNonLinear::SetName(const std::string& name_){
    fName = name_;
}

