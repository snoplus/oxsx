#include <Scale.h>
#include <sstream>
#include <DoubleParameter.h>
#include <Exceptions.h>
#include <ContainerTools.hpp>

void 
Scale::Construct(){
    if (fScaleFactor <= 0)
        throw ValueError("Scale factor must be >0 !");
    
    if(fTransObs.GetNObservables() != 1)
        throw RepresentationError("Scale systematic must have a 1D representation!");

    const AxisCollection& axes       = fAxes;
    // the axis to scale
    const std::string&  scaleAxisName = fTransObs.GetNames().at(0);
    const BinAxis& scaleAxis          = axes.GetAxis(fDistObs.GetIndex(scaleAxisName));


    const size_t nBins               = axes.GetNBins(); 
    const size_t scaleAxisNBins      = scaleAxis.GetNBins(); 
 
    for(size_t i = 0; i < nBins; i++){
        // For each old bin, work out the contributions into all of the new bins
        // indices in other components should be unaffected
        std::vector<size_t> oldIndices = axes.UnpackIndices(i);
        size_t scaleBin                = oldIndices.at(fDistObs.GetIndex(scaleAxisName));
        
        double scaledLow   = scaleAxis.GetBinLowEdge(scaleBin)  * fScaleFactor;
        double scaledHigh  = scaleAxis.GetBinHighEdge(scaleBin) * fScaleFactor;
        double scaledWidth = scaledHigh - scaledLow;

        // new bin to map into, mapping only happens if the indices are the same except the one to 
        // scale so, loop over the bins in the scale axes and leave other indices the same
        // the others are zero from initialisation 
        
        std::vector<size_t> newIndices = oldIndices;
        for(size_t j = 0; j < scaleAxisNBins; j++){
            newIndices[fDistObs.GetIndex(scaleAxisName)] = j;
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
                    contribution = 1/GetScaleFactor();

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
Scale::SetScaleFactor(double scaleFactor_){
    fScaleFactor = scaleFactor_;
}

double
Scale::GetScaleFactor() const{
    return fScaleFactor;
}


//////////////////////////////////////////////////////////////////
// Make this object fittable, so the scale factor is adjustable //
//////////////////////////////////////////////////////////////////

void
Scale::SetParameter(const std::string& name_, double value){
    if(name_ != fParamName)
        throw ParameterError("Scale: can't set " + name_ + ", " + fParamName + " is the only parameter" );
    fScaleFactor = value;
}

double
Scale::GetParameter(const std::string& name_) const{
   if(name_ != fParamName)
        throw ParameterError("Scale: can't get " + name_ + ", " + fParamName + " is the only parameter" );
   return fScaleFactor;
}

void
Scale::SetParameters(const ParameterDict& pd_){
    try{
        fScaleFactor = pd_.at(fParamName);
    }
    catch(const std::out_of_range& e_){
        throw ParameterError("Set dictionary is missing " + fParamName + ". I did contain: \n" + ContainerTools::ToString(ContainerTools::GetKeys(pd_)));
    }
}

ParameterDict
Scale::GetParameters() const{
    ParameterDict d;
    d[fParamName] = fScaleFactor;
    return d;
}

size_t
Scale::GetParameterCount() const{
    return 1;
}

std::set<std::string>
Scale::GetParameterNames() const{
    std::set<std::string> set;
    set.insert(fParamName);
    return set;
}

void
Scale::RenameParameter(const std::string& old_, const std::string& new_){
    if(old_ != fParamName)
        throw ParameterError("Scale: can't rename " + old_ + ", " + fParamName + " is the only parameter" );
    fParamName = new_;
}

std::string
Scale::GetName() const{
    return fName;
}

void
Scale::SetName(const std::string& name_){
    fName = name_;
}

