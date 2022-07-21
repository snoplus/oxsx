#include <Shift.h>
#include <sstream>
#include <DoubleParameter.h>
#include <Exceptions.h>
#include <ContainerTools.hpp>

void 
Shift::Construct(){
    
    if(fTransObs.GetNObservables() != 1)
        throw RepresentationError("Shift systematic must have a 1D representation!");

    const AxisCollection& axes       = fAxes;
    // the axis to shift
    const std::string&  shiftAxisName = fTransObs.GetNames().at(0);
    const BinAxis& shiftAxis          = axes.GetAxis(fDistObs.GetIndex(shiftAxisName));

    const size_t nBins               = axes.GetNBins(); 
    const size_t shiftAxisNBins      = shiftAxis.GetNBins(); 
 
    for(size_t i = 0; i < nBins; i++){
        // For each old bin, work out the contributions into all of the new bins
        // indices in other components should be unaffected
        std::vector<size_t> oldIndices = axes.UnpackIndices(i);
        size_t shiftBin                = oldIndices.at(fDistObs.GetIndex(shiftAxisName));
        
        double shiftedLow   = shiftAxis.GetBinLowEdge(shiftBin)  + fShift;
        double shiftedHigh  = shiftAxis.GetBinHighEdge(shiftBin) + fShift;
        double shiftedWidth = shiftedHigh - shiftedLow;

        // new bin to map into, mapping only happens if the indices are the same except the one to 
        // shift so, loop over the bins in the shift axes and leave other indices the same
        // the others are zero from initialisation 
        
        std::vector<size_t> newIndices = oldIndices;
        for(size_t j = 0; j < shiftAxisNBins; j++){
            newIndices[fDistObs.GetIndex(shiftAxisName)] = j;
            size_t newShiftBin = j;
                        
            double newLow  = shiftAxis.GetBinLowEdge(newShiftBin);
            double newHigh = shiftAxis.GetBinHighEdge(newShiftBin);

            double contribution;
            // Is it in the shift region at all?
            if (newLow > shiftedHigh || newHigh < shiftedLow) 
                contribution = 0;

            else{
                // Is it fully in the region?
                bool includedFromBelow = newLow > shiftedLow;
                bool includedFromAbove = newHigh < shiftedHigh;

                if (includedFromBelow && includedFromAbove) 
		    // fully inside
		    contribution = (newHigh - newLow)/shiftedWidth;

                else if (includedFromBelow) 
                    // spills partly over the top
                    contribution = (shiftedHigh - newLow)/shiftedWidth;

                else 
                    // spills partly over the bottom
                    contribution = (newHigh - shiftedLow)/shiftedWidth;
            }

            fResponse.SetComponent(axes.FlattenIndices(newIndices), i, contribution);
        }
               
    }
    return;
}

void
Shift::SetShift(double shift_){
    fShift = shift_;
}

double
Shift::GetShift() const{
    return fShift;
}


//////////////////////////////////////////////////////////////////
// Make this object fittable, so the shift is adjustable //
//////////////////////////////////////////////////////////////////

void
Shift::SetParameter(const std::string& name_, double value){
    if(name_ != fParamName)
        throw ParameterError("Shift: can't set " + name_ + ", " + fParamName + " is the only parameter" );
    fShift = value;
}

double
Shift::GetParameter(const std::string& name_) const{
   if(name_ != fParamName)
        throw ParameterError("Shift: can't get " + name_ + ", " + fParamName + " is the only parameter" );
   return fShift;
}

void
Shift::SetParameters(const ParameterDict& pd_){
    try{
        fShift = pd_.at(fParamName);
    }
    catch(const std::out_of_range& e_){
        throw ParameterError("Set dictionary is missing " + fParamName + ". I did contain: \n" + ContainerTools::ToString(ContainerTools::GetKeys(pd_)));
    }
}

ParameterDict
Shift::GetParameters() const{
    ParameterDict d;
    d[fParamName] = fShift;
    return d;
}

size_t
Shift::GetParameterCount() const{
    return 1;
}

std::set<std::string>
Shift::GetParameterNames() const{
    std::set<std::string> set;
    set.insert(fParamName);
    return set;
}

void
Shift::RenameParameter(const std::string& old_, const std::string& new_){
    if(old_ != fParamName)
        throw ParameterError("Shift: can't rename " + old_ + ", " + fParamName + " is the only parameter" );
    fParamName = new_;
}

std::string
Shift::GetName() const{
    return fName;
}

void
Shift::SetName(const std::string& name_){
    fName = name_;
}

