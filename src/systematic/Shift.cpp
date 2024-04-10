#include <Shift.h>
#include <sstream>
#include <DoubleParameter.h>
#include <Exceptions.h>
#include <ContainerTools.hpp>

void 
Shift::Construct(){
    
    if(fTransObs.GetNObservables() != 1)
        throw RepresentationError("Shift systematic must have a 1D representation!");
    if (!fAxes.GetNBins()) {
        throw LogicError("Shift::Construct(): Tried to construct response matrix without an axis collection!");
    }

    // If haven't already, generate mapping from full pdf bin IDs
    //--> transforming subspace bin IDs 
    if (!fCachedBinMapping) { CacheBinMapping(); }

    // Loop over bins of the shifting axis, and get the user-defined 
    // shift values at each bin centre.
    // the axis to shift
    const std::string&  shiftAxisName = fTransObs.GetNames().at(0);
    const BinAxis& shiftAxis          = fAxes.GetAxis(fDistObs.GetIndex(shiftAxisName));
    // Each (pre-shifted) bin gets a mapping from (post-shifted) bins to non-zero values
    std::vector<std::map<size_t,double>> shift_vals(shiftAxis.GetNBins(), std::map<size_t,double>{});
    for (size_t bin = 0; bin < shiftAxis.GetNBins(); bin++) {
        // First - work out edges of shifted bin interval along axis
        const double shiftLow   = shiftAxis.GetBinLowEdge(bin)  + fShift;
        const double shiftHigh  = shiftAxis.GetBinHighEdge(bin) + fShift;
        // Then - which unshifted bins do these edges lie in?
        const size_t shiftLowIndex = shiftAxis.FindBin(shiftLow);
        const size_t shiftHighIndex = shiftAxis.FindBin(shiftHigh);
        // Go through bins in this range, to find contribution to each
        for (size_t bin_test = shiftLowIndex; bin_test <= shiftHighIndex; bin_test++) {
            const double contribution = GetBinContribution(shiftAxis, shiftLow, shiftHigh, bin_test);
            if (contribution > 0) { shift_vals[bin][bin_test] = contribution; }
        }
    }

    // Finally, construct the full response matrix by setting the diagonal
    // elements to their appropriate value, using the bin ID mapping to help.
    fResponse.SetZeros();
    for (size_t obs_bin_id = 0; obs_bin_id < fAxes.GetNBins(); obs_bin_id++) {
        const size_t shiftAxisBin = fDistTransBinMapping.at(obs_bin_id);
        for (const auto& postShiftBinPair : shift_vals.at(shiftAxisBin)) {
            const size_t shifted_bin_id = fMappingDistAndTrans.GetComponent(obs_bin_id, postShiftBinPair.first);
            fResponse.SetComponent(shifted_bin_id, obs_bin_id, postShiftBinPair.second);
        }
    }
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


// PRIVATE METHODS

void Shift::CacheBinMapping() {
    /*
    * Because this shift systematic might only require a subset of the
    * observables to be defined, but still need to act on the whole event
    * distribution, we need a way of mapping from the full event distribution
    * to the subspace we actually want to calculate shift values over.
    */
    const size_t relativeIndex = fTransObs.GetRelativeIndices(fDistObs).at(0);

    // cache the equivalent index in the binning system of the systematic
    fDistTransBinMapping.resize(fAxes.GetNBins());
    for(size_t i = 0; i < fAxes.GetNBins(); i++) {
        fDistTransBinMapping[i] = fAxes.UnflattenIndex(i, relativeIndex);
    }
    /*
     * Given a bin idx in the full axis collection (1), and another on the
     * shifting axis (2), this matrix stores the bin idx of the bin with the
     * same position as (1), excepting for the shifting axis, in which it takes
     * (2)'s position.
    */
    const std::string&  shiftAxisName = fTransObs.GetNames().at(0);
    const BinAxis& shiftAxis          = fAxes.GetAxis(fDistObs.GetIndex(shiftAxisName));

    fMappingDistAndTrans = DenseMatrix(fAxes.GetNBins(),shiftAxis.GetNBins());

    for (size_t obs_bin_id = 0; obs_bin_id < fAxes.GetNBins(); obs_bin_id++) {
        const auto obs_bins_unflat = fAxes.UnpackIndices(obs_bin_id);

        for (size_t shift_bin_id = 0; shift_bin_id < shiftAxis.GetNBins(); shift_bin_id++) {
            auto shifted_bins_unflat = obs_bins_unflat;
            shifted_bins_unflat[fDistObs.GetIndex(shiftAxisName)] = shift_bin_id;
            const size_t idx = fAxes.FlattenIndices(shifted_bins_unflat);
            fMappingDistAndTrans.SetComponent(obs_bin_id, shift_bin_id, idx);
        }
    }

    fCachedBinMapping = true;
}

double Shift::GetBinContribution(const BinAxis& shiftAxis, double shiftedLow,
                                 double shiftedHigh, size_t bin_test) {
    /*
     * Some logic to work out the fraction of the shifted bin interval
     * within a given test bin interval
     */
    const double shiftedWidth = shiftedHigh - shiftedLow;
    const double testLow  = shiftAxis.GetBinLowEdge(bin_test);
    const double testHigh = shiftAxis.GetBinHighEdge(bin_test);
    // Guide:  | | = shifted bin interval, /  / = test bin interval
    // Is it in the shift region at all?
    // |   |   /   /        OR    /  /   |  |
    if (testLow > shiftedHigh || testHigh < shiftedLow) {
        return 0.;
    } else{
        // Is it fully in the region?
        bool includedFromBelow = testLow > shiftedLow;
        bool includedFromAbove = testHigh < shiftedHigh;

        if (includedFromBelow) {
            if (includedFromAbove) {
                // |  /   /  | : test fully inside shifted bin
                return (testHigh - testLow)/shiftedWidth;
            } else {
                // |  /   |  / : shifted hangs off low end of test only
                return (shiftedHigh - testLow)/shiftedWidth;
            }
        } else {
            if (includedFromAbove) {
                // /  |   /  | : shifted hangs off of high end only
                return (testHigh - shiftedLow)/shiftedWidth;
            } else {
                // /  |   |  / : shifted fully within test bin
                return 1;
            }
        }
    }
}
