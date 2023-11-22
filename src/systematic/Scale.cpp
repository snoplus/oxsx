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
    if (!fAxes.GetNBins()) {
        throw LogicError("Scale::Construct(): Tried to construct response matrix without an axis collection!");
    }
    // If haven't already, generate mapping from full pdf bin IDs
    //--> transforming subspace bin IDs 
    if (!fCachedBinMapping) { CacheBinMapping(); }

    // Loop over bins of the scaling axis, and get the user-defined 
    // shape scale factors at each bin centre.
    // the axis to scale
    const std::string&  scaleAxisName = fTransObs.GetNames().at(0);
    const BinAxis& scaleAxis          = fAxes.GetAxis(fDistObs.GetIndex(scaleAxisName));
    // Each (pre-scaled) bin gets a mapping from (post-scaled) bins to non-zero values
    std::vector<std::map<size_t,double>> scale_vals(scaleAxis.GetNBins(), std::map<size_t,double>{});
    for (size_t bin = 0; bin < scaleAxis.GetNBins(); bin++) {
        // First - work out edges of scaled bin interval along axis
        const double scaledLow   = scaleAxis.GetBinLowEdge(bin)  * fScaleFactor;
        const double scaledHigh  = scaleAxis.GetBinHighEdge(bin) * fScaleFactor;
        // Then - which unscaled bins do these edges lie in?
        const size_t scaledLowIndex = scaleAxis.FindBin(scaledLow);
        const size_t scaledHighIndex = scaleAxis.FindBin(scaledHigh);
        // Go through bins in this range, to find contribution to each
        for (size_t bin_test = scaledLowIndex; bin_test <= scaledHighIndex; bin_test++) {
            const double contribution = GetBinContribution(scaleAxis, scaledLow, scaledHigh, bin_test);
            if (contribution > 0) { scale_vals[bin][bin_test] = contribution; }
        }
    }

    // Finally, construct the full response matrix by setting the diagonal
    // elements to their appropriate value, using the bin ID mapping to help.
    fResponse.SetZeros();
    for (size_t obs_bin_id = 0; obs_bin_id < fAxes.GetNBins(); obs_bin_id++) {
        const size_t scaleAxisBin = fDistTransBinMapping.at(obs_bin_id);
        for (const auto& postScaleBinPair : scale_vals.at(scaleAxisBin)) {
            const size_t scaled_bin_id = fMappingDistAndTrans.GetComponent(obs_bin_id, postScaleBinPair.first);
            fResponse.SetComponent(scaled_bin_id, obs_bin_id, postScaleBinPair.second);
        }
    }
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

// PRIVATE METHODS

void Scale::CacheBinMapping() {
    /*
    * Because this shape systematic might only require a subset of the
    * observables to be defined, but still need to act on the whole event
    * distribution, we need a way of mapping from the full event distribution
    * to the subspace we actually want to calculate shape values over.
    */
    const size_t relativeIndex = fTransObs.GetRelativeIndices(fDistObs).at(0);

    // cache the equivalent index in the binning system of the systematic
    fDistTransBinMapping.resize(fAxes.GetNBins());
    for(size_t i = 0; i < fAxes.GetNBins(); i++) {
        fDistTransBinMapping[i] = fAxes.UnflattenIndex(i, relativeIndex);
    }
    /*
     * Given a bin idx in the full axis collection (1), and another on the
     * scaling axis (2), this matrix stores the bin idx of the bin with the
     * same position as (1), excepting for the scaling axis, in which it takes
     * (2)'s position.
    */
    const std::string&  scaleAxisName = fTransObs.GetNames().at(0);
    const BinAxis& scaleAxis          = fAxes.GetAxis(fDistObs.GetIndex(scaleAxisName));

    fMappingDistAndTrans = DenseMatrix(fAxes.GetNBins(),scaleAxis.GetNBins());

    for (size_t obs_bin_id = 0; obs_bin_id < fAxes.GetNBins(); obs_bin_id++) {
        const auto obs_bins_unflat = fAxes.UnpackIndices(obs_bin_id);

        for (size_t scale_bin_id = 0; scale_bin_id < scaleAxis.GetNBins(); scale_bin_id++) {
            auto scaled_bins_unflat = obs_bins_unflat;
            scaled_bins_unflat[fDistObs.GetIndex(scaleAxisName)] = scale_bin_id;
            const size_t idx = fAxes.FlattenIndices(scaled_bins_unflat);
            fMappingDistAndTrans.SetComponent(obs_bin_id, scale_bin_id, idx);
        }
    }

    fCachedBinMapping = true;
}

double Scale::GetBinContribution(const BinAxis& scaleAxis, double scaledLow,
                                 double scaledHigh, size_t bin_test) {
    /*
     * Some logic to work out the fraction of the scaled bin interval
     * within a given test bin interval
     */
    const double scaledWidth = scaledHigh - scaledLow;
    const double testLow  = scaleAxis.GetBinLowEdge(bin_test);
    const double testHigh = scaleAxis.GetBinHighEdge(bin_test);
    // Guide:  | | = scaled bin interval, /  / = test bin interval
    // Is it in the scale region at all?
    // |   |   /   /        OR    /  /   |  |
    if (testLow > scaledHigh || testHigh < scaledLow) {
        return 0.;
    } else{
        // Is it fully in the region?
        bool includedFromBelow = testLow > scaledLow;
        bool includedFromAbove = testHigh < scaledHigh;

        if (includedFromBelow) {
            if (includedFromAbove) {
                // |  /   /  | : test fully inside scaled bin
                return (testHigh - testLow)/scaledWidth;
            } else {
                // |  /   |  / : scaled hangs off low end of test only
                return (scaledHigh - testLow)/scaledWidth;
            }
        } else {
            if (includedFromAbove) {
                // /  |   /  | : scaled hangs off of high end only
                return (testHigh - scaledLow)/scaledWidth;
            } else {
                // /  |   |  / : scaled fully within test bin
                return 1;
            }
        }
    }
}
