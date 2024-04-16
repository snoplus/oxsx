#include <ScaleFunction.h>
#include <sstream>
#include <DoubleParameter.h>
#include <Exceptions.h>
#include <ContainerTools.hpp>

void ScaleFunction::SetScaleFunction(const ScaleFunc& scale_func,
                             const std::vector<std::string>& param_names) {
    /*
    * Setter for the scale function object; also the only way to add parameters
    * into the parameter dict: the FitComponent interface allows one to set 
    * values or rename params, but not add params.
    */
    fScaleFunc = scale_func;
    // Clear any existing param stuff, then add param names with default values
    fParamDict.clear();
    for (auto &&param : param_names) {
        fParamDict[param] = 0;
    }
}

void 
ScaleFunction::Construct(){

    if(fTransObs.GetNObservables() != 1)
        throw RepresentationError("ScaleFunction systematic must have a 1D representation!");
    if ( !fScaleFunc || !fAxes.GetNBins()) {
        throw LogicError("ScaleFunction::Construct(): Tried to construct response matrix without an axis collection!");
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
        const double scaledLow   = fScaleFunc( fParamDict, scaleAxis.GetBinLowEdge(bin) );
        const double scaledHigh  = fScaleFunc( fParamDict, scaleAxis.GetBinHighEdge(bin) );
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


////////////////////////////////////////////////////////////////////////
// Make this object fittable, so the scale func params are adjustable //
////////////////////////////////////////////////////////////////////////

void
ScaleFunction::SetParameter(const std::string& name_, double value){
    if(fParamDict.count(name_) == 0)
        throw ParameterError("ScaleFunction: can't set " + name_);
    fParamDict[name_] = value;
}

double
ScaleFunction::GetParameter(const std::string& name_) const{
    if(fParamDict.count(name_) == 0)
        throw ParameterError("ScaleFunction: parameter " + name_ + "does not exist");
    return fParamDict.at(name_);
}

void
ScaleFunction::SetParameters(const ParameterDict& pd_){
    for (auto &&pair : pd_) {
        try{ fParamDict[pair.first] = pair.second; }
        catch(const std::out_of_range& e_){
            throw ParameterError("ScaleFunction: cannot set parameter " + pair.first);
        }
    }   
}

std::set<std::string>
ScaleFunction::GetParameterNames() const {
    std::set<std::string> set;
    for (auto &&pair : fParamDict) {
        set.insert(pair.first);
    }
    return set;
}

void 
ScaleFunction::RenameParameter(const std::string& old_, const std::string& new_){
    if(fParamDict.count(old_) == 0)
        throw ParameterError("ScaleFunction: can't rename " + old_);
    fParamDict[new_] = fParamDict.at(old_);
    fParamDict.erase(old_);
}

// PRIVATE METHODS

void
ScaleFunction::CacheBinMapping() {
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

double ScaleFunction::GetBinContribution(const BinAxis& scaleAxis, double scaledLow,
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
