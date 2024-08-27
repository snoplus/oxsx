#include <Shape.h>
#include <Exceptions.h>

void Shape::SetShapeFunction(const ShapeFunction &shape_func,
                             const std::vector<std::string> &param_names)
{
    /*
     * Setter for the shape function object; also the only way to add parameters
     * into the parameter dict: the FitComponent interface allows one to set
     * values or rename params, but not add params.
     */
    fShapeFunc = shape_func;
    // Clear any existing param stuff, then add param names with default values
    fParamDict.clear();
    for (auto &&param : param_names)
    {
        fParamDict[param] = 0;
    }
}

void Shape::Construct()
{
    /*
     * Construct the response matrix for this shape systematic, given the
     * shape function given by the user.
     * Because impact of this systematic is to scale each bin independently,
     * the response matrix is just diagonal.
     */
    if (!fShapeFunc || !fAxes.GetNBins())
    {
        throw LogicError("Shape::Construct(): Tried to construct response matrix without either a shape function defined, or an axis collection!");
    }
    // If haven't already, generate mapping from full pdf bin IDs
    //--> transforming subspace bin IDs
    if (!fCachedBinMapping)
    {
        CacheBinMapping();
    }
    // Loop over bins in the transformation subspace, and get the user-defined
    // shape scale factors at each bin centre.
    std::vector<double> bin_centres(fTransAxes.GetNDimensions());
    std::vector<double> shape_vals;
    shape_vals.reserve(fTransAxes.GetNBins());

    for (size_t bin = 0; bin < fTransAxes.GetNBins(); bin++)
    {
        fTransAxes.GetBinCentres(bin, bin_centres);
        const double scale = fShapeFunc(fParamDict, bin_centres);
        if (scale < 0)
        {
            throw ValueError("User-defined shape function returned negative value: " + std::to_string(scale));
        }
        shape_vals.push_back(scale);
    }

    // Finally, construct the full response matrix by setting the diagonal
    // elements to their appropriate value, using the bin ID mapping to help.
    for (size_t obs_bin_id = 0; obs_bin_id < fAxes.GetNBins(); obs_bin_id++)
    {
        fResponse.SetComponent(obs_bin_id, obs_bin_id,
                               shape_vals.at(fDistTransBinMapping.at(obs_bin_id)));
    }
}

// FitComponent interface: handling systematic's parameters
// Note that once the parameter names have been set in SetShapeFunction(),
// The number of parameters cannot be changed with this interface.
// You can change the values and names of these parameters

void Shape::SetParameter(const std::string &name_, double value)
{
    if (fParamDict.count(name_) == 0)
        throw ParameterError("Shape: can't set " + name_);
    fParamDict[name_] = value;
}

double Shape::GetParameter(const std::string &name_) const
{
    if (fParamDict.count(name_) == 0)
        throw ParameterError("Shape: parameter " + name_ + "does not exist");
    return fParamDict.at(name_);
}

void Shape::SetParameters(const ParameterDict &pd_)
{
    for (auto &&pair : pd_)
    {
        try
        {
            fParamDict[pair.first] = pair.second;
        }
        catch (const std::out_of_range &e_)
        {
            throw ParameterError("Shift: cannot set parameter " + pair.first);
        }
    }
}

std::set<std::string> Shape::GetParameterNames() const
{
    std::set<std::string> set;
    for (auto &&pair : fParamDict)
    {
        set.insert(pair.first);
    }
    return set;
}

void Shape::RenameParameter(const std::string &old_, const std::string &new_)
{
    if (fParamDict.count(old_) == 0)
        throw ParameterError("Shift: can't rename " + old_);
    fParamDict[new_] = fParamDict.at(old_);
    fParamDict.erase(old_);
}

// PRIVATE METHODS

void Shape::CacheBinMapping()
{
    /*
     * Because this shape systematic might only require a subset of the
     * observables to be defined, but still need to act on the whole event
     * distribution, we need a way of mapping from the full event distribution
     * to the subspace we actually want to calculate shape values over.
     */
    std::vector<size_t> relativeIndices = fTransObs.GetRelativeIndices(fDistObs);
    const AxisCollection &axes = fAxes;

    //  get the axes that this systematic will act on
    fTransAxes = AxisCollection();
    for (size_t i = 0; i < relativeIndices.size(); i++)
    {
        fTransAxes.AddAxis(axes.GetAxis(relativeIndices.at(i)));
    }
    // cache the equivilent index in the binning system of the systematic
    fDistTransBinMapping.resize(fAxes.GetNBins());
    std::vector<size_t> sysIndices(relativeIndices.size(), 0);
    for (size_t i = 0; i < axes.GetNBins(); i++)
    {
        for (size_t dim = 0; dim < relativeIndices.size(); dim++)
        {
            sysIndices[dim] = axes.UnflattenIndex(i, relativeIndices.at(dim));
        }
        fDistTransBinMapping[i] = fTransAxes.FlattenIndices(sysIndices);
    }
    fCachedBinMapping = true;
}