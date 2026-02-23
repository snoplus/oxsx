#include <Renorm.h>
#include <Exceptions.h>

void Renorm::Construct()
{
    /*
     * Construct the response matrix for this renormalisation systematic.
     * Because impact of this systematic is to scale each bin independently,
     * the response matrix is just diagonal.
     */
    if (!fAxes.GetNBins())
    {
        throw LogicError("Renorm::Construct(): Tried to construct response matrix with an axis collection with wrong bins!");
    }

    fResponse.SetToIdentity();
    fResponse.Scale(fScaleFactor);
}

void Renorm::SetScaleFactor(double scaleFactor_)
{
    fScaleFactor = scaleFactor_;
}

double
Renorm::GetScaleFactor() const
{
    return fScaleFactor;
}

// FitComponent interface: handling systematic's parameters
// You can change the values and names of these parameters

void Renorm::SetParameter(const std::string &name_, double value)
{
    if (fParamDict.count(name_) == 0)
        throw ParameterError("Renorm: can't set " + name_);
    fParamDict[name_] = value;
}

double Renorm::GetParameter(const std::string &name_) const
{
    if (fParamDict.count(name_) == 0)
        throw ParameterError("Renorm: parameter " + name_ + "does not exist");
    return fParamDict.at(name_);
}

void Renorm::SetParameters(const ParameterDict &pd_)
{
    for (auto &&pair : pd_)
    {
        try
        {
            fParamDict[pair.first] = pair.second;
        }
        catch (const std::out_of_range &e_)
        {
            throw ParameterError("Renorm: cannot set parameter " + pair.first);
        }
    }
}

std::set<std::string> Renorm::GetParameterNames() const
{
    std::set<std::string> set;
    for (auto &&pair : fParamDict)
    {
        set.insert(pair.first);
    }
    return set;
}

void Renorm::RenameParameter(const std::string &old_, const std::string &new_)
{
    if (fParamDict.count(old_) == 0)
        throw ParameterError("Renorm: can't rename " + old_);
    fParamDict[new_] = fParamDict.at(old_);
    fParamDict.erase(old_);
}
