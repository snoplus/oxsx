#include <SquareRootScale.h>

Function *SquareRootScale::Clone() const
{
  return static_cast<Function *>(new SquareRootScale(*this));
}

double SquareRootScale::operator()(const std::vector<double> &vals_) const
{
  return fGradient * sqrt(abs(vals_.at(0)));
}

// FitComponent class interface
void SquareRootScale::SetParameter(const std::string &name_, double value_)
{
  if (name_ != fParamName)
    throw ParameterError("SquareRootScale: can't set " + name_ + ", " + fParamName + " is the only parameter");
  fGradient = value_;
}

double SquareRootScale::GetParameter(const std::string &name_) const
{
  if (name_ != fParamName)
    throw ParameterError("SquareRootScale: can't get " + name_ + ", " + fParamName + " is the only parameter");
  return fGradient;
}

void SquareRootScale::SetParameters(const ParameterDict &paraDict_)
{
  try
  {
    fGradient = paraDict_.at(fParamName);
  }
  catch (const std::out_of_range &e_)
  {
    throw ParameterError("Set dictionary is missing " + fParamName + ". I did contain: \n" + ContainerTools::ToString(ContainerTools::GetKeys(paraDict_)));
  }
}

ParameterDict SquareRootScale::GetParameters() const
{
  ParameterDict d;
  d[fParamName] = fGradient;
  return d;
}

std::set<std::string> SquareRootScale::GetParameterNames() const
{
  std::set<std::string> names_ = {fParamName};
  return names_;
}

void SquareRootScale::RenameParameter(const std::string &old_, const std::string &new_)
{
  if (old_ != fParamName)
    throw ParameterError("SquareRootScale: can't rename " + old_ + ", " + fParamName + " is the only parameter");
  fParamName = new_;
}
