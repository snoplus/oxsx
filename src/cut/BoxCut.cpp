#include <BoxCut.h>
#include <Exceptions.h>

bool BoxCut::PassesCut(const Event &ev_) const
{
  double val = 0;
  try
  {
    val = ev_.GetDatum(fObs);
  }
  catch (const NotFoundError &)
  {
    throw NotFoundError("Cut::Cut to non-existent data observable requested!");
  }

  return (val < fUpperLim && val > fLowerLim);
}

Cut *BoxCut::Clone() const
{
  return new BoxCut(*this);
}

std::string
BoxCut::GetName() const
{
  return fName;
}
