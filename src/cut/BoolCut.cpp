#include <BoolCut.h>
#include <Exceptions.h>
#include <Event.h>

bool BoolCut::PassesCut(const Event &ev_) const
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

  return val == fVal;
}

Cut *BoolCut::Clone() const
{
  return new BoolCut(*this);
}

std::string
BoolCut::GetName() const
{
  return fName;
}
