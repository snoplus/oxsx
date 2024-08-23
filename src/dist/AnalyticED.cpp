#include <AnalyticED.h>
#include <Exceptions.h>
#include <EventDistribution.h>
#include <PDF.h>

AnalyticED::AnalyticED(const std::string &name_, PDF *f_)
{
    fName = name_;
    fFunction = dynamic_cast<PDF *>(f_->Clone());
    fNorm = 1;
}

AnalyticED::~AnalyticED()
{
    delete fFunction;
}

AnalyticED::AnalyticED(const AnalyticED &other_)
{
    fName = other_.fName;
    fNorm = other_.fNorm;
    fObservables = other_.fObservables;
    fFunction = dynamic_cast<PDF *>(other_.fFunction->Clone());
}

EventDistribution *
AnalyticED::Clone() const
{
    return static_cast<EventDistribution *>(new AnalyticED(*this));
}

double
AnalyticED::Probability(const std::vector<double> &vals_) const
{
    try
    {
        return fFunction->operator()(vals_) / fNorm;
    }
    catch (const DimensionError &e_)
    {
        throw DimensionError(std::string("AnalyticED internal function ::") + e_.what());
    }
}

double
AnalyticED::Probability(const Event &event_) const
{
    try
    {
        return Probability(event_.ToObsSet(fObservables));
    }

    catch (const RepresentationError &e_)
    {
        throw RepresentationError("AnalyticED::Probability() failed with  " + std::string(e_.what()) +
                                  " is the rep set correctly?");
    }
}

double
AnalyticED::Integral() const
{
    return fNorm;
}

void AnalyticED::Normalise()
{
    fNorm = 1;
}

void AnalyticED::SetObservables(const ObsSet &rep_)
{
    fObservables = rep_;
}

ObsSet
AnalyticED::GetObservables() const
{
    return fObservables;
}

unsigned
AnalyticED::GetNDims() const
{
    return fFunction->GetNDims();
}

std::string
AnalyticED::GetName() const
{
    return fName;
}

void AnalyticED::SetName(const std::string &name_)
{
    fName = name_;
}

// Fitting this dist to data means adjusting the underlying function

void AnalyticED::RenameParameter(const std::string &old_, const std::string &new_)
{
    fFunction->RenameParameter(old_, new_);
}

void AnalyticED::SetParameter(const std::string &name_, double value_)
{
    fFunction->SetParameter(name_, value_);
}

double
AnalyticED::GetParameter(const std::string &name_) const
{
    return fFunction->GetParameter(name_);
}

void AnalyticED::SetParameters(const ParameterDict &ps_)
{
    try
    {
        fFunction->SetParameters(ps_);
    }
    catch (const ParameterError &e_)
    {
        throw ParameterError("AnalyticED internal function: " + std::string(e_.what()));
    }
}

ParameterDict
AnalyticED::GetParameters() const
{
    return fFunction->GetParameters();
}

size_t
AnalyticED::GetParameterCount() const
{
    return fFunction->GetParameterCount();
}

std::set<std::string>
AnalyticED::GetParameterNames() const
{
    return fFunction->GetParameterNames();
}
