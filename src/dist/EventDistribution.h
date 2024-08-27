#ifndef __OXSX_EVENT_DISTRIBUTION__
#define __OXSX_EVENT_DISTRIBUTION__
#include <ObsSet.h>
#include <Event.h>
#include <vector>
#include <string>

class EventDistribution
{
public:
    virtual ~EventDistribution() {}
    virtual EventDistribution *Clone() const = 0;

    virtual double Probability(const Event &) const = 0;

    virtual double Integral() const = 0;
    virtual void Normalise() = 0;

    virtual unsigned GetNDims() const = 0;

    virtual std::string GetName() const = 0;
    virtual void SetName(const std::string &) = 0;
};
#endif
