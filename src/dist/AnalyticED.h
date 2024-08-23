#ifndef __OXSX_ANALYTIC_ED__
#define __OXSX_ANALYTIC_ED__
#include <EventDistribution.h>
#include <FitComponent.h>
#include <ObsSet.h>
#include <string>

class PDF;
class AnalyticED : public EventDistribution, public FitComponent
{
public:
   AnalyticED(const std::string &, PDF *); // make a copy
   AnalyticED(const AnalyticED &other_);   // deep copy
   ~AnalyticED();                          // frees fFunction
   EventDistribution *Clone() const;

   double Probability(const Event &) const;
   double Probability(const std::vector<double> &) const;

   double Integral() const;
   void Normalise();

   void SetObservables(const ObsSet &);
   ObsSet GetObservables() const;

   unsigned GetNDims() const;

   std::string GetName() const;
   void SetName(const std::string &);

   // FitComponent interface : pass on calls to fFunction, change names
   void SetParameter(const std::string &name_, double value);
   double GetParameter(const std::string &name_) const;

   void SetParameters(const ParameterDict &);
   ParameterDict GetParameters() const;
   size_t GetParameterCount() const;

   std::set<std::string> GetParameterNames() const;
   void RenameParameter(const std::string &old_, const std::string &new_);

private:
   ObsSet fObservables;
   PDF *fFunction;
   double fNorm;
   std::string fName;
};
#endif
