#ifndef __OXSX_ED_MANAGER__
#define __OXSX_ED_MANAGER__
#include <FitComponent.h>
#include <ParameterManager.h>
#include <vector>

class EventDistribution;
class Event;

class EDManager : public FitComponent
{
public:
   EDManager() : fNDists(0), fNDims(0), fName("norms") {}
   ~EDManager(); // free pdf copies

   void AddDist(EventDistribution *); // take a copy
   void AddDists(const std::vector<EventDistribution *> &pdfs_);

   const std::vector<double> &GetNormalisations() const;
   void SetNormalisations(const std::vector<double> &);

   double Probability(const Event &) const;

   size_t GetNDists() const;
   size_t GetNDims() const;

   // Make a fittable component - i.e. rescale pdfs inside to fit
   void SetParameter(const std::string &name_, double value);
   double GetParameter(const std::string &name_) const;

   void SetParameters(const ParameterDict &);
   ParameterDict GetParameters() const;
   size_t GetParameterCount() const;

   std::set<std::string> GetParameterNames() const;
   void RenameParameter(const std::string &old_, const std::string &new_);

   std::string GetName() const;
   void SetName(const std::string &);

private:
   ParameterManager fParameterManager;
   size_t fNDists;
   size_t fNDims;
   std::vector<EventDistribution *> fDists;
   std::vector<double> fNormalisations;

   std::string fName; // component name
   void RegisterParameters();
};
#endif
