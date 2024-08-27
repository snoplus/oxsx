#ifndef __OXSX_PARAMETER_MANAGER__
#define __OXSX_PARAMETER_MANAGER__
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <ParameterDict.h>

class FitParameter;
class ParameterManager
{
public:
   ParameterManager() {}
   ~ParameterManager(); // frees fParamPtrs

   // don't copy the other objects parameter pointers, it has its own
   ParameterManager(const ParameterManager &) {}
   ParameterManager operator=(const ParameterManager &) { return ParameterManager(); }

   // Add parameters
   void Add(FitParameter *, const std::string &); // copies
   void AddDouble(double &, const std::string &);
   template <typename Container>
   void AddContainer(Container &,
                     const std::string &nameBase_);
   template <typename Container>
   void AddContainer(Container &,
                     const std::vector<std::string> &names_);
   void Clear();

   // Used to implement FitComponent interface in many places
   virtual void SetParameter(const std::string &name_, double value);
   virtual double GetParameter(const std::string &name_) const;

   virtual void SetParameters(const ParameterDict &);
   ParameterDict GetParameters() const;
   virtual size_t GetParameterCount() const;

   virtual std::set<std::string> GetParameterNames() const;
   virtual void RenameParameter(const std::string &old_, const std::string &new_);

private:
   std::map<std::string, FitParameter *> fParamPtrs;
};

// implements AddContainerOfParameters(-)
#include <ParameterManager.hpp>
#endif
