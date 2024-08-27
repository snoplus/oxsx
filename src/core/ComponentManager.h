#ifndef __OXSX_COMPONENT_MANAGER__
#define __OXSX_COMPONENT_MANAGER__
#include <FitComponent.h>
#include <vector>
#include <set>
#include <string>

class ComponentManager
{
public:
   ComponentManager() : fComponentCount(0) {}
   void AddComponent(FitComponent *);

   void SetParameters(const ParameterDict &);
   std::set<std::string> GetParameterNames() const;
   ParameterDict GetParameters() const;
   int GetTotalParameterCount() const;
   size_t GetComponentCount() const;
   void Clear();

   double GetParameter(const std::string &) const;

private:
   std::vector<FitComponent *> fComponents;
   int fComponentCount;
};
#endif
