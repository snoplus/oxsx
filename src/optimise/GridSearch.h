#ifndef __GRID_SEARCH__
#define __GRID_SEARCH__
#include <Optimiser.h>
#include <vector>
#include <FitResult.h>
#include <ParameterDict.h>

class TestStatistic;
class GridSearch : public Optimiser
{
public:
   GridSearch() : fMaximising(false) {}
   virtual const FitResult &Optimise(TestStatistic *);

   void SetMinima(const ParameterDict &);
   void SetMaxima(const ParameterDict &);
   void SetStepSizes(const ParameterDict &);

   ParameterDict GetMinima() const;
   ParameterDict GetMaxima() const;
   ParameterDict GetStepSizes() const;

   FitResult GetFitResult() const;
   void SetMaximising(bool);
   bool GetMaximising() const;

private:
   ParameterDict fParamVals;
   FitResult fFitResult;

   bool fMaximising;

   ParameterDict fMinima;
   ParameterDict fMaxima;
   ParameterDict fStepSizes;
   bool Increment(ParameterDict::iterator it_, ParameterDict::iterator end_);
   double fMinVal;
};
#endif
