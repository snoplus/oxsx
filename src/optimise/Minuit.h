/************************************************************************************/
/* Largely inspired by the similar class in RAT, written by P.G.Jones and M.Mottram */
/************************************************************************************/

#ifndef __OXSX_MINUIT__
#define __OXSX_MINUIT__
#include <Optimiser.h>
#include <string>
#include <vector>
#include <MinuitFCN.h>
#include <Minuit2/MnApplication.h>
#include <Minuit2/FunctionMinimum.h>
#include <FitResult.h>
#include <DenseMatrix.h>
#include <set>

class TestStatistic;

class Minuit : public Optimiser
{
public:
   Minuit() : fMaxCalls(0),
              fTolerance(0.1),
              fMethod("Migrad"),
              fMinimiser(NULL),
              fMaximising(false) {}
   ~Minuit();

   virtual const FitResult &Optimise(TestStatistic *);

   void Fix(const std::string &param_);
   void Release(const std::string &param_);

   void SetMethod(const std::string &);
   std::string GetMethod() const;

   void SetInitialValues(const ParameterDict &);
   void SetInitialErrors(const ParameterDict &);

   void SetUpperContourEdge(double);
   double GetUpperContourEdge() const;

   void SetMinima(const ParameterDict &minima_);
   ParameterDict GetMinima() const;

   void SetMaxima(const ParameterDict &maxima_);
   ParameterDict GetMaxima() const;

   void SetMaxCalls(unsigned);
   unsigned GetMaxCalls() const;

   void SetTolerance(double);
   double GetTolerance() const;

   void SetStrategy(int);
   int GetStrategy() const;
  
   void SetMaximising(bool b_) { fMaximising = b_; }
   bool GetMaximising() const { return fMaximising; }

   FitResult GetFitResult() const;

private:
   void Initialise(TestStatistic *);
   MinuitFCN fMinuitFCN; // wrapper on evaluator so migrad can call it
   ParameterDict fInitialValues;
   ParameterDict fInitialErrors;

   ParameterDict fMinima;
   ParameterDict fMaxima;
   std::set<std::string> fFixedParameters;

   unsigned fMaxCalls;
   double fTolerance;
   int fStrategy = 2;

   std::string fMethod;
   ROOT::Minuit2::MnApplication *fMinimiser;

   std::set<std::string> fParameterNames;
   // the order they are held in vectors for ROOT

   FitResult fFitResult;
   bool fMaximising;

   DenseMatrix CalcCovarianceMatrix(ROOT::Minuit2::FunctionMinimum);
};
#endif
