/*************************************************/
/* Wrap an evaluation routine for minuit to call */
/*************************************************/
#ifndef __MINUIT_FCN__
#define __MINUIT_FCN__
#include <TestStatistic.h>
#include <Minuit2/FCNBase.h>
#include <ContainerTools.hpp>
#include <Exceptions.h>

class MinuitFCN : public ROOT::Minuit2::FCNBase
{
public:
   MinuitFCN() : pTestStatistic(NULL), fUp(0.5), fFlipSign(false) {}
   MinuitFCN(TestStatistic *statistic_, const std::set<std::string> &names_) : pTestStatistic(statistic_),
                                                                               fUp(0.5), fFlipSign(false),
                                                                               fParameterNames(names_) {}
   MinuitFCN(const MinuitFCN &other_);
   MinuitFCN &operator=(const MinuitFCN &other_);
   // the second argument above is to absolutely ensure that the parameters are interpreted in the right order

   // these two required by minit
   double Up() const { return fUp; }
   double operator()(const std::vector<double> &paramVals_) const;

   void SetUp(double up_) { fUp = up_; }
   double GetUp() const { return Up(); }

   void SetSignFlip(bool b_) { fFlipSign = b_; }
   bool GetSignFlip() const { return fFlipSign; }

private:
   TestStatistic *pTestStatistic;
   double fUp;
   bool fFlipSign;                       // if true, result is multiplied by -ve 1. changes minimisation to maximisation
   mutable ParameterDict fSetParameters; // edit this map in place to avoid recreating it constantly
   std::set<std::string> fParameterNames;
};
#endif
