/*****************************************************************************************************/
/* Central place for control of random numbers, at the moment dumb c++ rand() but switch it out here */
/* to change everywhere                                                                              */
/*****************************************************************************************************/
#ifndef __OXSX_RAND__
#define __OXSX_RAND__
#include <TRandom3.h>
class Function;
class Rand
{
public:
   static double Uniform(double max_ = 1);
   static double UniformRange(double min_ = 0, double max_ = 1);
   static double Gaus(double mean_ = 0, double sigma_ = 1);
   static void SetSeed(unsigned seed_ = 0);
   static unsigned GetSeed();
   static int Poisson(double rate_);
   static int Shoot(int max_);
   // Fix me : currently just 1D
   static double VonNeumannSample(Function *f_, double xMin_, double xMax_, double yMax_);

private:
   static TRandom3 fRandomGen;
};
#endif
