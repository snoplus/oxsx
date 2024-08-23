#include <Rand.h>
#include <stdlib.h>
#include <TRandom3.h>
#include <iostream>
#include <Exceptions.h>
#include <Function.h>
TRandom3
    Rand::fRandomGen = TRandom3();

double
Rand::Uniform(double max_)
{
  return max_ * fRandomGen.Rndm();
}

double
Rand::UniformRange(double min_, double max_)
{
  return max_ - (max_ - min_) * fRandomGen.Rndm();
}

int Rand::Shoot(int max_)
{
  return int((max_)*Uniform());
}

double
Rand::Gaus(double mean_, double sigma_)
{
  return fRandomGen.Gaus(mean_, sigma_);
}

void Rand::SetSeed(unsigned seed_)
{
  fRandomGen.SetSeed(seed_);
}

unsigned
Rand::GetSeed()
{
  return fRandomGen.GetSeed();
}

int Rand::Poisson(double rate_)
{
  return fRandomGen.Poisson(rate_);
}
double
Rand::VonNeumannSample(Function *f_, double xMin_, double xMax_, double yMax_)
{
  if (xMin_ >= xMax_)
    throw ValueError(Formatter() << "Rand::Sample called with xMin_ = "
                                 << xMin_ << " xMax_ = " << xMax_);

  if (yMax_ <= 0)
    throw ValueError(Formatter() << "Rand::Sample called with yMax_ = "
                                 << yMax_ << ". This is the maximum of the PDF (>0)");
  if (f_->GetNDims() != 1)
    throw ValueError("Von Neumann sampling only implemented for 1D histograms");

  // Von-neuman sampling
  double draw = 0;
  double fOfX = 1;
  double val = 0;

  while (val > fOfX)
  {
    draw = Uniform() * (xMax_ - xMin_) + xMin_;
    fOfX = f_->operator()(std::vector<double>(1, draw));
    val = Uniform() * yMax_;
  }
  return draw;
}
