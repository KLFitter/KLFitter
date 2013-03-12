#include "ResDoubleGaussBase.h"
#include <iostream>

#include <cmath>

// --------------------------------------------------------- 
KLFitter::ResDoubleGaussBase::ResDoubleGaussBase(const char * filename) : KLFitter::ResolutionBase(10)
{
  // read parameters from file 
  ReadParameters(filename, 10); 
}

// ---------------------------------------------------------
KLFitter::ResDoubleGaussBase::ResDoubleGaussBase(std::vector<double> const& parameters) : KLFitter::ResolutionBase(parameters)
{
  // check number of parameters 
  if (parameters.size() != 10)
    {
      std::cout << "KLFitter::ResDoubleGaussBase::ResDoubleGaussBase(). Number of parameters != 10." << std::endl;
      return;
    } 
}

// --------------------------------------------------------- 
KLFitter::ResDoubleGaussBase::~ResDoubleGaussBase()
{
  ;
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussBase::GetSigma(double xmeas){
  /* Calculate mean width of both gaussians; weight the width of the 2nd one with its amplitude */
  double sigma1 = GetSigma1(xmeas);
  double sigma2 = GetSigma2(xmeas); 
  double amplitude2 = GetAmplitude2(xmeas); 
  double sigma = (sigma1 + amplitude2*sigma2) / (1+amplitude2);

  /* sigma estimates the fractional resolution, but we want absolute */
  return sigma*xmeas;
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussBase::p(double x, double xmeas, bool &good)
{
  double m1 = GetMean1(x); 
  double s1 = GetSigma1(x);
  double a2 = GetAmplitude2(x); 
  double m2 = GetMean2(x); 
  double s2 = GetSigma2(x); 

  // sanity checks for p2, p3 and p5
  good = CheckDoubleGaussianSanity(s1, a2, s2);

  double dx = (x - xmeas) / x; 
 
  // calculate double-Gaussian 
  return 1./sqrt(2.*M_PI) / (s1 + a2 * s2) * ( exp( -(dx-m1)*(dx-m1)/(2 * s1*s1) ) + a2 * exp( -(dx-m2)*(dx-m2)/(2 * s2 * s2) ) ); 
}

// --------------------------------------------------------- 

