#include "ResDoubleGaussE_1.h" 
#include <iostream>

#include <cmath>

// --------------------------------------------------------- 
KLFitter::ResDoubleGaussE_1::ResDoubleGaussE_1(const char * filename) : KLFitter::ResolutionBase(10)
{
  // read parameters from file 
  ReadParameters(filename, 10); 
}

// ---------------------------------------------------------
KLFitter::ResDoubleGaussE_1::ResDoubleGaussE_1(std::vector<double> const& parameters) : KLFitter::ResolutionBase(parameters)
{
  // check number of parameters 
  if (parameters.size() != 10)
    {
      std::cout << "KLFitter::ResDoubleGaussE_1::ResDoubleGaussE_1(). Number of parameters != 10." << std::endl;
      return;
    } 
}

// --------------------------------------------------------- 
KLFitter::ResDoubleGaussE_1::~ResDoubleGaussE_1()
{
  ;
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_1::p(double x, double xmeas)
{
  double p1 = fParameters[0] + fParameters[1] * x; 
  double p2 = fParameters[2] / sqrt(x) + fParameters[3];
  double p3 = fParameters[4] + fParameters[5] * x; 
  double p4 = fParameters[6] + fParameters[7] * x; 
  double p5 = fParameters[8] + fParameters[9] * x; 

  // sanity checks for p3 and p5
  ResDoubleGaussE_1::CheckDoubleGaussianSanity(p2, p3, p5);

  double dx = (x - xmeas) / x; 
 
  // calculate double-Gaussian 
  double prob= 0.398942280401 / (p2 + p3 * p5) * ( exp( -(dx-p1)*(dx-p1)/(2 * p2*p2) ) + p3 * exp( -(dx-p4)*(dx-p4)/(2 * p5 * p5) ) );

  return prob; 
}

// --------------------------------------------------------- 

