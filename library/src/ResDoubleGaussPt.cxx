#include "ResDoubleGaussPt.h" 
#include "ResDoubleGaussE_1.h"
#include <iostream>
#include <cmath>

// --------------------------------------------------------- 
KLFitter::ResDoubleGaussPt::ResDoubleGaussPt(const char * filename) : KLFitter::ResolutionBase(10)
{
  // read parameters from file 
  ReadParameters(filename, 10); 
}

// ---------------------------------------------------------
KLFitter::ResDoubleGaussPt::ResDoubleGaussPt(std::vector<double> const& parameters) : KLFitter::ResolutionBase(parameters)
{
  // check number of parameters 
  if (parameters.size() != 10)
    {
      std::cout << "KLFitter::ResDoubleGaussPt::ResDoubleGaussPt(). Number of parameters != 10." << std::endl;
      return;
    } 
}

// --------------------------------------------------------- 
KLFitter::ResDoubleGaussPt::~ResDoubleGaussPt()
{
  ;
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussPt::GetSigma(double xmeas){
  /* Calculate mean width of both gaussians; weight the width of the 2nd one with its amplitude */
  double sigma1 = fParameters[2] + xmeas * fParameters[3];
  double sigma2 = fParameters[8] + xmeas * fParameters[9];
  double amplitude2 = fParameters[4] + xmeas * fParameters[5];
  double sigma = (sigma1 + amplitude2*sigma2) / (1+amplitude2);

  /* sigma estimates the fractional resolution, but we want absolute */
  return sigma*xmeas;
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussPt::p(double x, double xmeas, bool &good)
{
  double p1 = fParameters[0] + x * fParameters[1];
  double p2 = fParameters[2] + x * fParameters[3];
  double p3 = fParameters[4] + x * fParameters[5];
  double p4 = fParameters[6] + x * fParameters[7];
  double p5 = fParameters[8] + x * fParameters[9];

  // sanity checks for p2, p3 and p5
  good = ResDoubleGaussE_1::CheckDoubleGaussianSanity(p2, p3, p5);

  double dx = (x - xmeas) / x; 

  // calculate double-Gaussian 
  double prob= 1./sqrt(2.*M_PI) / (p2 + p3 * p5) * ( exp( -(dx - p1)*(dx - p1)/(2 * p2*p2) ) + p3 * exp( -(dx - p4)*(dx -p4)/(2 * p5 * p5) ) );

  return prob; 
}

// --------------------------------------------------------- 

