#include "ResDoubleGaussE_4.h" 
#include "ResDoubleGaussE_1.h" 
#include <iostream>

#include <cmath>
//used for lJets/bjets in mc12
// --------------------------------------------------------- 
KLFitter::ResDoubleGaussE_4::ResDoubleGaussE_4(const char * filename) : KLFitter::ResolutionBase(10)
{
  // read parameters from file 
  ReadParameters(filename, 10); 
}

// ---------------------------------------------------------
KLFitter::ResDoubleGaussE_4::ResDoubleGaussE_4(std::vector<double> const& parameters) : KLFitter::ResolutionBase(parameters)
{
  // check number of parameters 
  if (parameters.size() != 10)
    {
      std::cout << "KLFitter::ResDoubleGaussE_4::ResDoubleGaussE_4(). Number of parameters != 10." << std::endl;
      return;
    } 
}

// --------------------------------------------------------- 
KLFitter::ResDoubleGaussE_4::~ResDoubleGaussE_4()
{
  ;
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_4::p(double x, double xmeas, bool &good)
{
  double p1  = fParameters[0] + fParameters[1] / x; 
  double p2 = fParameters[2]  + fParameters[3] / sqrt(x);
  double p3 = fParameters[4] + fParameters[5] / x; 
  double p4  = fParameters[6]  + fParameters[7]  / sqrt(x); 
  double p5 = fParameters[8] + fParameters[9] * x; 

  // sanity checks for p2, p3 and p5
  good = ResDoubleGaussE_1::CheckDoubleGaussianSanity(p2, p3, p5);

  double dx = (x - xmeas) / x; 
 
  // calculate double-Gaussian 
  return 1./sqrt(2.*M_PI) / (p2 + p3 * p5) * ( exp( -(dx-p1)*(dx-p1)/(2 * p2*p2) ) + p3 * exp( -(dx-p4)*(dx-p4)/(2 * p5 * p5) ) ); 
}

// --------------------------------------------------------- 

