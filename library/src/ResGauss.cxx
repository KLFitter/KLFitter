#include "ResGauss.h" 
#include <TMath.h>
#include <iostream>

// --------------------------------------------------------- 
KLFitter::ResGauss::ResGauss(const char * filename) : KLFitter::ResolutionBase(1)
{
  // read parameters from file 
  ReadParameters(filename, 1); 
}

// --------------------------------------------------------- 
KLFitter::ResGauss::ResGauss(double sigma) : KLFitter::ResolutionBase(1)
{
  // set parameter
  SetPar(0, sigma);
}

// --------------------------------------------------------- 
KLFitter::ResGauss::~ResGauss()
{

}

// --------------------------------------------------------- 
double KLFitter::ResGauss::GetSigma(double dummy){
  return fParameters[0];
}

// --------------------------------------------------------- 
double KLFitter::ResGauss::p(double x, double xmeas, bool &good)
{
  good = true;
  return TMath::Gaus(xmeas, x, fParameters[0], true); 
}

// --------------------------------------------------------- 

