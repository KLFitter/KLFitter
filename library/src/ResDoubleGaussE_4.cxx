#include "ResDoubleGaussE_4.h" 
#include <iostream>

#include <cmath>
//used for lJets/bjets in mc12
// --------------------------------------------------------- 
KLFitter::ResDoubleGaussE_4::ResDoubleGaussE_4(const char * filename) : KLFitter::ResDoubleGaussBase(filename)
{ }

// ---------------------------------------------------------
KLFitter::ResDoubleGaussE_4::ResDoubleGaussE_4(std::vector<double> const& parameters) : KLFitter::ResDoubleGaussBase(parameters)
{ }

// --------------------------------------------------------- 
KLFitter::ResDoubleGaussE_4::~ResDoubleGaussE_4()
{
  ;
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_4::GetMean1(double x){
  return fParameters[0] + fParameters[1] / x; 
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_4::GetSigma1(double x){
  return fParameters[2] + fParameters[3] / sqrt(x);
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_4::GetAmplitude2(double x){
  return fParameters[4] + fParameters[5] / x; 
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_4::GetMean2(double x){
  return fParameters[6] + fParameters[7]  / sqrt(x); 
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_4::GetSigma2(double x){
  return fParameters[8] + fParameters[9] * x; 
}

// --------------------------------------------------------- 
