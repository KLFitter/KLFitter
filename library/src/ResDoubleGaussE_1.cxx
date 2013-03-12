#include "ResDoubleGaussE_1.h" 
#include <iostream>

#include <cmath>

// --------------------------------------------------------- 
KLFitter::ResDoubleGaussE_1::ResDoubleGaussE_1(const char * filename) : KLFitter::ResDoubleGaussBase(filename)
{ }

// ---------------------------------------------------------
KLFitter::ResDoubleGaussE_1::ResDoubleGaussE_1(std::vector<double> const& parameters) : KLFitter::ResDoubleGaussBase(parameters)
{ }

// --------------------------------------------------------- 
KLFitter::ResDoubleGaussE_1::~ResDoubleGaussE_1()
{
  ;
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_1::GetMean1(double x){
  return fParameters[0] + fParameters[1] * x; 
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_1::GetSigma1(double x){
  return fParameters[2] / sqrt(x) + fParameters[3];
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_1::GetAmplitude2(double x){
  return fParameters[4] + fParameters[5] * x; 
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_1::GetMean2(double x){
  return fParameters[6] + fParameters[7] * x; 
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_1::GetSigma2(double x){
  return fParameters[8] + fParameters[9] * x; 
}

// --------------------------------------------------------- 
