#include "ResDoubleGaussE_2.h" 
#include <iostream>

#include <cmath>

// --------------------------------------------------------- 
KLFitter::ResDoubleGaussE_2::ResDoubleGaussE_2(const char * filename) : KLFitter::ResDoubleGaussBase(filename)
{ }

// ---------------------------------------------------------
KLFitter::ResDoubleGaussE_2::ResDoubleGaussE_2(std::vector<double> const& parameters) : KLFitter::ResDoubleGaussBase(parameters)
{ }

// --------------------------------------------------------- 
KLFitter::ResDoubleGaussE_2::~ResDoubleGaussE_2()
{
  ;
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_2::GetMean1(double x){
  return fParameters[0] / sqrt(x) + fParameters[1] * x; 
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_2::GetSigma1(double x){
  return fParameters[2] / sqrt(x) + fParameters[3];
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_2::GetAmplitude2(double x){
  return fParameters[4] / sqrt(x) + fParameters[5] * x; 
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_2::GetMean2(double x){
  return fParameters[6] + fParameters[7] * x; 
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_2::GetSigma2(double x){
  return fParameters[8] + fParameters[9] * x; 
}

// --------------------------------------------------------- 
