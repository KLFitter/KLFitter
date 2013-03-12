#include "ResDoubleGaussE_5.h" 
#include <iostream>

#include <cmath>
//used for electrons in mc12
// --------------------------------------------------------- 
KLFitter::ResDoubleGaussE_5::ResDoubleGaussE_5(const char * filename) : KLFitter::ResDoubleGaussBase(filename)
{ }

// ---------------------------------------------------------
KLFitter::ResDoubleGaussE_5::ResDoubleGaussE_5(std::vector<double> const& parameters) : KLFitter::ResDoubleGaussBase(parameters)
{ }

// --------------------------------------------------------- 
KLFitter::ResDoubleGaussE_5::~ResDoubleGaussE_5()
{
  ;
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_5::GetMean1(double x){
  return fParameters[0]  + fParameters[1] * x; 
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_5::GetSigma1(double x){
  return fParameters[2]  + fParameters[3]/sqrt(x);
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_5::GetAmplitude2(double x){
  return fParameters[4] + fParameters[5] * x; 
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_5::GetMean2(double x){
  return fParameters[6] + fParameters[7]/sqrt(x); 
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussE_5::GetSigma2(double x){
  return fParameters[8] + fParameters[9] * x; 
}

// --------------------------------------------------------- 
