#include "ResDoubleGaussPt.h" 
#include <iostream>
#include <cmath>

// --------------------------------------------------------- 
KLFitter::ResDoubleGaussPt::ResDoubleGaussPt(const char * filename) : KLFitter::ResDoubleGaussBase(filename)
{ }

// ---------------------------------------------------------
KLFitter::ResDoubleGaussPt::ResDoubleGaussPt(std::vector<double> const& parameters) : KLFitter::ResDoubleGaussBase(parameters)
{ }

// --------------------------------------------------------- 
KLFitter::ResDoubleGaussPt::~ResDoubleGaussPt()
{
  ;
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussPt::GetMean1(double x){
  return fParameters[0] + x * fParameters[1];
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussPt::GetSigma1(double x){
  return fParameters[2] + x * fParameters[3];
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussPt::GetAmplitude2(double x){
  return fParameters[4] + x * fParameters[5];
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussPt::GetMean2(double x){
  return fParameters[6] + x * fParameters[7];
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGaussPt::GetSigma2(double x){
  return fParameters[8] + x * fParameters[9];
}

// --------------------------------------------------------- 
