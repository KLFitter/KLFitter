#include "ResGauss_MET.h" 
#include <TMath.h>
#include <iostream>

// --------------------------------------------------------- 
KLFitter::ResGauss_MET::ResGauss_MET(const char * filename) : KLFitter::ResolutionBase(4)
{
  // read parameters from file 
  ReadParameters(filename, 4); 
}

// --------------------------------------------------------- 
KLFitter::ResGauss_MET::ResGauss_MET(std::vector<double> const& parameters) :KLFitter::ResolutionBase(parameters)
{
  // check number of parameters 
  if (parameters.size() != 4)
    {
      std::cout << "KLFitter::ResGauss_MET::ResGauss_MET(). Number of parameters != 4." << std::endl;
      return;
    } 
}

// --------------------------------------------------------- 
KLFitter::ResGauss_MET::~ResGauss_MET()
{

}

// --------------------------------------------------------- 
double KLFitter::ResGauss_MET::GetSigma(double sumet){
  return fParameters[0]+fParameters[1]/(1+exp(-fParameters[2]*(sumet-fParameters[3])));
}

// --------------------------------------------------------- 
double KLFitter::ResGauss_MET::p(double x, double xmeas, bool &good, double sumet)
{
  good = true;
	// calculate MET TF with 4 parameters (MC10b or later)
	double sigma = GetSigma(sumet);
	return TMath::Gaus(xmeas, x, sigma, true); 
}

// --------------------------------------------------------- 
