#include "ResDoubleGauss.h" 

// --------------------------------------------------------- 
KLFitter::ResDoubleGauss::ResDoubleGauss(const char * filename) : KLFitter::ResolutionBase(10)
{
	// read parameters from file 
	this -> ReadParameters(filename, 10); 
}

// ---------------------------------------------------------
KLFitter::ResDoubleGauss::ResDoubleGauss(const std::vector<double> parameters) : KLFitter::ResolutionBase(parameters)
{
	// check number of parameters 
	if (parameters.size() != 10)
		{
			std::cout << "KLFitter::ResDoubleGauss::ResDoubleGauss(). Number of parameters != 10." << std::endl;
			return;
		} 
}

// --------------------------------------------------------- 
KLFitter::ResDoubleGauss::~ResDoubleGauss()
{
	;
}

// --------------------------------------------------------- 
double KLFitter::ResDoubleGauss::p(double x, double xmeas)
{
	double mean1  = fParameters.at(0) / x + fParameters.at(1); 
	double sigma1 = sqrt( fParameters.at(2)*fParameters.at(2) / x + fParameters.at(3)*fParameters.at(3) );
	double mixing = fParameters.at(4) + fParameters.at(5) * x; 
	double mean2  = fParameters.at(6) / x + fParameters.at(7); 
	double sigma2 = fParameters.at(8)/ x + fParameters.at(9); 

	double dx = (x - xmeas) / x; 

	// calculate double-Gaussian 
	double prob = 0.398942 / (sigma1 + mixing * sigma2) * 
	(exp(- (dx - mean1) * (dx - mean1) / 2.0 / sigma1 / sigma1)
	 + mixing * exp(- (dx - mean2) * (dx - mean2) / 2.0 / sigma2 / sigma2)); 

	return prob; 
}

// --------------------------------------------------------- 

