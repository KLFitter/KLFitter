/*!
 * \class KLFitter::ResDoubleGauss
 * \brief A class describing a resolution parameterized with a double Gaussian. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.3
 * \date 03.12.2009
 *
 * This class offers a simple parameterization of a resolution. The
 * parameterization is a double Gaussian with energy dependent
 * parameters.
 */

// --------------------------------------------------------- 

#ifndef __RESDOUBLEGAUSS__H
#define __RESDOUBLEGAUSS__H

#include <vector>
#include "ResolutionBase.h" 
#include "TMath.h"

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

	class ResDoubleGauss : public ResolutionBase
	{
		
	public: 
		
		/** \name Constructors and destructors */ 
		/* @{ */ 
		
		/** 
		 * The default constructor. 
		 */ 
		ResDoubleGauss(const char * filename); 

		/**
		 * A constructor.
		 * @param parameters The parameters of the parameterization. 
		 */
		ResDoubleGauss(const std::vector<double> parameters);

		/**
		 * The default destructor.
		 */
		virtual ~ResDoubleGauss(); 

		/* @} */
		/** \name Member functions (Get)  */
		/* @{ */

		/**
		 * Return the probability of the true value of x given the
		 * measured value, xmeas.
		 * @param x The true value of x.
		 * @param xmeas The measured value of x.
		 * @return The probability. 
		 */ 
		double p(double x, double xmeas); 

		/* @} */
		/** \name Member functions (Set)  */
		/* @{ */

		/**
		 * Set the width of the Gaussian 
		 * @param sigma The width of the Gaussian. 
		 */ 
		void SetSigma(double sigma)
		{ if (sigma < 0) sigma = - sigma; this -> SetPar(0, sigma); }; 

		/* @} */
		/** \name Member functions (misc)  */
		/* @{ */
		
		/* @} */

 private: 

	}; 
	
} // namespace KLFitter 

// --------------------------------------------------------- 

#endif // __RESDOUBLEGAUSS__H

