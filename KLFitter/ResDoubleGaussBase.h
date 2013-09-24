/*!
 * \class KLFitter::ResDoubleGaussBase
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

#ifndef RESDOUBLEGAUSSBASE
#define RESDOUBLEGAUSSBASE

#include <vector>
#include "ResolutionBase.h" 
#include <iostream>

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class ResDoubleGaussBase : public ResolutionBase
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    ResDoubleGaussBase(const char * filename); 

    /**
     * A constructor.
     * @param parameters The parameters of the parameterization. 
     */
    ResDoubleGaussBase(std::vector<double> const& parameters);

    /**
     * The default destructor.
     */
    virtual ~ResDoubleGaussBase(); 

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /**
     * Calculate the mean of the first Gaussian from the TF parameters and the value of x.
     * @param x The value of x.
     * @return The width. 
     */ 
    virtual double GetMean1(double x) = 0;

    /**
     * Calculate the width of the first Gaussian from the TF parameters and the value of x.
     * @param x The value of x.
     * @return The width. 
     */ 
    virtual double GetSigma1(double x) = 0;

    /**
     * Calculate the amplitude of the second Gaussian from the TF parameters and the value of x.
     * @param x The value of x.
     * @return The width. 
     */ 
    virtual double GetAmplitude2(double x) = 0;

    /**
     * Calculate the mean of the second Gaussian from the TF parameters and the value of x.
     * @param x The value of x.
     * @return The width. 
     */ 
    virtual double GetMean2(double x) = 0;

    /**
     * Calculate the width of the sedcond Gaussian from the TF parameters and the value of x.
     * @param x The value of x.
     * @return The width. 
     */ 
    virtual double GetSigma2(double x) = 0;

    /**
     * Return the approximate width of the TF depending on the measured value of x.
     * Used to adjust the range of the fit parameter that correspond to the TF.
     * @param xmeas The measured value of x.
     * @return The width. 
     */ 
    virtual double GetSigma(double xmeas);

    /**
     * Return the probability of the true value of x given the
     * measured value, xmeas.
     * @param x The true value of x.
     * @param xmeas The measured value of x.
     * @param good False if problem with TF.
     * @return The probability. 
     */ 
    double p(double x, double xmeas, bool &good); 

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */
        
    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */
                
    /* @} */

    /**
     * Sanity check for double gaussian parameters p2, p3 and p5 (1st sigma, scale and 2nd sigma).
     * @param p2 (the 1st sigma).
     * @param p3 (the scale parameter).
     * @param p5 (the 2nd sigma).
     * @return False if problem with TF.
     */
    inline static bool CheckDoubleGaussianSanity(double &sigma1, double &amplitude2, double &sigma2) {
      if (amplitude2 < 0.) amplitude2 = 0.;      
      if (sigma1 < 0.) {
//        std::cout << "KLFitter::ResDoubleGauss::CheckDoubleGaussianSanity() ERROR IN TRANSFERFUNCTIONS the sigma of the 1st Gaussian is < 0  -  FIT RESULT MAY NOT BE RELIABLE" << std::endl;
//        std::cout << "--> Fitter is in a bad condition! Please check your input files (E out of validation scope?)." << std::endl;
        sigma1 = 0.00000001;
        return false;
	}
      if (sigma2 < 0.) {
//        std::cout << "KLFitter::ResDoubleGauss::CheckDoubleGaussianSanity() ERROR IN TRANSFERFUNCTIONS the sigma of the 2nd Gaussian is < 0  -  FIT RESULT MAY NOT BE RELIABLE" << std::endl;
//        std::cout << "--> Fitter is in a bad condition! Please check your input files (E out of validation scope?)." << std::endl;
        sigma2 = 0.000000001;
        return false;
	}

      return true;
    }

  private: 

  }; 
        
} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

