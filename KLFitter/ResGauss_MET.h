/*!
 * \class KLFitter::ResGauss_MET
 * \brief A class describing a Gaussian resolution, parametrized for MET. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.4
 * \date 24.06.2011
 *
 * This class offers a simple parameterization of a resolution. The
 * parameterization is a Gaussian with a width of a constant times the
 * square root of the true parameter.
 */

// --------------------------------------------------------- 

#ifndef RESGAUSS_MET
#define RESGAUSS_MET

#include "ResolutionBase.h" 

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class ResGauss_MET : public ResolutionBase
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    ResGauss_MET(const char * filename); 

    /** 
     * A constructor. 
     * @param parameters The parameters of the parameterization. 
     */ 
    ResGauss_MET(std::vector<double> const& parameters); 
                
    /**
     * The default destructor.
     */
    virtual ~ResGauss_MET(); 

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /**
     * Return the width of the TF depending on the value of sumET.
     * Used to adjust the range of the fit parameter that correspond to the TF.
     * @param sumet SumET as parameter for the MET TF.
     * @return The width. 
     */ 
    virtual double GetSigma(double sumet);

    /**
     * Return the probability of the true value of x given the
     * measured value, xmeas.
     * @param x The true value of x.
     * @param xmeas The measured value of x.
     * @param good False if problem with TF.
     * @param sumet SumET, as the width of the TF depends on this.
     * @return The probability. 
     */ 
    double p(double x, double xmeas, bool &good, double sumet);

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

  private: 

  }; 
        
} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

