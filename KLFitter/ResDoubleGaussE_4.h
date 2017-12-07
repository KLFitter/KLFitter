/*!
 * \class KLFitter::ResDoubleGaussE_4
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

#ifndef RESDOUBLEGAUSSE_4
#define RESDOUBLEGAUSSE_4

#include <vector>
#include "ResDoubleGaussBase.h" 

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class ResDoubleGaussE_4 : public ResDoubleGaussBase
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    ResDoubleGaussE_4(const char * filename); 

    /**
     * A constructor.
     * @param parameters The parameters of the parameterization. 
     */
    ResDoubleGaussE_4(std::vector<double> const& parameters);

    /**
     * The default destructor.
     */
    virtual ~ResDoubleGaussE_4(); 

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /**
     * Calculate the mean of the first Gaussian from the TF parameters and the value of x.
     * @param x The value of x.
     * @return The width. 
     */ 
    virtual double GetMean1(double x);

    /**
     * Calculate the width of the first Gaussian from the TF parameters and the value of x.
     * @param x The value of x.
     * @return The width. 
     */ 
    virtual double GetSigma1(double x);

    /**
     * Calculate the amplitude of the second Gaussian from the TF parameters and the value of x.
     * @param x The value of x.
     * @return The width. 
     */ 
    virtual double GetAmplitude2(double x);

    /**
     * Calculate the mean of the second Gaussian from the TF parameters and the value of x.
     * @param x The value of x.
     * @return The width. 
     */ 
    virtual double GetMean2(double x);

    /**
     * Calculate the width of the sedcond Gaussian from the TF parameters and the value of x.
     * @param x The value of x.
     * @return The width. 
     */ 
    virtual double GetSigma2(double x);

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */
        
    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */
                
    /* @} */

  private: 

  }; 
        
} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

