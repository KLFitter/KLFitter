/*!
 * \class KLFitter::ResDoubleGaussPt
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

#ifndef RESDOUBLEGAUSSPT
#define RESDOUBLEGAUSSPT

#include <vector>
#include "ResolutionBase.h" 

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class ResDoubleGaussPt : public ResolutionBase
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    ResDoubleGaussPt(const char * filename); 

    /**
     * A constructor.
     * @param parameters The parameters of the parameterization. 
     */
    ResDoubleGaussPt(std::vector<double> const& parameters);

    /**
     * The default destructor.
     */
    virtual ~ResDoubleGaussPt(); 

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

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

  private: 

  }; 
        
} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

