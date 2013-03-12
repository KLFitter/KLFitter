/*!
 * \class KLFitter::ResolutionBase
 * \brief A base class for describing resolutions. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.3
 * \date 03.12.2009
 *
 * This base class can be used to decribe resolutions. 
 */

// --------------------------------------------------------- 

#ifndef RESOLUTIONBASE
#define RESOLUTIONBASE

// --------------------------------------------------------- 

#include "PREPROC.h"
#include <vector>

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class ResolutionBase
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     * @param npar The number of parameters needed for the parameterization.
     */ 
    ResolutionBase(int npar); 

    /** 
     * A constructor. 
     * @param parameters A vector with the parameters. 
     */ 
    ResolutionBase(std::vector<double> parameters); 
                
    /**
     * The default destructor.
     */
    virtual ~ResolutionBase(); 

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /**
     * Return the (approximate) width of the TF depending on the relevant parameter 
     * (mostly xmeas, but sumET in case of ResGauss_MET).
     * Use to adjust the range of the fit parameter that correspond to the TF.
     * @param par Parameter on which the width depends
     * @return The width. 
     */ 
    virtual double GetSigma(double KLFITTER_UNUSED(par)) = 0;

    /**
     * Return the probability of the true value of x given the
     * measured value, xmeas.
     * @param x The true value of x.
     * @param xmeas The measured value of x.
     * @param good False if problem with TF.
     * @return The probability. 
     */ 
    virtual double p(double KLFITTER_UNUSED(x), double KLFITTER_UNUSED(xmeas), bool &good)
    { good = true; return 0; } 

    /**
     * Return the probability of the true value of x given the
     * measured value, xmeas.
     * @param x The true value of x.
     * @param xmeas The measured value of x.
     * @param good False if problem with TF.
     * @param par Optional additional parameter (SumET in case of MET TF).
     * @return The probability. 
     */ 
    virtual double p(double KLFITTER_UNUSED(x), double KLFITTER_UNUSED(xmeas), bool &good, double KLFITTER_UNUSED(par))
    { good = true; return 0; } 

    /**
     * Return a parameter of the parameterization. 
     * @param index The parameter index. 
     * @param par The parameter value. 
     * @return An error flag. 
     */ 
    int Par(int index, double &par); 

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */

    /**
     * Set a parameter value. 
     * @param index The parameter index.
     * @param value The parameter value. 
     * @return An error code. 
     */ 
    int SetPar(int index, double value); 

    /**
     * Set all parameter values.
     * @param parameters A vector of parameters.
     * @return An error code. 
     */ 
    int SetPar(std::vector<double> parameters); 

    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */

    /**
     * Read parameter values from ASCII file. 
     * @param filename The name of the file. 
     * @param nparameters The number of parameters. 
     * @return An error code. 
     */ 
    int ReadParameters(const char * filename, int nparameters); 

    /**
     * Return a status code. 
     * @return A status code (1: ok, 0: error). 
     */ 
    int Status(); 
                
    /* @} */

  protected: 

    /**
     * The number of parameters. 
     */ 
    int fNParameters; 

    /**
     * The parameter values. 
     */ 
    std::vector <double> fParameters; 

  }; 
        
} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

