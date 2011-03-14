/*!
 * \class KLFitter::LikelihoodTTGamma_HadWRadDecay
 * \brief A class implementing a likelihood for radiative decay of the hadronic W
 * \author Johannes Erdmann
 * \version 1.3
 * \date 26.02.2010
 *
 * This class represents a likelihood for radiative decay of the hadronic W
 */

// --------------------------------------------------------- 

#ifndef LIKELIHOODTTGAMMA_HADWRADDECAY
#define LIKELIHOODTTGAMMA_HADWRADDECAY

// --------------------------------------------------------- 

#include "LikelihoodTTGamma.h" 

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class LikelihoodTTGamma_HadWRadDecay : public KLFitter::LikelihoodTTGamma
  {
                
  public: 
                
    /** 
     * The default constructor. 
     */ 
    LikelihoodTTGamma_HadWRadDecay(); 
                
    /**
     * The default destructor.
     */
    virtual ~LikelihoodTTGamma_HadWRadDecay(); 

  protected: 

    /**
     * Update 4-vectors of model particles. 
     * @return An error flag.
     */ 
    virtual int CalculateLorentzVectors(std::vector <double> parameters); 

  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 
