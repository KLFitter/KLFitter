/*!
 * \class KLFitter::LikelihoodTTGamma_RadTopProd
 * \brief A class implementing a likelihood for radiative decay of the hadronic top
 * \author Johannes Erdmann
 * \version 1.3
 * \date 26.02.2010
 *
 * This class represents a likelihood for radiative decay of the hadronic top
 */

// --------------------------------------------------------- 

#ifndef LIKELIHOODTTGAMMA_HADTOPRADDECAY
#define LIKELIHOODTTGAMMA_HADTOPRADDECAY

// --------------------------------------------------------- 

#include "LikelihoodTTGamma.h" 

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class LikelihoodTTGamma_HadTopRadDecay : public KLFitter::LikelihoodTTGamma
  {
                
  public: 
                
    /** 
     * The default constructor. 
     */ 
    LikelihoodTTGamma_HadTopRadDecay(); 
                
    /**
     * The default destructor.
     */
    virtual ~LikelihoodTTGamma_HadTopRadDecay(); 

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
