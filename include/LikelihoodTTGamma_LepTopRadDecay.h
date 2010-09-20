/*!
 * \class KLFitter::LikelihoodTTGamma_LepTopRadDecay
 * \brief A class implementing a likelihood for radiative decay of the leptonic top
 * \author Johannes Erdmann
 * \version 1.3
 * \date 26.02.2010
 *
 * This class represents a likelihood for radiative decay of the leptonic top
 */

// --------------------------------------------------------- 

#ifndef LIKELIHOODTTGAMMA_LEPTOPRADDECAY
#define LIKELIHOODTTGAMMA_LEPTOPRADDECAY

// --------------------------------------------------------- 

#include "LikelihoodTTGamma.h" 

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class LikelihoodTTGamma_LepTopRadDecay : public KLFitter::LikelihoodTTGamma
  {
                
  public: 
                
    /** 
     * The default constructor. 
     */ 
    LikelihoodTTGamma_LepTopRadDecay(); 
                
    /**
     * The default destructor.
     */
    virtual ~LikelihoodTTGamma_LepTopRadDecay(); 

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
