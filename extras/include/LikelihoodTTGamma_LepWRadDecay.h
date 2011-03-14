/*!
 * \class KLFitter::LikelihoodTTGamma_LepWRadDecay
 * \brief A class implementing a likelihood for radiative decay of the leptonic W
 * \author Johannes Erdmann
 * \version 1.3
 * \date 26.02.2010
 *
 * This class represents a likelihood for radiative decay of the leptonic W
 */

// --------------------------------------------------------- 

#ifndef LIKELIHOODTTGAMMA_LEPWRADDECAY
#define LIKELIHOODTTGAMMA_LEPWRADDECAY

// --------------------------------------------------------- 

#include "LikelihoodTTGamma.h" 

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class LikelihoodTTGamma_LepWRadDecay : public KLFitter::LikelihoodTTGamma
  {
                
  public: 
                
    /** 
     * The default constructor. 
     */ 
    LikelihoodTTGamma_LepWRadDecay(); 
                
    /**
     * The default destructor.
     */
    virtual ~LikelihoodTTGamma_LepWRadDecay(); 

  protected: 

    /**
     * Update 4-vectors of model particles. 
     * @return An error flag.
     */ 
    virtual int CalculateLorentzVectors(std::vector <double> parameters); 

    /**
     * Calculates the neutrino pz solutions from the measured values
     * and the W mass.
     * @return A vector with 0, 1 or 2 neutrino pz solutions.
     */
    virtual std::vector<double> GetNeutrinoPzSolutions();

  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 
