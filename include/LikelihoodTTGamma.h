/*!
 * \class KLFitter::LikelihoodTTGamma_RadTopProd
 * \brief A class implementing a base likelihood for the ttbar+gamma events (semileptonic channel)
 * \author Johannes Erdmann
 * \version 1.3
 * \date 26.02.2010
 *
 * This class represents a base likelihood for ttbar+gamma (semileptonic)
 */

// --------------------------------------------------------- 

#ifndef LIKELIHOODTTGAMMA
#define LIKELIHOODTTGAMMA

// --------------------------------------------------------- 

#include "PREPROC.h"
#include "LikelihoodTopLeptonJets.h" 

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class LikelihoodTTGamma : public KLFitter::LikelihoodTopLeptonJets
  {
                
  public: 
                
    /** 
     * The default constructor. 
     */ 
    LikelihoodTTGamma(); 
                
    /**
     * The default destructor.
     */
    virtual ~LikelihoodTTGamma(); 

    /**
     * Update 4-vectors of model particles. 
     * @return An error flag.
     */ 
    virtual int CalculateLorentzVectors(std::vector <double> parameters); 

    /** 
     * The posterior probability definition, overloaded from BCModel. 
     * @param parameters A vector of parameters (double values). 
     * @return The logarithm of the prior probability. 
     */
    virtual double LogLikelihood(std::vector <double> parameters);


    /**
     * Get initial values for the parameters with a dummy of "0.0" for the neutrino pz.
     * The decision on the initial value for the neutrino pz then needs to be done in
     * GetInitialParameters().
     * @return vector of initial values.
     */
    std::vector<double> GetInitialParametersWoNeutrinoPz();

  protected: 

    /**
     * Adjust parameter ranges 
     */ 
    int AdjustParameterRanges(); 

    /**
     * Define the parameters of the fit. 
     */ 
    // this is protected, because of implicit calls from base class constructor
    // -> has become just a helper for the constructor
    void DefineParameters();

    /**
     * Define the model particles
     * @return An error code.
     */ 
    // this is protected, because of implicit calls from base class constructor
    // -> has become just a helper for the constructor
    int DefineModelParticles();

    /**
     * Calculates the neutrino pz solutions from the measured values
     * and the W mass.
     * @return A vector with 0, 1 or 2 neutrino pz solutions.
     */
    virtual std::vector<double> GetNeutrinoPzSolutions() { return KLFitter::LikelihoodTopLeptonJets::GetNeutrinoPzSolutions(); };

  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 
