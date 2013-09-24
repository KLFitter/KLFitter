/*!
 * \class KLFitter::LikelihoodTopLeptonJets_JetAngles
 * \brief A class implementing a likelihood for the ttbar lepton+jets channel. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.3
 * \date 03.12.2009
 *
 * This class represents a likelihood for the ttbar into lepton+jets.
 */

// --------------------------------------------------------- 

#ifndef LIKELIHOODTOPLEPTONJETS_JETANGLES
#define LIKELIHOODTOPLEPTONJETS_JETANGLES

// --------------------------------------------------------- 

#include "PREPROC.h"
#include "LikelihoodBase.h" 
#include "TLorentzVector.h"

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class LikelihoodTopLeptonJets_JetAngles : public KLFitter::LikelihoodBase
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    LikelihoodTopLeptonJets_JetAngles(); 
                
    /**
     * The default destructor.
     */
    virtual ~LikelihoodTopLeptonJets_JetAngles(); 

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */

    /**
     * Enumerator for the lepton type.
     */
    enum LeptonType { kElectron, kMuon };

    /**
     * Set the values for the missing ET x and y components.
     * @param etx missing ET x component.
     * @param ety missing ET y component.
     * @return An error flag.
     */
    int SetET_miss_XY(double etx, double ety);

    /**
     * Set a flag. If flag is true the invariant top quark mass is
     * fixed to the pole mass. 
     * @param flag The flag. 
     */ 
    void SetFlagTopMassFixed(bool flag)
    { fFlagTopMassFixed = flag; }; 

    void SetFlagUseJetMass(bool flag)
    { fFlagUseJetMass = flag; }; 

    /**
     * Set the type of lepton 
     * @param leptontype The type of lepton: kElectron or kMuon
     */ 
    void SetLeptonType(LeptonType leptontype); 

    /**
     * Set the type of lepton 
     * @param leptontype The type of lepton: electron(1) or muon (2) 
     */ 
    void SetLeptonType(int leptontype); 

    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */

    /* @} */
    /** \name Member functions (BAT)  */
    /* @{ */

    /**
     * Define the parameters of the fit. 
     */ 
    virtual void DefineParameters();

    /** 
     * The prior probability definition, overloaded from BCModel. 
     * @param parameters A vector of parameters (double values). 
     * @return The logarithm of the prior probability. 
     */
    virtual double LogAPrioriProbability(const std::vector <double> & KLFITTER_UNUSED(parameters))
    { return 0; }; 
                
    /** 
     * The posterior probability definition, overloaded from BCModel. 
     * @param parameters A vector of parameters (double values). 
     * @return The logarithm of the prior probability. 
     */
    virtual double LogLikelihood(const std::vector <double> & parameters); 

    /**
     * Get initial values for the parameters.
     * @return vector of initial values.
     */
    virtual std::vector<double> GetInitialParameters();

    /**
     * Get initial values for the parameters with a dummy of "0.0" for the neutrino pz.
     * The decision on the initial value for the neutrino pz then needs to be done in
     * GetInitialParameters().
     * @return vector of initial values.
     */
    virtual std::vector<double> GetInitialParametersWoNeutrinoPz();

    /**
     * Check if there are TF problems.
     * @return Return false if TF problem.
     */
    virtual bool NoTFProblem(std::vector<double> parameters);

    /* @} */

  protected: 

    /** \name Member functions (misc)  */
    /* @{ */

    /**
     * Update 4-vectors of model particles. 
     * @return An error flag.
     */ 
    virtual int CalculateLorentzVectors(std::vector <double> const& parameters); 

    /**
     * Initialize the likelihood for the event
     */ 
    virtual int Initialize(); 

    /**
     * Adjust parameter ranges 
     */ 
    virtual int AdjustParameterRanges(); 

    /**
     * Define the model particles
     * @return An error code.
     */ 
    virtual int DefineModelParticles();

    /**
     * Remove invariant particle permutations.
     * @return An error code. 
     */ 
    int RemoveInvariantParticlePermutations(); 
    
     /**
     * Remove forbidden particle permutations.
     * @return An error code. 
     */ 
    int RemoveForbiddenParticlePermutations();    

    /* @} */

  protected: 
                
    /**
     * A flag for using a fixed top mass (true) or not (false).
     */ 
    bool fFlagTopMassFixed; 

    /**
     * A flag for using the measured jet masses (true) instead of
     * parton masses (false);
     */ 
    bool fFlagUseJetMass; 

    /**
     * Return the neutrino pz solutions from the measured values
     * and the W mass.
     * @return A vector with 0, 1 or 2 neutrino pz solutions.
     */
    virtual std::vector<double> GetNeutrinoPzSolutions();

    /**
     * Calculates the neutrino pz solutions from the measured values
     * and the W mass. An additional particle to be added to the 
     * charged lepton may be specified, for example a photon
     * in ttbargamma, which is radiated from the leptonic W
     * or the charged lepton;
     * @param Pointer to a 4-vector of a particle which is added to the charged lepton in the calculation
     * @return A vector with 0, 1 or 2 neutrino pz solutions.
     */
    std::vector<double> CalculateNeutrinoPzSolutions(TLorentzVector * additionalParticle = 0x0);

    /**
     * The values of the x component of the missing ET.
     */
    double ETmiss_x;

    /**
     * The values of the y component of the missing ET.
     */
    double ETmiss_y;

    /**
     * An index deciding if the event is electron (1) or muon (2) plus
     * jets.
     */ 
    LeptonType fTypeLepton; 

    /**
     * The value of pi stored here for CPU time reasons.
     */ 
    double fPi;

    /**
     * The value of 2*pi stored here for CPU time reasons.
     */ 
    double fTwoPi;

    /**
     * Global variable for TF problems.
     */
    bool fTFgood;

  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

