/*!
 * \class KLFitter::LikelihoodBase
 * \brief A base class for likelihoods. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.3
 * \date 03.12.2009
 */

// --------------------------------------------------------- 

#ifndef LIKELIHOODBASE
#define LIKELIHOODBASE

// --------------------------------------------------------- 

#include "PREPROC.h"
#include "Particles.h" 

#include "BAT/BCModel.h"

#include <TLorentzVector.h>

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class PhysicsConstants;
  class Permutations;
  class DetectorBase;  

  class LikelihoodBase : public BCModel 
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /**
     * The default constructor.
     * @param particles A pointer to the measured particles. 
     */ 
    LikelihoodBase(Particles ** particles = 0); 
                
    /**
     * The default destructor.
     */
    virtual ~LikelihoodBase(); 

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /**
     * Return the table of physics constants.
     * @return A pointer to the physics constants.
     */ 
    KLFitter::PhysicsConstants* PhysicsConstants()
      { return fPhysicsConstants; }; 

    /**
     * Return the detector.
     * @return A pointer to the detector.
     */ 
    KLFitter::DetectorBase* Detector()
      { return *fDetector; }; 

    /**
     * Return the set of measured particles. 
     * @return A pointer to the particles. 
     */ 
    KLFitter::Particles** PParticlesPermuted()
      { return fParticlesPermuted; }; 

    /**
     * Return the set of model particles. 
     * @return A pointer to the particles. 
     */ 
    KLFitter::Particles* ParticlesModel()
      { return fParticlesModel; }; 
    KLFitter::Particles** PParticlesModel()
      { return &fParticlesModel; }; 

    /**
     * Return model particles at some index. 
     * @param index The particle index.
     * @return A pointer to the TLorentzVector of the particle. 
     */ 
    TLorentzVector* ParticleModel(int index); 

    /**
     * Return model particles with a name.
     * @param name The particle name. 
     * @return A pointer to the TLorentzVector of the particle. 
     */                 
    TLorentzVector* ParticleModel(const char * name); 

    /**
     * Return the number of model particles. 
     * @return The number of model particles.
     */ 
    int NParticlesModel()
    { return int(fParticlesModel -> NParticles()); }; 

    /**
     * Return the number of parameters. 
     * @return The number of parameters of the model. 
     */ 
    int NParameters()
    { return this -> GetNParameters(); };

    /**
     * Return the lower boundary of a parameter
     * @param index The index of the parameter.
     * @return The lower boundary. 
     */ 
    double ParMin(int index); 

    /**
     * Return the upper boundary of a parameter
     * @param index The index of the parameter.
     * @return The upper boundary. 
     */ 
    double ParMax(int index); 

    double CutBTag()
    { return fCutBTag; }; 

    bool FlagBTagging() { return fFlagBTagging; } 

    bool FlagIntegrate() { return fFlagIntegrate; } 

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */

    /**
     * Set the physics constants
     * @param physicsconstants A pointer to physics constants.
     * @return An error flag
     */
    int SetPhysicsConstants(KLFitter::PhysicsConstants* physicsconstants); 

    /**
     * Set the detector
     * @param detector A pointer to a pointer of the detector. 
     * @return An error flag
     */
    int SetDetector(KLFitter::DetectorBase** detector); 

    /**
     * Set the measured particles. 
     * @param particles The measured particles. 
     * @return An error flag. 
     */ 
    int SetParticlesPermuted(KLFitter::Particles** particles); 

    /**
     * Set the values for the missing ET x and y components.
     * @param etx missing ET x component.
     * @param ety missing ET y component.
     * @return An error flag.
     */
    virtual int SetET_miss_XY(double KLFITTER_UNUSED(etx), double KLFITTER_UNUSED(ety)) { return 0; }

    /**
     * Set the permutation object.
     * @param permutations The permutation object.
     * @return An error flag. 
     */ 
    int SetPermutations(KLFitter::Permutations ** permutations); 

    /**
     * Set the range of a model parameter. 
     * @param index The index of the parameter.
     * @param parmin The minimum value.
     * @param parmax The maximum value. 
     */ 
    int SetParameterRange(int index, double parmin, double parmax); 

    /**
     * Set the initial values for the minimization, etc. 
     * @param parameters The initial values.
     * @return An error flag. 
     */ 
    int SetInitialParameters(std::vector<double> const& parameters); 

    /**
     * Set flag to use b-tagging or not.
     * @param flag The flag.
     * @return An error flag. 
     */ 
    int SetFlagBTagging(bool flag) { fFlagBTagging = flag; return 1; } 

    /**
     * Set flag FlagIsNan. This Flag should be true if Minuit gave parameters with NaN values to LogLikelihood.
     * @param flag The flag.
     * @return An error flag. 
     */ 
    int SetFlagIsNan(bool flag) { fFlagIsNan=flag; return 1; }

    /**
     * Get flag FlagIsNan. This Flag should be true if Minuit gave parameters with NaN values to LogLikelihood.
     * @return The flag.
     */ 
    bool GetFlagIsNan(void) { return fFlagIsNan; }

    void SetCutBTagging(double cut)
    { fCutBTag = cut; }; 

    /**
     * Set flag to integrate or not.
     * @param flag The flag. 
     * @return An error flag.
     */ 
    int SetFlagIntegrate(bool flag) { fFlagIntegrate = flag; return 1; } 

    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */

    /**
     * Initialize the likelihood for the event
     * @return An error code 
     */ 
    virtual int Initialize()
    { return 1; }; 

    /**
     * Sets the flavor tag based on the b-tagging probability.
     */ 
    void CalculateFlavorTags(); 

    /* @} */
    /** \name Member functions (BAT)  */
    /* @{ */

    /**
     * Define the parameters of the fit. 
     */ 
    virtual void DefineParameters()
    { ; };

    /** 
     * The prior probability definition, overloaded from BCModel. 
     * @param parameters A vector of parameters (double values). 
     * @return The logarithm of the prior probability. 
     */
    virtual double LogAPrioriProbability(std::vector <double> KLFITTER_UNUSED(parameters))
    { return 0; } 
                
    /** 
     * The posterior probability definition, overloaded from BCModel. 
     * @param parameters A vector of parameters (double values). 
     * @return The logarithm of the prior probability. 
     */
    virtual double LogLikelihood(std::vector <double> KLFITTER_UNUSED(parameters))
    { return 0; } 

    /**
     * Return the log of the event probability fof the current
     * combination
     * @return The event probability 
     */ 
    double LogEventProbability(); 

    /**
     * Return the combined b-tagging probability for all jets
     * @return The probability. 
     */ 
    double BTaggingProbability();

    /**
     * Remove invariant particle permutations.
     * @return An error code. 
     */ 
    virtual int RemoveInvariantParticlePermutations()
    { return 1; }; 

    /**
     * Get initial values for the parameters.
     * @return vector of initial values.
     */
    virtual std::vector<double> GetInitialParameters()
    { std::vector<double> v; return v; };

    /* @} */

  protected: 

    /**
     * A pointer to the measured particles. 
     */ 
    KLFitter::Particles** fParticlesPermuted; 

    /**
     * A pointer to the permutation object. 
     */ 
    KLFitter::Permutations** fPermutations; 

    /**
     * A pointer to the model particles. 
     */ 
    KLFitter::Particles* fParticlesModel; 

    /**
     * A pointer to the table of physics constants 
     */ 
    KLFitter::PhysicsConstants* fPhysicsConstants; 

    /**
     * A pointer to the detector 
     */ 
    KLFitter::DetectorBase** fDetector; 

    /**
     * The event probabilities for the different permutations 
     */ 
    std::vector<double> fEventProbability; 

    /**
     * A flag to use b-tagging or not 
     */ 
    bool fFlagBTagging; 

    /**
     * Cut value for b-tagging. 
     */ 
    double fCutBTag; 

    /**
     * A flag to integrate over the likelihood or not 
     */ 
    bool fFlagIntegrate; 

    /**
     * A flag for knowing that Minuit gave parameters with NaN values to LogLikelihood
     */ 
    bool fFlagIsNan;

  private: 

  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

