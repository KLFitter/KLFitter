/*!
 * \class KLFitter::Fitter
 * \brief The fitter class. 
 * \author Kevin Kr&ouml;ninger
 * \version 1.3
 * \date 03.12.2009
 *
 * This class owns all particles, the detector description, the
 * likelihood, etc. This is the class seen by the user. 
 */

// --------------------------------------------------------- 

#ifndef FITTER
#define FITTER

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

  class Particles;
  class DetectorBase;
  class LikelihoodBase;
  class Permutations;

  class Fitter
  {
                
  public: 
                
    /** \name Constructors and destructors */ 
    /* @{ */ 
                
    /** 
     * The default constructor. 
     */ 
    Fitter(); 
                
    /**
     * The default destructor.
     */
    virtual ~Fitter(); 

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /**
     * Return the detector.
     * @return A pointer to the detector.
     */ 
    KLFitter::DetectorBase * Detector()
      { return fDetector; }; 

    /**
     * Return the measured particles. 
     * @return A pointer to the particles.
     */ 
    KLFitter::Particles * Particles()
      { return fParticles; }; 
    KLFitter::Particles ** PParticles()
      { return &fParticles; }; 

    /**
     * Return the permutation object. 
     * @return A pointer to the permutation object. 
     **/ 
    KLFitter::Permutations * Permutations()
      { return fPermutations; }; 
    KLFitter::Permutations ** PPermutations()
      { return &fPermutations; }; 

    /**
     * Return the lieklihood . 
     * @return A pointer to the likelihood object. 
     **/ 
    KLFitter::LikelihoodBase * Likelihood()
      { return fLikelihood; }; 
    KLFitter::LikelihoodBase ** PLikelihood()
      { return &fLikelihood; }; 

    /**
     * Return the Minuit status
     * @return The Minuit stats
     */ 
    int MinuitStatus()
    { return fMinuitStatus; }; 

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */

    /**
     * Set the detector description. 
     * @param detector A pointer to the detector. 
     * @return An error code. 
     */ 
    int SetDetector(KLFitter::DetectorBase * detector); 

    /**
     * Set the particles. 
     * @param particles A pointer to a set of particles. 
     * @return An error flag. 
     */ 
    int SetParticles(KLFitter::Particles * particles); 

    /**
     * Set x and y component of the missing ET.
     * @param etx component of the missing ET.
     * @param ety component of the missing ET.
     * @return An error flag.
     */
    int SetET_miss_XY(double etx, double ety);

    /**
     * Set the likelihood for the actual fit. 
     * @param likelihood A pointer to the likelihood. 
     * @return An error code. 
     */ 
    int SetLikelihood(KLFitter::LikelihoodBase * likelihood); 

    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */

    /**
     * Perform the fit for a single permutation of jets and leptons. 
     * @param index The permutation index.
     * @return An error code.
     */ 
    int Fit(int index); 

    /**
     * Perform the fit for all permutations of jets and leptons. 
     * @return An error code. 
     */ 
    int Fit(); 

    /**
     * Check if detector, likelihood, etc. are all set. 
     * @return A status code. 
     */ 
    int Status(); 

    /**
     * Turn of simulated annealing.
     */
    void TurnOffSA() { fTurnOffSA = true; };

    /* @} */
                
  private: 

    /**
     * A pointer to the detector. 
     */ 
    KLFitter::DetectorBase * fDetector; 

    /**
     * A pointer to the set of original particles. 
     */ 
    KLFitter::Particles * fParticles; 

    /**
     * The x and y component of the missing ET.
     */
    double ETmiss_x;
    double ETmiss_y;

    /**
     * A pointer to the set of permuted particles. 
     */ 
    KLFitter::Particles * fParticlesPermuted; 

    /**
     * A pointer to the likelihood. 
     */ 
    KLFitter::LikelihoodBase * fLikelihood; 

    /**
     * A pointer to the permutation object.
     */ 
    KLFitter::Permutations * fPermutations; 

    /**
     * The TMinuit status
     */ 
    int fMinuitStatus; 

    /**
     * Flag for turning off simulated annealing.
     */
    int fTurnOffSA;

  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

