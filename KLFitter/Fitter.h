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

#include <vector>

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

    /**
     * Return the convergence status bit.
     * @return The convergence status bit.
     */
    unsigned int ConvergenceStatus()
    { return fConvergenceStatus; }

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
     * Set truth particles. 
     * @param particles A pointer to a set of particles. 
     * @return An error flag. 
     */ 
    int SetMyParticlesTruth(KLFitter::Particles * particles); 

    /**
     * Set x and y component of the missing ET.
     * @param etx component of the missing ET.
     * @param ety component of the missing ET.
     * @param sumet The measured scalar sum of transverse energy. 
     * @return An error flag.
     */
    int SetET_miss_XY_SumET(double etx, double ety, double sumet);

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

    /**
     * Enumerator for convergence errors.
     */
    enum ConvergenceErrors {
      MinuitDidNotConverge = 1,
      FitAbortedDueToNaN = 2,
      AtLeastOneFitParameterAtItsLimit = 3,
      InvalidTransferFunctionAtConvergence = 4
    };

    /**
     * Bit masks for convergence errors.
     */
    static const unsigned int MinuitDidNotConvergeMask = 0x1 << MinuitDidNotConverge;
    static const unsigned int FitAbortedDueToNaNMask = 0x1 << FitAbortedDueToNaN;
    static const unsigned int AtLeastOneFitParameterAtItsLimitMask = 0x1 << AtLeastOneFitParameterAtItsLimit;
    static const unsigned int InvalidTransferFunctionAtConvergenceMask = 0x1 << InvalidTransferFunctionAtConvergence;

    /**
     * Enumerator for the minimization methods.
     */
    enum kMinimizationMethod { kMinuit, kSimulatedAnnealing, kMarkovChainMC };

    /**
     * Set the minimization method.
     * @param method The minimization method.
     */
    void SetMinimizationMethod(kMinimizationMethod method) { fMinimizationMethod = method; }

    /**
     * Write fCachedMinuitStatus and fCachedConvergenceStatus to 
     * fCachedMinuitStatusVector.at(iperm) 
     * and fCachedConvergenceStatusVector.at(iperm)
     * @param iperm Current permutation
     * @return An error code. 
     */
   int SetFitStatusToCache(int iperm, int nperms);

    /**
     * Write parameters from fCachedMinuitStatusVector.at(iperm) 
     * and fCachedConvergenceStatusVector.at(iperm) to fCachedMinuitStatus
     * and fCachedConvergenceStatus
     * @param iperm Current permutation
     * @return An error code. 
     */
   int GetFitStatusFromCache(int iperm);

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
     * The x and y component of the missing ET and the sumET.
     */
    double ETmiss_x;
    double ETmiss_y;
    double SumET;

    /**
     * A pointer to the set of permuted particles. 
     */ 
    KLFitter::Particles * fParticlesPermuted; 

    /**
     * A pointer to the set of truth particles. 
     */ 
    KLFitter::Particles * fMyParticlesTruth;

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
     * The convergence status bit
     */ 
    unsigned int fConvergenceStatus;

    /**
     * Flag for turning off simulated annealing.
     */
    int fTurnOffSA;

    /**
     * The minimization method.
     */
    kMinimizationMethod fMinimizationMethod;

    /**
     * A vector of cached Minuit status
     */ 
    std::vector<int>  fCachedMinuitStatusVector; 

    /**
     * A vector of cached convergence status
     */ 
    std::vector<unsigned int>  fCachedConvergenceStatusVector; 

    /**
     * Resets the values of all parameter cache vectors
     * @return An error code. 
     */
     int ResetCache();

  }; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif 

