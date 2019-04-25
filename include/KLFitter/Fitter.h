/*
 * Copyright (c) 2009--2018, the KLFitter developer team
 *
 * This file is part of KLFitter.
 *
 * KLFitter is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * KLFitter is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with KLFitter. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef KLFITTER_FITTER_H_
#define KLFITTER_FITTER_H_

#include <memory>
#include <vector>

// ---------------------------------------------------------

/**
 * The main KLFitter namespace. This namespace should contain all
 * functions, classes and sub-namespaces of the project.
 */
namespace KLFitter {
class ParticleCollection;
class DetectorBase;
class LikelihoodBase;
class PermutationHandler;

/**
 * \class KLFitter::Fitter
 * \brief The fitter class.
 *
 * This class owns all particles, the detector description, the
 * likelihood, etc. This is the class seen by the user.
 */
class Fitter final {
 public:
  /// The default constructor.
  Fitter();

  /// The (defaulted) destructor.
  ~Fitter();

  /** \name Member functions (Get)  */
  /** @{ */

  /// Return the detector.
  KLFitter::DetectorBase * Detector() { return fDetector; }

  /// Return the measured particles.
  KLFitter::ParticleCollection * Particles() { return fParticles; }

  /// Return the permutation object.
  KLFitter::PermutationHandler * Permutations() { return fPermutations.get(); }

  /// Return the lieklihood .
  KLFitter::LikelihoodBase * Likelihood() { return fLikelihood; }

  /// Return the Minuit status
  int MinuitStatus() { return fMinuitStatus; }

  /// Return the convergence status bit.
  unsigned int ConvergenceStatus() { return fConvergenceStatus; }

  /** @} */
  /** \name Member functions (Set)  */
  /** @{ */

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
  int SetParticles(KLFitter::ParticleCollection * particles);

  /**
   * Set truth particles.
   * @param particles A pointer to a set of particles.
   * @return An error flag.
   */
  int SetMyParticlesTruth(KLFitter::ParticleCollection * particles);

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

  /** @} */
  /** \name Member functions (misc)  */
  /** @{ */

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

  /// Turn off simulated annealing.
  void TurnOffSA() { fTurnOffSA = true; }

  /// Enumerator for convergence errors.
  enum ConvergenceErrors {
    MinuitDidNotConverge = 1,
    FitAbortedDueToNaN = 2,
    AtLeastOneFitParameterAtItsLimit = 3,
    InvalidTransferFunctionAtConvergence = 4
  };

  /// Bit masks for convergence errors.
  static const unsigned int MinuitDidNotConvergeMask = 0x1 << MinuitDidNotConverge;
  static const unsigned int FitAbortedDueToNaNMask = 0x1 << FitAbortedDueToNaN;
  static const unsigned int AtLeastOneFitParameterAtItsLimitMask = 0x1 << AtLeastOneFitParameterAtItsLimit;
  static const unsigned int InvalidTransferFunctionAtConvergenceMask = 0x1 << InvalidTransferFunctionAtConvergence;

  /// Enumerator for the minimization methods.
  enum kMinimizationMethod { kMinuit, kSimulatedAnnealing, kMarkovChainMC };

  /// Set the minimization method.
  void SetMinimizationMethod(kMinimizationMethod method) { fMinimizationMethod = method; }

  /**
   * Write fCachedMinuitStatus and fCachedConvergenceStatus to
   * fCachedMinuitStatusVector.at(iperm)
   * and fCachedConvergenceStatusVector.at(iperm)
   * @param iperm Current permutation
   * @param nperms The number of permutations
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

  /** @} */

 private:
  /// A pointer to the detector.
  KLFitter::DetectorBase * fDetector;

  /// A pointer to the set of original particles.
  KLFitter::ParticleCollection * fParticles;

  /// The x and y component of the missing ET and the sumET.
  double ETmiss_x;
  double ETmiss_y;
  double SumET;

  /// A pointer to the set of permuted particles.
  KLFitter::ParticleCollection * fParticlesPermuted;

  /// A pointer to the set of truth particles.
  KLFitter::ParticleCollection * fMyParticlesTruth;

  /// A pointer to the likelihood.
  KLFitter::LikelihoodBase * fLikelihood;

  /// A pointer to the permutation object.
  std::unique_ptr<KLFitter::PermutationHandler> fPermutations;

  /// The TMinuit status
  int fMinuitStatus;

  /// The convergence status bit
  unsigned int fConvergenceStatus;

  /// Flag for turning off simulated annealing.
  int fTurnOffSA;

  /// The minimization method.
  kMinimizationMethod fMinimizationMethod;

  /// A vector of cached Minuit status
  std::vector<int>  fCachedMinuitStatusVector;

  /// A vector of cached convergence status
  std::vector<unsigned int>  fCachedConvergenceStatusVector;

  /**
   * Resets the values of all parameter cache vectors
   * @return An error code.
   */
  int ResetCache();
};
}  // namespace KLFitter

#endif  // KLFITTER_FITTER_H_

