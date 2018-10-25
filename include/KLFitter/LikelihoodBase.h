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

#ifndef KLFITTER_LIKELIHOODBASE_H_
#define KLFITTER_LIKELIHOODBASE_H_

#include <iostream>
#include <memory>
#include <vector>

#include "BAT/BCLog.h"
#include "BAT/BCModel.h"
#include "KLFitter/ParticleCollection.h"
#include "KLFitter/PhysicsConstants.h"

// ---------------------------------------------------------

namespace KLFitter {
class Permutations;
class DetectorBase;

/**
 * A base class for likelihoods. DETAILED DESCRIPTION TO BE ADDED
 */
class LikelihoodBase : public BCModel {
 public:
  /// Enumerate for b-tagging possibilities
  enum BtaggingMethod{
    kNotag,            ///< flavor-tagging information ignored
    kVeto,             ///< b-tagged jets on light positions vetoed
    kVetoNoFit,        ///< b-tagged jets on light positions vetoed, no fit performed
    kVetoNoFitLight,   ///< non-tagged jets on b-jet positions vetoed, no fit performed
    kVetoNoFitBoth,    ///< combination of kVetoNoFit and kVetoNoFitLight
    kVetoHybridNoFit,  ///< like kVetoNoFit, unless all permutations vetoed, then kVetoNoFitLight
    kWorkingPoint,     ///< event probability considers b-tagging efficiency values
    kVetoLight,        ///< non-tagged jets on b-jet positions vetoed
    kVetoBoth          ///< combination of kVeto and kVetoLight
  };

  /**
   * The default constructor.
   * @param particles A pointer to the measured particles.
   */
  explicit LikelihoodBase(ParticleCollection ** particles = 0);

  /// The default destructor.
  virtual ~LikelihoodBase();

  /** \name Member functions (Get)  */
  /** @{ */

  /**
   * Return the table of physics constants.
   * @return A pointer to the physics constants.
   */
  KLFitter::PhysicsConstants* PhysicsConstants() { return &fPhysicsConstants; }

  /// Return the detector.
  KLFitter::DetectorBase* Detector() { return *fDetector; }

  /// Return the set of measured particles.
  KLFitter::ParticleCollection** PParticlesPermuted() { return fParticlesPermuted; }

  /// Return the set of model particles.
  KLFitter::ParticleCollection* ParticlesModel() {
    BuildModelParticles();
    return fParticlesModel.get();
  }

  /// Return the number of model particles.
  int NParticlesModel() { return static_cast<int>(fParticlesModel -> NParticles()); }

  /// Return the number of parameters.
  int NParameters() { return this -> GetNParameters(); }

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

  /// Get flag to use b-tagging or not.
  BtaggingMethod GetBTagging() { return fBTagMethod;}

  bool FlagIntegrate() { return fFlagIntegrate; }

  /** @} */
  /** \name Member functions (Set)  */
  /** @{ */

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
  virtual int SetDetector(KLFitter::DetectorBase** detector);

  /**
   * Set the measured particles.
   * @param particles The measured particles.
   * @return An error flag.
   */
  int SetParticlesPermuted(KLFitter::ParticleCollection** particles);

  /**
   * Set the values for the missing ET x and y components and the
   * SumET.
   * @param etx missing ET x component.
   * @param ety missing ET y component.
   * @param sumet total scalar ET.
   * @return An error flag.
   */
  virtual int SetET_miss_XY_SumET(double etx, double ety, double sumet) = 0;

  /**
   * Set the permutation object.
   * @param permutations The permutation object.
   * @return An error flag.
   */
  int SetPermutations(std::unique_ptr<KLFitter::Permutations>* permutations);

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
   * Set the initial values for the minimization, etc. for each
   * chain.
   * @param parameters The initial values.
   * @param nchains The number of chains.
   * @return An error flag.
   */
  int SetInitialParametersNChains(std::vector<double> const& parameters, unsigned int nchains);

  /**
   * Set which b-tagging you wish to use.
   * @param btagmethod The enum of btagging method.
   * @return An error flag.
   */
  int SetBTagging(BtaggingMethod btagmethod) { fBTagMethod = btagmethod; return 1; }

  /**
   * THIS IS AN OUTDATED METHOD - JUST HERE FOR BACKWARD COMPATIBILITY.
   * Set flag to use b-tagging or not.
   * @param flag The flag.
   * @return An error flag.
   */
  int SetFlagBTagging(bool flag) {
    std::cout << "LikelihoodBase::SetFlagBTagging(bool flag) is an outdated method - please use SetBTagging(BtaggingMethod btagmethod, double cutvalue, double btageff, double btagrej)." << std::endl;
    fBTagMethod = flag ? kVeto : kNotag;
    return 1;
  }

  /**
   * Set flag FlagIsNan. This Flag should be true if Minuit gave
   * parameters with NaN values to LogLikelihood.
   * @param flag The flag.
   * @return An error flag.
   */
  int SetFlagIsNan(bool flag) { fFlagIsNan = flag; return 1; }

  /**
   * Get flag FlagIsNan. This Flag should be true if Minuit gave
   * parameters with NaN values to LogLikelihood.
   */
  bool GetFlagIsNan(void) { return fFlagIsNan; }

  /**
   * Set flag to integrate or not.
   * @param flag The flag.
   * @return An error flag.
   */
  int SetFlagIntegrate(bool flag) { fFlagIntegrate = flag; return 1; }

  /**
   * Set flag to use measured jet masses (true) instead of
   * parton masses (false);
   */
  void SetFlagUseJetMass(bool flag) { fFlagUseJetMass = flag; }

  /** @} */
  /** \name Member functions (misc)  */
  /** @{ */

  /**
   * Update 4-vectors of model particles.
   * @return An error flag.
   */
  virtual int CalculateLorentzVectors(std::vector <double> const& parameters) = 0;

  /**
   * Initialize the likelihood for the event
   * @return An error code
   */
  virtual int Initialize();

  /// Adjust parameter ranges
  virtual int AdjustParameterRanges() = 0;

  /**
   * Define the model particles
   * @return An error code.
   */
  virtual int DefineModelParticles() = 0;

  /**
   * Propagate the b-tagging information from the permuted
   * (measured) particles to the model particles.
   */
  void PropagateBTaggingInformation();

  /**
   * Request necessary resolution functions from the detector.
   * @return An error code.
   */
  virtual void RequestResolutionFunctions() = 0;

  /** @} */
  /** \name Member functions (BAT)  */
  /** @{ */

  /// Define the parameters of the fit.
  virtual void DefineParameters() = 0;

  /**
   * The prior probability definition, overloaded from BCModel.
   * @param parameters A vector of parameters (double values).
   * @return The logarithm of the prior probability.
   */
  virtual double LogAPrioriProbability(const std::vector <double> & parameters) override { return 0; }

  /**
   * The posterior probability definition, overloaded from BCModel.
   * @param parameters A vector of parameters (double values).
   * @return The logarithm of the prior probability.
   */
  virtual double LogLikelihood(const std::vector <double> & parameters) = 0;

  /**
   * The posterior probability definition, overloaded from
   * BCModel. Split up into several subcomponents
   * @param parameters A vector of parameters (double values).
   * @return A vector with the components of the logarithm of the
   * prior probability.
   */
  virtual std::vector<double> LogLikelihoodComponents(std::vector <double> parameters) = 0;

  /**
   * Return the log of the event probability fof the current
   * combination
   * @return The event probability
   */
  virtual double LogEventProbability();

  /**
   * Return the contribution from b tagging to the log of the
   * event probability for the current combination
   * @return The event probability contribution
   */
  virtual double LogEventProbabilityBTag();

  /**
   * Remove invariant particle permutations.
   * @return An error code.
   */
  virtual int RemoveInvariantParticlePermutations() = 0;

  /**
   * Remove forbidden particle permutations.
   * @return An error code.
   */
  virtual int RemoveForbiddenParticlePermutations();

  /**
   * Build the model particles from the best fit parameters.
   * @return An error code.
   */
  virtual int BuildModelParticles() = 0;

  /**
   * Get initial values for the parameters.
   * @return vector of initial values.
   */
  virtual std::vector<double> GetInitialParameters() = 0;

  /**
   * Check if there are TF problems.
   * @return Return false if TF problem.
   */
  virtual bool NoTFProblem(std::vector<double> parameters);

  /// Best fit parameters, overloaded from BCModel
  std::vector <double> GetBestFitParameters();
  using BCModel::GetBestFitParameters;

  /// Errors of the best fit parameters, overloaded from BCModel
  std::vector <double> GetBestFitParameterErrors();
  using BCModel::GetBestFitParameterErrors;

  /// Best fit parameter at position (i), overloaded from BCModel
  double GetBestFitParameter(unsigned int index);
  using BCModel::GetBestFitParameter;

  /// Errors of the best fit parameter i, overloaded from BCModel
  double GetBestFitParameterError(unsigned int index);
  using BCModel::GetBestFitParameterError;

  /**
   * Check if the permutation is LH invariant.
   * @return Permutation of the invariant partner, -1 if there is no one.
   */
  virtual int LHInvariantPermutationPartner(int /*iperm*/, int /*nperms*/, int* /*switchpar1*/, int* /*switchpar2*/) { return -1; }

  /**
   * Write parameters from fCachedParametersVector.at(iperm) to
   * fCachedParameters.
   * @param iperm Current permutation
   * @return An error code.
   */
  int GetParametersFromCache(int iperm);

  /**
   * Write parameters to fCachedParametersVector.at(iperm) from
   * GetBestFitParameter().
   * @param iperm Current permutation
   * @param nperms Number of permutations
   * @return An error code.
   */
  int SetParametersToCache(int iperm, int nperms);

  /// Normalization factor of the probability, overloaded from BCModel
  double GetIntegral();
  using  BCIntegrate::GetIntegral;

  /**
   * Resets the values of all parameter cache vectors
   * @return An error code.
   */
  int ResetCache();

  /** @} */

 protected:
  /// Save permuted particles.
  virtual int SavePermutedParticles() = 0;

  /// Save resolution functions.
  virtual int SaveResolutionFunctions() = 0;

  /**
   * Set model parton mass according to fFlagUseJetMass.
   * @param jetmass The jet mass.
   * @param quarkmass The quark mass.
   * @param px The parton px (will be modified, if necessary).
   * @param py The parton py (will be modified, if necessary).
   * @param pz The parton pz (will be modified, if necessary).
   * @param e The parton energy (not modified).
   * @return The parton mass.
   */
  double SetPartonMass(double jetmass, double quarkmass, double *px, double *py, double *pz, double e);

  /// A pointer to the measured particles.
  KLFitter::ParticleCollection** fParticlesPermuted;

  /// A pointer to the permutation object.
  std::unique_ptr<KLFitter::Permutations>* fPermutations;

  /// A pointer to the model particles.
  std::unique_ptr<KLFitter::ParticleCollection> fParticlesModel;

  /// A pointer to the table of physics constants
  KLFitter::PhysicsConstants fPhysicsConstants;

  /// A pointer to the detector
  KLFitter::DetectorBase** fDetector;

  /// The event probabilities for the different permutations
  std::vector<double> fEventProbability;

  /// A flag to integrate over the likelihood or not
  bool fFlagIntegrate;

  /// A flag for knowing that Minuit gave parameters with NaN values to LogLikelihood
  bool fFlagIsNan;

  /// A flag for using the measured jet masses instead of parton masses
  bool fFlagUseJetMass;

  /// Global variable for TF problems.
  bool fTFgood;

  /// Name of btagging enum
  BtaggingMethod fBTagMethod;

  /// The cached parameters used for the current permutation
  std::vector<double>  fCachedParameters;

  /// The cached parameter errors used for the current permutation
  std::vector<double>  fCachedParameterErrors;

  /// A vector of cached parameters, one for each permutation. Has to be set via fitter.
  std::vector<std::vector<double> >  fCachedParametersVector;

  /// A vector of cached parameter errors, one for each permutation. Has to be set via fitter.
  std::vector<std::vector<double> >  fCachedParameterErrorsVector;

  /// The cached normalization, needed for the overloaded BCIntegrate::GetIntegral
  double  fCachedNormalization;

  /// A vector of cached parameters, one for each permutation. Has to be set via fitter.
  std::vector<double>  fCachedNormalizationVector;
};
}  // namespace KLFitter

#endif  // KLFITTER_LIKELIHOODBASE_H_
