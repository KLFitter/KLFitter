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

#ifndef KLFITTER_LIKELIHOODTOPALLHADRONIC_H_
#define KLFITTER_LIKELIHOODTOPALLHADRONIC_H_

#include <iostream>
#include <vector>

#include "KLFitter/LikelihoodBase.h"

// ---------------------------------------------------------

namespace KLFitter {
class ResolutionBase;

/**
 * Likelihood for the ttbar allhadronic channel, where both tops
 * decay hadronically. This likelihood takes six jets, out of
 * which two are expected to come from b-quarks.
 */
class LikelihoodTopAllHadronic : public KLFitter::LikelihoodBase {
 public:
  /// The default constructor.
  LikelihoodTopAllHadronic();

  /// The (defaulted) destructor.
  ~LikelihoodTopAllHadronic();

  /** \name Member functions (Set)  */
  /** @{ */

  /// Enumerator for the parameters.
  enum Parameters { parBhad1E, parBhad2E, parLQ1E, parLQ2E, parLQ3E, parLQ4E, parTopM };

  /**
   * Set the values for the missing ET x and y components and the SumET.
   * Reimplemented with dummy implementation to overwrite purely virtual
   * function in base class.
   * @return An error flag.
   */
  int SetET_miss_XY_SumET(double /*etx*/, double /*ety*/, double /*sumet*/) override { return 1; }

  /// Request the necessary resolution functions from the detector.
  void RequestResolutionFunctions() override;

  /**
   * Set a flag. If flag is true the invariant top quark mass is
   * fixed to the pole mass.
   * @param flag The flag.
   */
  void SetFlagTopMassFixed(bool flag) { fFlagTopMassFixed = flag; }

  void SetFlagGetParSigmasFromTFs(bool flag) { fFlagGetParSigmasFromTFs = flag; }

  /** @} */
  /** \name Member functions (BAT)  */
  /** @{ */

  /// Define the parameters of the fit.
  void DefineParameters() override;

  /**
   * The posterior probability definition, overloaded from BCModel.
   * @param parameters A vector of parameters (double values).
   * @return The logarithm of the prior probability.
   */
  double LogLikelihood(const std::vector <double> &  parameters) override;

  /**
   * The posterior probability definition, overloaded from BCModel. Split up into several subcomponents
   * @param parameters A vector of parameters (double values).
   * @return A vector with the components of the logarithm of the prior probability. Its components are:
   * 0:  TF_bhad1
   * 1:  TF_bhad2
   * 2:  TF_lq1
   * 3:  TF_lq2
   * 4:  TF_lq3
   * 5:  TF_lq4
   * 6:  BW_Whad1
   * 7:  BW_Whad2
   * 8:  BW_Thad1
   * 9:  BW_Thad2
   */
  std::vector<double> LogLikelihoodComponents(std::vector <double> parameters) override;

  /// Get initial values for the parameters.
  std::vector<double> GetInitialParameters() override;

  /** @} */

 protected:
  /** \name Member functions (misc)  */
  /** @{ */

  /**
   * Update 4-vectors of model particles.
   * @return An error flag.
   */
  int CalculateLorentzVectors(std::vector <double> const& parameters) override;

  /// Adjust parameter ranges
  int AdjustParameterRanges() override;

  /**
   * Define the model particles
   * @return An error code.
   */
  int DefineModelParticles() override;

  /**
   * Remove invariant particle permutations.
   * @return An error code.
   */
  int RemoveInvariantParticlePermutations() override;

  /**
   * Build the model particles from the best fit parameters.
   * @return An error code.
   */
  int BuildModelParticles() override;

  /// Save permuted particles.
  int SavePermutedParticles() override;

  /// Save resolution functions.
  int SaveResolutionFunctions() override;

  /** @} */
  /** \name Member attributes */
  /** @{ */

  /// A flag for using a fixed top mass (true) or not (false).
  bool fFlagTopMassFixed;

  ///  Flag for using ResolutionBase::GetSigma() to retrieve the parameter ranges
  bool fFlagGetParSigmasFromTFs;

  /// Pointer to resolution function for hadronic b quark.
  ResolutionBase * fResEnergyBhad1;

  /// Pointer to resolution function for leptonic b quark.
  ResolutionBase * fResEnergyBhad2;

  /// Pointer to resolution function for first light quark jet.
  ResolutionBase * fResEnergyLQ1;

  /// Pointer to resolution function for second light quark jet.
  ResolutionBase * fResEnergyLQ2;

  /// Pointer to resolution function for third light quark jet.
  ResolutionBase * fResEnergyLQ3;

  /// Pointer to resolution function for fourth light quark jet.
  ResolutionBase * fResEnergyLQ4;

  /** @} */
  /** \name Member attributes (measured parameters) */
  /** @{ */

  double bhad1_meas_e;
  double bhad1_meas_p;
  double bhad1_meas_m;
  double bhad1_meas_deteta;
  double bhad1_meas_eta;
  double bhad1_meas_phi;
  double bhad1_meas_px;
  double bhad1_meas_py;
  double bhad1_meas_pz;

  double bhad2_meas_e;
  double bhad2_meas_p;
  double bhad2_meas_m;
  double bhad2_meas_deteta;
  double bhad2_meas_eta;
  double bhad2_meas_phi;
  double bhad2_meas_px;
  double bhad2_meas_py;
  double bhad2_meas_pz;

  double lq1_meas_e;
  double lq1_meas_p;
  double lq1_meas_m;
  double lq1_meas_deteta;
  double lq1_meas_eta;
  double lq1_meas_phi;
  double lq1_meas_px;
  double lq1_meas_py;
  double lq1_meas_pz;

  double lq2_meas_e;
  double lq2_meas_p;
  double lq2_meas_m;
  double lq2_meas_deteta;
  double lq2_meas_eta;
  double lq2_meas_phi;
  double lq2_meas_px;
  double lq2_meas_py;
  double lq2_meas_pz;

  double lq3_meas_e;
  double lq3_meas_p;
  double lq3_meas_m;
  double lq3_meas_deteta;
  double lq3_meas_eta;
  double lq3_meas_phi;
  double lq3_meas_px;
  double lq3_meas_py;
  double lq3_meas_pz;

  double lq4_meas_e;
  double lq4_meas_p;
  double lq4_meas_m;
  double lq4_meas_deteta;
  double lq4_meas_eta;
  double lq4_meas_phi;
  double lq4_meas_px;
  double lq4_meas_py;
  double lq4_meas_pz;

  /** @} */
  /** \name Member attributes (fitted parameters) */
  /** @{ */

  double bhad1_fit_e;
  double bhad1_fit_px;
  double bhad1_fit_py;
  double bhad1_fit_pz;

  double bhad2_fit_e;
  double bhad2_fit_px;
  double bhad2_fit_py;
  double bhad2_fit_pz;

  double lq1_fit_e;
  double lq1_fit_px;
  double lq1_fit_py;
  double lq1_fit_pz;

  double lq2_fit_e;
  double lq2_fit_px;
  double lq2_fit_py;
  double lq2_fit_pz;

  double lq3_fit_e;
  double lq3_fit_px;
  double lq3_fit_py;
  double lq3_fit_pz;

  double lq4_fit_e;
  double lq4_fit_px;
  double lq4_fit_py;
  double lq4_fit_pz;

  double whad1_fit_m;
  double whad2_fit_m;
  double thad1_fit_m;
  double thad2_fit_m;

  /** @} */
};
}  // namespace KLFitter

#endif  // KLFITTER_LIKELIHOODTOPALLHADRONIC_H_
