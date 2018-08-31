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

#ifndef KLFITTER_LIKELIHOODSGTOPWTLJ_H_
#define KLFITTER_LIKELIHOODSGTOPWTLJ_H_

#include <vector>

#include "KLFitter/LikelihoodBase.h"

// ---------------------------------------------------------

namespace KLFitter {
class ResolutionBase;

/**
 * Likelihood for the production of a single top quark in
 * association with a W boson (tW). The top quark is assumed to
 * decay into jets, whereas the W is assumed to decay into lepton
 * and neutrino.
 */
class LikelihoodSgTopWtLJ : public KLFitter::LikelihoodBase {
 public:
  /// The default constructor.
  LikelihoodSgTopWtLJ();

  /// The (defaulted) destructor.
  ~LikelihoodSgTopWtLJ();

  /** \name Member functions (Get)  */
  /** @{ */

  /// Get a pointer to the four-vector of the lepton
  TLorentzVector* GetLepton(KLFitter::Particles* particles);

  /// Return the lepton type
  int GetLeptonType();

  /// Return the top decay hypothesis
  bool GetHadronicTop() { return fHadronicTop; }

  /** @} */
  /** \name Member functions (Set)  */
  /** @{ */

  /// Enumerator for the lepton type.
  enum LeptonType { kElectron, kMuon };

  /// Enumerator for the parameters.
  enum Parameters { parBE, parLQ1E, parLQ2E, parLepE, parNuPx, parNuPy, parNuPz };

  /**
   * Set the values for the missing ET x and y components and the SumET.
   * @param etx missing ET x component.
   * @param ety missing ET y component.
   * @param sumet total scalar ET.
   * @return An error flag.
   */
  int SetET_miss_XY_SumET(double etx, double ety, double sumet) override;

  /// Request the necessary resolution functions from the detector.
  void RequestResolutionFunctions() override;

  /// Associate the hadronic leg of the event to the top quark for likelihood calculation.
  void SetHadronicTop() { fHadronicTop = true; }

  /// Associate the leptonic leg of the event to the top quark for likelihood calculation.
  void SetLeptonicTop() { fHadronicTop = false; }

  /// Set the type of lepton
  void SetLeptonType(int leptontype);

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
  double LogLikelihood(const std::vector<double> & parameters) override;

  /// Get initial values for the parameters.
  std::vector<double> GetInitialParameters() override;

  /// Dummy reimplementation from the abstract base class.
  std::vector<double> LogLikelihoodComponents(std::vector<double>) override { return std::vector<double>{}; }

  /** @} */

 protected:
  /** \name Member functions (misc)  */
  /** @{ */

  /**
   * Update 4-vectors of model particles.
   * @return An error flag.
   */
  int CalculateLorentzVectors(const std::vector <double>& parameters) override;

  /// Initialize the likelihood for the event
  int Initialize() override;

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

  /** @} */

  /**
   * A flag for associating either the hadronic (true)
   * or the leponic (false) leg of the event to the top quark
   * for likelihood calculation.
   */
  bool fHadronicTop;

  /**
   * Calculates the neutrino pz solutions from the measured values
   * and the W mass.
   * @return A vector with 0, 1 or 2 neutrino pz solutions.
   */
  std::vector<double> GetNeutrinoPzSolutions();

  /// Save permuted particles.
  int SavePermutedParticles() override;

  /// Save resolution functions.
  int SaveResolutionFunctions() override;

  /// The values of the x component of the missing ET.
  double ETmiss_x;

  /// The values of the y component of the missing ET.
  double ETmiss_y;

  /// The values of the total scalar ET.
  double SumET;

  /// Index whether l+jets event is electron (1) or muon (2).
  int fTypeLepton;

  /// Save resolution functions since the eta of the partons is not fitted.
  ResolutionBase * fResEnergyB;
  ResolutionBase * fResEnergyLQ1;
  ResolutionBase * fResEnergyLQ2;
  ResolutionBase * fResLepton;
  ResolutionBase * fResMET;

  /// Save measured particle values for frequent calls
  double b_meas_e;
  double b_meas_p;
  double b_meas_m;
  double b_meas_deteta;
  double b_meas_eta;
  double b_meas_phi;
  double b_meas_px;
  double b_meas_py;
  double b_meas_pz;

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

  double lep_meas_e;
  double lep_meas_deteta;
  double lep_meas_sintheta;
  double lep_meas_pt;
  double lep_meas_px;
  double lep_meas_py;
  double lep_meas_pz;

  /// Save fit particle values for frequent calls
  double b_fit_e;
  double b_fit_px;
  double b_fit_py;
  double b_fit_pz;

  double lq1_fit_e;
  double lq1_fit_px;
  double lq1_fit_py;
  double lq1_fit_pz;

  double lq2_fit_e;
  double lq2_fit_px;
  double lq2_fit_py;
  double lq2_fit_pz;

  double lep_fit_e;
  double lep_fit_px;
  double lep_fit_py;
  double lep_fit_pz;

  double nu_fit_e;
  double nu_fit_px;
  double nu_fit_py;
  double nu_fit_pz;

  double whad_fit_m;
  double wlep_fit_m;
  double thad_fit_m;
  double tlep_fit_m;
};
}  // namespace KLFitter

#endif  // KLFITTER_LIKELIHOODSGTOPWTLJ_H_
