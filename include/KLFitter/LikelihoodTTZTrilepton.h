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

#ifndef KLFITTER_LIKELIHOODTTZTRILEPTON_H_
#define KLFITTER_LIKELIHOODTTZTRILEPTON_H_

#include <iostream>
#include <vector>

#include "KLFitter/LikelihoodBase.h"

// ---------------------------------------------------------

namespace KLFitter {
namespace Resolutions {
  class ResolutionBase;
}

/**
 * Likelihood for the ttZ trilepton channel, where ttbar decays
 * into lepton+jets and the Z boson into a lepton pair. It is
 * largely based on LikelihoodTopLeptonJets, with the necessary
 * adjustments to include the extra Z boson.
 */
class LikelihoodTTZTrilepton : public KLFitter::LikelihoodBase {
 public:
  /// The default constructor.
  LikelihoodTTZTrilepton();

  /// The (defaulted) destructor.
  ~LikelihoodTTZTrilepton();

  /** \name Member functions (Get)  */
  /** @{ */

  /// Get the cut-off value of the 1/E^2 distribution.
  double GetInvMassCutoff() { return fInvMassCutoff; }

  /// Get the fraction of on-shell events.
  float GetOnShellFraction() { return fOnShellFraction; }

  /** @} */
  /** \name Member functions (Set)  */
  /** @{ */

  /// Enumerator for the lepton type.
  enum LeptonType { kElectron, kMuon };

  /// Enumerator for the parameters.
  enum Parameters { parBhadE, parBlepE, parLQ1E, parLQ2E, parLepE, parNuPx, parNuPy, parNuPz, parTopM, parLepZ1E, parLepZ2E, parZM};

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

  /// Set the cut-off value of the 1/E^2 distribution.
  void SetInvMassCutoff(double cutoff) { fInvMassCutoff = cutoff; }

  /// Set the fraction of on-shell events.
  void SetOnShellFraction(double fraction) { fOnShellFraction = fraction; }

  /**
   * Set a flag. If flag is true the invariant top quark mass is
   * fixed to the pole mass.
   * @param flag The flag.
   */
  void SetFlagTopMassFixed(bool flag) { fFlagTopMassFixed = flag; }

  void SetFlagGetParSigmasFromTFs(bool flag) { fFlagGetParSigmasFromTFs = flag; }

  /// Set the type of lepton according to LeptonType.
  void SetLeptonType(LeptonType leptontype);

  /// Set the type of lepton: (1) electron, or (2) muon.
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
  double LogLikelihood(const std::vector <double> & parameters) override;

  /**
   * The posterior probability definition, overloaded from BCModel. Split up into several subcomponents
   * @param parameters A vector of parameters (double values).
   * @return A vector with the components of the logarithm of the prior probability. Its components are:
   * 0:  TF_bhad
   * 1:  TF_blep
   * 2:  TF_lq1
   * 3:  TF_lq2
   * 4:  TF_lepZ1
   * 5:  TF_lepZ2
   * 6:  TF_lep
   * 7:  TF_METx
   * 8:  TF_METy
   * 9:  BW_Whad
   * 10: BW_Wlep
   * 11: BW_Thad
   * 12: BW_Tlep
   * 13: BW_Z
   */
  std::vector<double> LogLikelihoodComponents(std::vector <double> parameters) override;

  /// Get initial values for the parameters.
  std::vector<double> GetInitialParameters() override;

  /**
   * Get initial values for the parameters with a dummy of "0.0" for the neutrino pz.
   * The decision on the initial value for the neutrino pz then needs to be done in
   * GetInitialParameters().
   * @return vector of initial values.
   */
  std::vector<double> GetInitialParametersWoNeutrinoPz();

  /** @} */

 protected:
  /** \name Member functions (misc)  */
  /** @{ */

  /**
   * Provide a local modification of BCMath::LogBreitWignerRel such
   * that the relativistic Breit-Wigner distribution is normalised
   * to 1. The function then returns the log of this distribution.
   *
   * @param x Value to be evaluated.
   * @param mean The mean of the distribution, i.e. Z pole mass.
   * @param gamma The FWHM of the distribution, i.e. the Z decay width.
   * @return Log of the relativistic B-W.
   */
  double LogBreitWignerRelNorm(const double& x, const double& mean, const double& gamma);

  /**
   * Evaluate a combined Z/y invariant mass distribution. The B-W
   * function and the 1/E^2 distribution are weighted according to
   * fOnShellFraction.
   *
   * @param x Value to be evaluated.
   * @param mean The mean of the distribution, i.e. Z pole mass.
   * @param gamma The FWHM of the distribution, i.e. the Z decay width.
   * @return Log of combined mass distribution.
   */
  double LogZCombinedDistribution(const double& x, const double& mean, const double& gamma);

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

  /** @} */

  /**
   * Return the neutrino pz solutions from the measured values
   * and the W mass.
   * @return A vector with 0, 1 or 2 neutrino pz solutions.
   */
  std::vector<double> GetNeutrinoPzSolutions();

  /**
   * Calculates the neutrino pz solutions from the measured values
   * and the W mass. An additional particle to be added to the
   * charged lepton may be specified, for example a photon
   * in ttbargamma, which is radiated from the leptonic W
   * or the charged lepton;
   * @param additionalParticle Pointer to a 4-vector of a particle which is
   * added to the charged lepton in the calculation
   * @return A vector with 0, 1 or 2 neutrino pz solutions.
   */
  std::vector<double> CalculateNeutrinoPzSolutions(TLorentzVector * additionalParticle = 0x0);

  /// Save permuted particles.
  int SavePermutedParticles() override;

  /// Save resolution functions.
  int SaveResolutionFunctions() override;

  /** \name Member attributes */
  /** @{ */

  /// A flag for using a fixed top mass (true) or not (false).
  bool fFlagTopMassFixed;

  ///  Flag for using Resolutions::ResolutionBase::GetSigma() to retrieve the parameter ranges
  bool fFlagGetParSigmasFromTFs;

  /// The values of the x component of the missing ET.
  double ETmiss_x;

  /// The values of the y component of the missing ET.
  double ETmiss_y;

  /// The values of the total scalar ET.
  double SumET;

  /// Index whether l+jets event is electron (1) or muon (2).
  LeptonType fTypeLepton;

  /// Cut-off value for the 1/E^2 distribution (in GeV).
  double fInvMassCutoff;

  /**
   * Fraction of on-shell events, i.e. weighting factor between
   * on-shell and off-shell distribution. A value of 1 corresponds
   * to a pure on-shell distribution, 0 to a pure off-shell
   * distribution.
   */
  double fOnShellFraction;

  /// Pointer to resolution function for hadronic b quark.
  Resolutions::ResolutionBase * fResEnergyBhad;

  /// Pointer to resolution function for leptonic b quark.
  Resolutions::ResolutionBase * fResEnergyBlep;

  /// Pointer to resolution function for first light quark jet.
  Resolutions::ResolutionBase * fResEnergyLQ1;

  /// Pointer to resolution function for second light quark jet.
  Resolutions::ResolutionBase * fResEnergyLQ2;

  /// Pointer to resolution function for the first lepton from Z.
  Resolutions::ResolutionBase * fResLeptonZ1;

  /// Pointer to resolution function for the second lepton from Z.
  Resolutions::ResolutionBase * fResLeptonZ2;

  /// Pointer to resolution function for the lepton.
  Resolutions::ResolutionBase * fResLepton;

  /// Pointer to resolution function for MET.
  Resolutions::ResolutionBase * fResMET;

  /** @{ */
  /** \name Member attributes (measured parameters) */
  /** @{ */

  double bhad_meas_e;
  double bhad_meas_p;
  double bhad_meas_m;
  double bhad_meas_deteta;
  double bhad_meas_eta;
  double bhad_meas_phi;
  double bhad_meas_px;
  double bhad_meas_py;
  double bhad_meas_pz;

  double blep_meas_e;
  double blep_meas_p;
  double blep_meas_m;
  double blep_meas_deteta;
  double blep_meas_eta;
  double blep_meas_phi;
  double blep_meas_px;
  double blep_meas_py;
  double blep_meas_pz;

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

  double lepZ1_meas_e;
  double lepZ1_meas_deteta;
  double lepZ1_meas_sintheta;
  double lepZ1_meas_pt;
  double lepZ1_meas_px;
  double lepZ1_meas_py;
  double lepZ1_meas_pz;

  double lepZ2_meas_e;
  double lepZ2_meas_deteta;
  double lepZ2_meas_sintheta;
  double lepZ2_meas_pt;
  double lepZ2_meas_px;
  double lepZ2_meas_py;
  double lepZ2_meas_pz;

  /** @} */
  /** \name Member attributes (fitted parameters) */
  /** @{ */

  double bhad_fit_e;
  double bhad_fit_px;
  double bhad_fit_py;
  double bhad_fit_pz;

  double blep_fit_e;
  double blep_fit_px;
  double blep_fit_py;
  double blep_fit_pz;

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

  double lepZ1_fit_e;
  double lepZ1_fit_px;
  double lepZ1_fit_py;
  double lepZ1_fit_pz;

  double lepZ2_fit_e;
  double lepZ2_fit_px;
  double lepZ2_fit_py;
  double lepZ2_fit_pz;

  double Z_fit_m;

  /** @} */
};
}  // namespace KLFitter

#endif  // KLFITTER_LIKELIHOODTTZTRILEPTON_H_
