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

#ifndef KLFITTER_BOOSTEDLIKELIHOODTOPLEPTONJETS_H_
#define KLFITTER_BOOSTEDLIKELIHOODTOPLEPTONJETS_H_

#include <iostream>
#include <vector>

class TLorentzVector;

namespace KLFitter {
  class ResolutionBase;
}

#include "KLFitter/LikelihoodBase.h"

// ---------------------------------------------------------

namespace KLFitter {
/**
 * Likelihood to reconstruct an l+jets ttbar pair in the boosted
 * top quark scenario. ADD DETAILED DESCRIPTION HERE
 */
class BoostedLikelihoodTopLeptonJets : public KLFitter::LikelihoodBase {
 public:
  /// The default constructor.
  BoostedLikelihoodTopLeptonJets();

  /// The (defaulted) destructor.
  ~BoostedLikelihoodTopLeptonJets();

  /** \name Member functions (Set)  */
  /** @{ */

  /// Enumerator for the lepton type.
  enum LeptonType { kElectron, kMuon };

  /// Enumerator for the parameters.
  enum Parameters { parBhadE, parBlepE, parLQE, parLepE, parNuPx, parNuPy, parNuPz, parTopM };

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
   * 2:  TF_lq
   * 3:  TF_lep
   * 4:  TF_METx
   * 5:  TF_METy
   * 6:  BW_Wlep
   * 7:  BW_Thad
   * 8: BW_Tlep
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

  ///  Flag for using ResolutionBase::GetSigma() to retrieve the parameter ranges
  bool fFlagGetParSigmasFromTFs;

  /// The values of the x component of the missing ET.
  double ETmiss_x;

  /// The values of the y component of the missing ET.
  double ETmiss_y;

  /// The values of the total scalar ET.
  double SumET;

  /// Index whether l+jets event is electron (1) or muon (2).
  LeptonType fTypeLepton;

  /// Pointer to resolution function for hadronic b quark.
  ResolutionBase * fResEnergyBhad;

  /// Pointer to resolution function for leptonic b quark.
  ResolutionBase * fResEnergyBlep;

  /// Pointer to resolution function for first light quark jet.
  ResolutionBase * fResEnergyLQ;

  /// Pointer to resolution function for the lepton.
  ResolutionBase * fResLepton;

  /// Pointer to resolution function for MET.
  ResolutionBase * fResMET;

  /** @} */
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

  double lq_meas_e;
  double lq_meas_p;
  double lq_meas_m;
  double lq_meas_deteta;
  double lq_meas_eta;
  double lq_meas_phi;
  double lq_meas_px;
  double lq_meas_py;
  double lq_meas_pz;

  double lep_meas_e;
  double lep_meas_deteta;
  double lep_meas_sintheta;
  double lep_meas_pt;
  double lep_meas_px;
  double lep_meas_py;
  double lep_meas_pz;

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

  double lq_fit_e;
  double lq_fit_px;
  double lq_fit_py;
  double lq_fit_pz;

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

  /** @} */
};
}  // namespace KLFitter

#endif  // KLFITTER_BOOSTEDLIKELIHOODTOPLEPTONJETS_H_
