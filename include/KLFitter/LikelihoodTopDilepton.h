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

#ifndef KLFITTER_LIKELIHOODTOPDILEPTON_H_
#define KLFITTER_LIKELIHOODTOPDILEPTON_H_

#include <assert.h>

#include <iostream>
#include <memory>
#include <utility>
#include <vector>

class TLorentzVector;

#include "BAT/BCH1D.h"
#include "BAT/BCModel.h"
#include "KLFitter/LikelihoodBase.h"

// ---------------------------------------------------------

namespace KLFitter {
class NuSolutions;  // defined in implementation file
class ResolutionBase;

/**
 * Likelihood for the ttbar dilepton channel, where both tops
 * decay leptonically. ADD MORE DETAILED DESCRIPTION HERE
 * (particularly about the techniques used to determine the
 * neutrino solutions)
 */
class LikelihoodTopDilepton : public KLFitter::LikelihoodBase {
 public:
  /// The default constructor.
  LikelihoodTopDilepton();

  /// The (defaulted) destructor.
  ~LikelihoodTopDilepton();

  /** \name Member functions (Set)  */
  /** @{ */

  /// Enumerator for the lepton type.
  enum LeptonType { kElectron, kMuon };

  /// Enumerator for the parameters.
  enum Parameters { parTopM = 0, parB1E, parB2E, parLep1E, parLep2E, parAntiNuEta, parNuEta };

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

  /**
   * Set the neutrino pseudorapidity sigma linear dependency on mtop
   * according to SM expectations
   */
  void SetEtaNuParams(std::vector<double> etanuparam) { nueta_params = etanuparam; }

  /**
   * Set the type of lepton
   * @param leptontype_1 The type of the first lepton: kElectron or kMuon
   * @param leptontype_2 The type of the second lepton: kElectron or kMuon
   */
  void SetLeptonType(LeptonType leptontype_1, LeptonType leptontype_2);

  /**
   * Set the type of lepton
   * @param leptontype_1 The type of the first lepton: electron(1) or muon (2)
   * @param leptontype_2 The type of the second lepton: electron(1) or muon (2)
   */
  void SetLeptonType(int leptontype_1, int leptontype_2);

  /** @} */
  /** \name Member functions (BAT)  */
  /** @{ */

  /// Define the parameters of the fit.
  void DefineParameters() override;

  /// Define sharp gauss prior for mtop par if mtop fixed
  void DefinePrior();

  /// Define BCH1D and TH1D histograms to be filled in MCMCIterationInterface
  void DefineHistograms();

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
   * 0:  NuWT
   * 1:  TF_b1
   * 2:  TF_b2
   * 3:  TF_lep1
   * 4:  TF_lep2
   * 5:  AntiNu_Eta
   * 6:  Nu_Eta
   * 7:  Minv(lep,jet)
   */
  std::vector<double> LogLikelihoodComponents(std::vector <double> parameters) override;

  /// Get initial values for the parameters.
  std::vector<double> GetInitialParameters() override;

  /**
   * Return Gaussian term for neutrino
   * pseudorapidity.
   * @return A double.
   */
  double GaussNuEta(std::vector<double> parameters);

  /// Return Gaussian term for antineutrino pseudorapidity.
  double GaussAntiNuEta(std::vector<double> parameters);

  /// Return NuWT weight
  double CalculateWeight(const std::vector<double> & parameters);

  /// Return NuWT weight for a set of jet1, jet2, lep1, lep2
  double CalculateWeightPerm(TLorentzVector * l1, TLorentzVector * l2, TLorentzVector * j1, TLorentzVector * j2, const std::vector<double> & parameters);

  /// Return set of neutrino/antineutrino kinematic solutions (up to 2).
  KLFitter::NuSolutions SolveForNuMom(TLorentzVector * l, TLorentzVector * b, double mtop, double nueta);

  /// Return neutrino weight for a given nu solution and antinu solution
  double neutrino_weight(TLorentzVector nu, TLorentzVector nubar);

  /// Return sum of invariant masses of each (lep,jet) pair, including a tuning factor alpha.
  double CalculateMLepJet();

  /**
   * Set a flag. If flag is true the sumloglikelihood
   * option is used, instead of the default best-permutation
   * @param flag The flag.
   */
  void SetDoSumLogLik(bool flag) { doSumloglik = flag; }

  /** @} */

 protected:
  /** \name Member functions (misc)  */
  /** @{ */

  /**
   * Update 4-vectors of model particles.
   * @return An error flag.
   */
  int CalculateLorentzVectors(std::vector <double> const& parameters) override;

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

  /// Calculate other variables out of the KLFitter parameters for each MCMCiteration
  void MCMCIterationInterface() override;

  /// Get BAT BCH1D histograms of Mttbar
  BCH1D * GetHistMttbar() { return fHistMttbar.get(); }

  /// Get BAT BCH1D histograms of CosTheta
  BCH1D * GetHistCosTheta() { return fHistCosTheta.get(); }

  /// Calculate cos(theta*) for both top and antitop
  std::pair<float, float> CalculateCosTheta(std::vector <TLorentzVector> *particles);

  /** @} */

  /// Save permuted particles.
  int SavePermutedParticles() override;

  /// Save resolution functions.
  int SaveResolutionFunctions() override;

  /** \name Member attributes */
  /** @{ */

  /// A flag for using a fixed top mass (true) or not (false).
  bool fFlagTopMassFixed;

  /// The values of the x component of the missing ET.
  double ETmiss_x;

  /// The values of the y component of the missing ET.
  double ETmiss_y;

  /// The values of the total scalar ET.
  double SumET;

  /// Lepton 1 Type (electron or muon)
  LeptonType fTypeLepton_1;

  /// Lepton 2 Type (electron or muon)
  LeptonType fTypeLepton_2;

  /**
   * vector including nu pseudorapidity sigma
   * dependency on mtop ( if sigma = a + b*mtop = >
   * vector[0] = a, vector[1] = b)
   */
  std::vector<double> nueta_params;

  /// A flag for using sumloglikelihood option
  bool doSumloglik;

  /**
   * TH1D histograms to be filled
   * with functions of interest, e.g.: mttbar,
   * costheta*,etc. for each MCNCiteration
   */
  TH1D *  hist_mttbar;
  TH1D *  hist_costheta;

  /// BAT BCH1D Histogram for mttbar
  std::unique_ptr<BCH1D> fHistMttbar;

  /// BAT BCH1D Histogram cos(theta*)
  std::unique_ptr<BCH1D> fHistCosTheta;

  /// BAT BCH1D Histogram for dR(truth top, fit top)
  BCH1D * fHistdRTop;

  /// BAT BCH1D Histogram for dR(truth antitop, fit antitop)
  BCH1D * fHistdRAntiTop;

  /// BAT BCH1D Histogram for dR(truth nu, fit nu)
  BCH1D * fHistdRNu;

  /// BAT BCH1D Histogram for dR(truth antinu, fit antinu)
  BCH1D * fHistdRAntiNu;

  /// Pointer to resolution function for the first b quark.
  ResolutionBase * fResEnergyB1;

  /// Pointer to resolution function for the second b quark.
  ResolutionBase * fResEnergyB2;

  /// Pointer to resolution function for the first lepton.
  ResolutionBase * fResLepton1;

  /// Pointer to resolution function for the second lepton.
  ResolutionBase * fResLepton2;

  /// Pointer to resolution function for MET.
  ResolutionBase * fResMET;

  /** @} */
  /** \name Member attributes (measured parameters) */
  /** @{ */

  double b1_meas_e;
  double b1_meas_p;
  double b1_meas_m;
  double b1_meas_deteta;
  double b1_meas_eta;
  double b1_meas_phi;
  double b1_meas_px;
  double b1_meas_py;
  double b1_meas_pz;

  double b2_meas_e;
  double b2_meas_p;
  double b2_meas_m;
  double b2_meas_deteta;
  double b2_meas_eta;
  double b2_meas_phi;
  double b2_meas_px;
  double b2_meas_py;
  double b2_meas_pz;

  double lep1_meas_e;
  double lep1_meas_deteta;
  float  lep1_meas_charge;
  double lep1_meas_sintheta;
  double lep1_meas_pt;
  double lep1_meas_px;
  double lep1_meas_py;
  double lep1_meas_pz;

  double lep2_meas_e;
  double lep2_meas_deteta;
  float  lep2_meas_charge;
  double lep2_meas_sintheta;
  double lep2_meas_pt;
  double lep2_meas_px;
  double lep2_meas_py;
  double lep2_meas_pz;

  /** @} */
  /** \name Member attributes (fitted parameters) */
  /** @{ */

  double lep1_fit_e;
  double lep1_fit_px;
  double lep1_fit_py;
  double lep1_fit_pz;

  double lep2_fit_e;
  double lep2_fit_px;
  double lep2_fit_py;
  double lep2_fit_pz;

  double b1_fit_e;
  double b1_fit_px;
  double b1_fit_py;
  double b1_fit_pz;

  double b2_fit_e;
  double b2_fit_px;
  double b2_fit_py;
  double b2_fit_pz;

  /** @} */
};
}  // namespace KLFitter

#endif  // KLFITTER_LIKELIHOODTOPDILEPTON_H_
