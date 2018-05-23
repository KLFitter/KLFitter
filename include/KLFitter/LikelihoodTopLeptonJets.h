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

#ifndef KLFITTER_LIKELIHOODTOPLEPTONJETS_H_
#define KLFITTER_LIKELIHOODTOPLEPTONJETS_H_

#include <iostream>
#include <vector>

#include "KLFitter/LikelihoodBase.h"

class TLorentzVector;

// ---------------------------------------------------------

namespace KLFitter {
class ResolutionBase;

/**
 * A class implementing a likelihood for the ttbar lepton+jets channel. This
 * class represents a likelihood for the ttbar into lepton+jets.
 */
class LikelihoodTopLeptonJets : public LikelihoodBase {
 public:
  /// \name Constructors and destructors
  /// @{

  /**
   * The default constructor. This initializes all member attributes and calls
   * the functions DefineModelParticles() and DefineParameters().
   */
  LikelihoodTopLeptonJets();

  /// The (defaulted) destructor.
  ~LikelihoodTopLeptonJets();

  /// @}

  /// Enumerator for the lepton type.
  enum LeptonType { kElectron,  ///< Leptons of type electron
                    kMuon       ///< Leptons of type muon
  };

  /// Enumerator for the fitted parameters of this likelihood.
  enum Parameters { parBhadE,   ///< Energy of the hadronic b quark
                    parBlepE,   ///< Energy of the leptonic b quark
                    parLQ1E,    ///< Energy of the light quark 1
                    parLQ2E,    ///< Energy of the light quark 2
                    parLepE,    ///< Energy of the lepton
                    parNuPx,    ///< p_x of the neutrino
                    parNuPy,    ///< p_y of the neutrino
                    parNuPz,    ///< p_z of the neutrino
                    parTopM     ///< Mass of the top quark
  };

  /// \name Member functions (Set)
  /// @{

  /**
   * Set a flag. If flag is true the invariant top quark mass is
   * fixed to the pole mass.
   * @param flag The flag.
   */
  void SetFlagTopMassFixed(bool flag) { m_flag_top_mass_fixed = flag; }

  /**
   * Set a flag. If flag is true, take the parameter sigma values from the
   * transfer functions.
   * @param flag The flag.
   */
  void SetFlagGetParSigmasFromTFs(bool flag) { m_flag_get_par_sigmas_from_TFs = flag; }

  /**
   * Set the type of lepton
   * @param leptontype The type of lepton: kElectron or kMuon
   */
  void SetLeptonType(LeptonType leptontype);

  /**
   * Set the type of lepton
   * @param leptontype The type of lepton: electron(1) or muon (2)
   */
  void SetLeptonType(int leptontype);

  /// @}
  /// \name Member functions (BAT)
  /// @{

  /**
   * Define the parameters of the fit. This calls BCModel::AddParameter() for
   * all parameters in the enum #Parameters.
   */
  void DefineParameters() override;

  /**
   * Get initial values for the parameters. This calls
   * GetInitialParametersWoNeutrinoPz() and retrieves the neutrino solutions
   * via GetNeutrinoPzSolutions(). The initial values are then set accordingly.
   * @return vector of initial values.
   */
  std::vector<double> GetInitialParameters() override;

  /**
   * Get initial values for the parameters with a dummy of "0.0" for the neutrino pz.
   * The decision on the initial value for the neutrino pz then needs to be done in
   * GetInitialParameters().
   * @return vector of initial values.
   */
  virtual std::vector<double> GetInitialParametersWoNeutrinoPz();

  /**
   * The posterior probability definition, overloaded from BCModel.
   * @param parameters A vector of parameters (double values).
   * @return The logarithm of the prior probability.
   */
  double LogLikelihood(const std::vector <double> & parameters) override;

  /**
   * The posterior probability definition, overloaded from BCModel. Instead of
   * the final log likelihood value as in LogLikelihood(), this returns all
   * subcomponents.
   * @param parameters A vector of parameters (double values).
   * @return A vector with the components of the logarithm of the prior
   * probability. Its components are: \n
   *   0) TF_bhad \n
   *   1) TF_blep \n
   *   2) TF_lq1 \n
   *   3) TF_lq2 \n
   *   4) TF_lep \n
   *   5) TF_METx \n
   *   6) TF_METy \n
   *   7) BW_Whad \n
   *   8) BW_Wlep \n
   *   9) BW_Thad \n
   *   10) BW_Tlep
   */
  std::vector<double> LogLikelihoodComponents(std::vector <double> parameters) override;

  /// @}
  /// \name Member functions (misc)
  /// @{

  /**
   * Set the values for the missing ET x and y components and the m_et_miss_sum. This
   * sets the internal variables #m_et_miss_x, #m_et_miss_y and #m_et_miss_sum to the
   * given values.
   * @param etx missing ET x component.
   * @param ety missing ET y component.
   * @param sumet total scalar ET.
   * @return An error flag.
   */
  int SetET_miss_XY_SumET(double etx, double ety, double sumet) override;

  /// @}

 protected:
  /**
   * Adjust ranges of the parameters. This either takes the parameter sigmas
   * from the transfer functions, if #m_flag_get_par_sigmas_from_TFs is set.
   * Otherwise, a fixed sigma is used to calculate the parameter ranges. Then,
   * SetParameterRange() is called to set the range. If #m_flag_top_mass_fixed is
   * set, the top mass will be fixed to the pole mass.
   * @return An error code.
   */
  int AdjustParameterRanges() override;

  /**
   * Build the model particles from the best fit parameters. This sets the
   * particles of #fParticlesModel to the fitted parameter values. The W boson
   * and top quark parameters are combined from the other particles.
   * @return An error code.
   */
  int BuildModelParticles() override;

  /**
   * Update 4-vector values of the model particles. This updates the internal
   * variables, such as #m_bhad_fit_px or #m_thad_fit_m to the latest values from
   * the fit. These variables are for example used in the log likelihood
   * calculation.
   * @return An error flag.
   */
  int CalculateLorentzVectors(std::vector <double> const& parameters) override;

  /**
   * Calculates the neutrino pz solutions from the measured values and the W
   * mass. An additional particle to be added to the charged lepton may be
   * specified, for example a photon in ttbargamma, which is radiated from the
   * leptonic W or the charged lepton;
   * @param additionalParticle Pointer to a 4-vector of a particle which is
   * added to the charged lepton in the calculation.
   * @return A vector with 0, 1 or 2 neutrino pz solutions.
   */
  std::vector<double> CalculateNeutrinoPzSolutions(TLorentzVector* additionalParticle = nullptr);

  /**
   * Define the model particles. Create the object #fParticlesModel and add all
   * particles of this likelihood. The 4-vector components are set to zero.
   * @return An error code.
   */
  int DefineModelParticles() override;

  /**
   * Return the neutrino pz solutions from the measured values and the W mass.
   * This calls CalculateNeutrinoPzSolutions() and returns the results.
   * @return A vector with 0, 1 or 2 neutrino pz solutions.
   */
  std::vector<double> GetNeutrinoPzSolutions();

  /**
   * Remove the invariant particle permutations. In particular, this removes
   * permutations where the two light jets are swapped. In addition, in case
   * two types of leptons were added, only permutations with leptons of type
   * #m_lepton_type are kept.
   * @return An error code.
   */
  int RemoveInvariantParticlePermutations() override;

  /**
   * Save permuted particles. This takes the permuted particles, stored in
   * #fParticlesPermuted, and saves their values in the internal variables,
   * such as #m_bhad_meas_eta.
   * @return An error code.
   */
  int SavePermutedParticles() override;

  /**
   * Save the resolution functions from the detector to the internal pointers.
   * This sets the internal pointers, such as #m_res_energy_bhad to the resolution
   * objects stored in the detector.
   * @return An error code.
   */
  int SaveResolutionFunctions() override;

  /// \name Member attributes
  /// @{

  /// A flag for using a fixed top mass (true) or not (false).
  bool m_flag_top_mass_fixed;

  /// Flag for using ResolutionBase::GetSigma() to retrieve the parameter ranges
  bool m_flag_get_par_sigmas_from_TFs;

  /// The values of the x component of the missing ET.
  double m_et_miss_x;

  /// The values of the y component of the missing ET.
  double m_et_miss_y;

  /// The values of the total scalar ET.
  double m_et_miss_sum;

  /// An index deciding if the event is electron (1) or muon (2) plus jets.
  LeptonType m_lepton_type;

  /// Pointer to resolution function for hadronic b quark.
  ResolutionBase* m_res_energy_bhad;

  /// Pointer to resolution function for leptonic b quark.
  ResolutionBase* m_res_energy_blep;

  /// Pointer to resolution function for first light quark jet.
  ResolutionBase* m_res_energy_lq1;

  /// Pointer to resolution function for second light quark jet.
  ResolutionBase* m_res_energy_lq2;

  /// Pointer to resolution function for the lepton.
  ResolutionBase* m_res_lepton;

  /// Pointer to resolution function for MET.
  ResolutionBase* m_res_met;

  /// @}
  /// \name Member attributes (measured parameters)
  /// @{

  double m_bhad_meas_e;
  double m_bhad_meas_p;
  double m_bhad_meas_m;
  double m_bhad_meas_deteta;
  double m_bhad_meas_eta;
  double m_bhad_meas_phi;
  double m_bhad_meas_px;
  double m_bhad_meas_py;
  double m_bhad_meas_pz;

  double m_blep_meas_e;
  double m_blep_meas_p;
  double m_blep_meas_m;
  double m_blep_meas_deteta;
  double m_blep_meas_eta;
  double m_blep_meas_phi;
  double m_blep_meas_px;
  double m_blep_meas_py;
  double m_blep_meas_pz;

  double m_lq1_meas_e;
  double m_lq1_meas_p;
  double m_lq1_meas_m;
  double m_lq1_meas_deteta;
  double m_lq1_meas_eta;
  double m_lq1_meas_phi;
  double m_lq1_meas_px;
  double m_lq1_meas_py;
  double m_lq1_meas_pz;

  double m_lq2_meas_e;
  double m_lq2_meas_p;
  double m_lq2_meas_m;
  double m_lq2_meas_deteta;
  double m_lq2_meas_eta;
  double m_lq2_meas_phi;
  double m_lq2_meas_px;
  double m_lq2_meas_py;
  double m_lq2_meas_pz;

  double m_lep_meas_e;
  double m_lep_meas_deteta;
  double m_lep_meas_sintheta;
  double m_lep_meas_pt;
  double m_lep_meas_px;
  double m_lep_meas_py;
  double m_lep_meas_pz;

  /// @}
  /// \name Member attributes (fitted parameters)
  /// @{

  double m_bhad_fit_e;
  double m_bhad_fit_px;
  double m_bhad_fit_py;
  double m_bhad_fit_pz;

  double m_blep_fit_e;
  double m_blep_fit_px;
  double m_blep_fit_py;
  double m_blep_fit_pz;

  double m_lq1_fit_e;
  double m_lq1_fit_px;
  double m_lq1_fit_py;
  double m_lq1_fit_pz;

  double m_lq2_fit_e;
  double m_lq2_fit_px;
  double m_lq2_fit_py;
  double m_lq2_fit_pz;

  double m_lep_fit_e;
  double m_lep_fit_px;
  double m_lep_fit_py;
  double m_lep_fit_pz;

  double m_nu_fit_e;
  double m_nu_fit_px;
  double m_nu_fit_py;
  double m_nu_fit_pz;

  double m_wlep_fit_e;
  double m_wlep_fit_px;
  double m_wlep_fit_py;
  double m_wlep_fit_pz;

  double m_whad_fit_e;
  double m_whad_fit_px;
  double m_whad_fit_py;
  double m_whad_fit_pz;

  double m_whad_fit_m;
  double m_wlep_fit_m;
  double m_thad_fit_m;
  double m_tlep_fit_m;

  /// @}
};
}  // namespace KLFitter

#endif  // KLFITTER_LIKELIHOODTOPLEPTONJETS_H_
