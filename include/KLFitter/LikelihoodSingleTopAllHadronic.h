/*
 * Copyright (c) 2009--2021, the KLFitter developer team
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

#ifndef KLFITTER_LIKELIHOODSINGLETOPALLHADRONIC_H_
#define KLFITTER_LIKELIHOODSINGLETOPALLHADRONIC_H_

#include <vector>

#include "KLFitter/LikelihoodBase.h"

// ---------------------------------------------------------

namespace KLFitter {
class ResolutionBase;

/**
 * A class implementing a likelihood for the single-top production in the
 * hadronic decay channel.
 */
class LikelihoodSingleTopAllHadronic : public LikelihoodBase {
 public:
  /** \name Constructors and destructors */
  /* @{ */

  /**
   * The default constructor. This initializes all member attributes and calls
   * the functions DefineModelParticles() and DefineParameters().
   */
  LikelihoodSingleTopAllHadronic();

  /// The (defaulted) destructor.
  ~LikelihoodSingleTopAllHadronic();

  /* @} */

  /// Enumerator for the fitted parameters of this likelihood.
  enum Parameters { parBhadE,   ///< Energy of the hadronic b quark
                    parLQ1E,    ///< Energy of the light quark 1
                    parLQ2E,    ///< Energy of the light quark 2
                    parTopM     ///< Mass of the top quark
  };

  /** \name Member functions (Set)  */
  /* @{ */

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

  /* @} */
  /** \name Member functions (BAT)  */
  /* @{ */

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
   *   1) TF_lq1 \n
   *   2) TF_lq2 \n
   *   3) BW_Whad \n
   *   4) BW_Thad
   */
  std::vector<double> LogLikelihoodComponents(std::vector <double> parameters) override;

  /* @} */
  /** \name Member functions (misc)  */
  /* @{ */

  /**
   * Set the values for the missing ET x and y components and the m_et_miss_sum. This
   * sets the internal variables #m_et_miss_x, #m_et_miss_y and #m_et_miss_sum to the
   * given values.
   * @param etx missing ET x component.
   * @param ety missing ET y component.
   * @param sumet total scalar ET.
   * @return An error flag.
   */
  int SetET_miss_XY_SumET(double /*etx*/, double /*ety*/, double /*sumet*/) override { return 1; }

  /// Request the necessary resolution functions from the detector.
  void RequestResolutionFunctions() override;

  /* @} */

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
   * Define the model particles. Create the object #fParticlesModel and add all
   * particles of this likelihood. The 4-vector components are set to zero.
   * @return An error code.
   */
  int DefineModelParticles() override;

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

  /** \name Member attributes */
  /* @{ */

  /// A flag for using a fixed top mass (true) or not (false).
  bool m_flag_top_mass_fixed;

  /// Flag for using ResolutionBase::GetSigma() to retrieve the parameter ranges
  bool m_flag_get_par_sigmas_from_TFs;

  /// Pointer to resolution function for hadronic b quark.
  ResolutionBase* m_res_energy_bhad;

  /// Pointer to resolution function for first light quark jet.
  ResolutionBase* m_res_energy_lq1;

  /// Pointer to resolution function for second light quark jet.
  ResolutionBase* m_res_energy_lq2;

  /* @} */
  /** \name Member attributes (measured parameters) */
  /* @{ */

  double m_bhad_meas_e;
  double m_bhad_meas_p;
  double m_bhad_meas_m;
  double m_bhad_meas_deteta;
  double m_bhad_meas_eta;
  double m_bhad_meas_phi;
  double m_bhad_meas_px;
  double m_bhad_meas_py;
  double m_bhad_meas_pz;

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

  /* @} */
  /** \name Member attributes (fitted parameters) */
  /* @{ */

  double m_bhad_fit_e;
  double m_bhad_fit_px;
  double m_bhad_fit_py;
  double m_bhad_fit_pz;

  double m_lq1_fit_e;
  double m_lq1_fit_px;
  double m_lq1_fit_py;
  double m_lq1_fit_pz;

  double m_lq2_fit_e;
  double m_lq2_fit_px;
  double m_lq2_fit_py;
  double m_lq2_fit_pz;

  double m_whad_fit_e;
  double m_whad_fit_px;
  double m_whad_fit_py;
  double m_whad_fit_pz;

  double m_whad_fit_m;
  double m_thad_fit_m;

  /* @} */
};
}  // namespace KLFitter

#endif  // KLFITTER_LIKELIHOODSINGLETOPALLHADRONIC_H_
