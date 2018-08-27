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

#ifndef KLFITTER_LIKELIHOODTOPLEPTONJETS_JETANGLES_H_
#define KLFITTER_LIKELIHOODTOPLEPTONJETS_JETANGLES_H_

#include <vector>

#include "KLFitter/LikelihoodTopLeptonJets.h"

// ---------------------------------------------------------

namespace KLFitter {
/**
  * This is an extension of LikelihoodTopLeptonJets, which adds eta and phi
  * values of the jets to the fit. The additional parameters need the
  * corresponding parameterizations in the transfer functions.
  */
class LikelihoodTopLeptonJets_JetAngles : public KLFitter::LikelihoodTopLeptonJets {
 public:
  /// Enumerator for the fitted parameters of this likelihood.
  enum Parameters { parBhadE,     ///< Energy of the hadronic b quark
                    parBlepE,     ///< Energy of the leptonic b quark
                    parLQ1E,      ///< Energy of the light quark 1
                    parLQ2E,      ///< Energy of the light quark 2
                    parLepE,      ///< Energy of the lepton
                    parNuPx,      ///< p_x of the neutrino
                    parNuPy,      ///< p_y of the neutrino
                    parNuPz,      ///< p_z of the neutrino
                    parBhadEta,   ///< Eta of the hadronic b quark
                    parBlepEta,   ///< Eta of the leptonic b quark
                    parLQ1Eta,    ///< Eta of the light quark 1
                    parLQ2Eta,    ///< Eta of the light quark 2
                    parBhadPhi,   ///< Phi of the hadronic b quark
                    parBlepPhi,   ///< Phi of the leptonic b quark
                    parLQ1Phi,    ///< Phi of the light quark 1
                    parLQ2Phi,    ///< Phi of the light quark 2
                    parTopM       ///< Mass of the top quark
  };

  /** \name Constructors and destructors */
  /** @{ */

  /// The (defaulted) constructor.
  LikelihoodTopLeptonJets_JetAngles();

  /// The (defaulted) destructor.
  ~LikelihoodTopLeptonJets_JetAngles();

  /** @} */
  /** \name Member functions (misc)  */
  /** @{ */

  /**
   * Define the parameters of the fit. This calls BCModel::AddParameter() for
   * all parameters in the enum #Parameters. Reimplemented with respect to
   * LikelihoodTopLeptonJets::DefineParameters() due to the additional
   * parameters in #Parameters.
   */
  void DefineParameters() override;

  /**
   * Identical to LikelihoodTopLeptonJets::GetInitialParametersWoNeutrinoPz(),
   * but adds the initial values for the additional angular parameters.
   * @return vector of initial values.
   */
  std::vector<double> GetInitialParametersWoNeutrinoPz() override;

  /**
   * The posterior probability definition, overloaded from BCModel. With respect
   * to LikelihoodTopLeptonJets::LogLikelihood(), this also adds resolutions
   * terms for the jet eta and phi values.
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
   *   7)  TFeta_bhad \n
   *   8)  TFeta_blep \n
   *   9)  TFeta_lq1 \n
   *  10)  TFeta_lq2 \n
   *  11)  TFphi_bhad \n
   *  12)  TFphi_blep \n
   *  13)  TFphi_lq1 \n
   *  14)  TFphi_lq2 \n
   *  15)  BW_Whad \n
   *  16)  BW_Wlep \n
   *  17)  BW_Thad \n
   *  18)  BW_Tlep \n
   */
  std::vector<double> LogLikelihoodComponents(std::vector <double> parameters) override;

  /// Request the necessary resolution functions from the detector.
  void RequestResolutionFunctions() override;

  /** @} */

 protected:
  /**
   * Adjust ranges of the parameters. This either takes the parameter sigmas
   * from the transfer functions, if #m_flag_get_par_sigmas_from_TFs is set.
   * Otherwise, a fixed sigma is used to calculate the parameter ranges. Then,
   * SetParameterRange() is called to set the range. If #m_flag_top_mass_fixed
   * is set, the top mass will be fixed to the pole mass. Reimplemented with
   * respect to LikelihoodTopLeptonJets::AdjustParameterRanges() due to
   * additional parameters in #Parameters.
   * @return An error code.
   */
  int AdjustParameterRanges() override;

  /**
   * Update 4-vector values of the model particles. This updates the internal
   * variables, such as #bhad_fit_px or #thad_fit_m to the latest values from
   * the fit. These variables are for example used in the log likelihood
   * calculation. As opposed to
   * LikelihoodTopLeptonJets::CalculateLorentzVectors(), the 4-vector values for
   * the jets include the fitted angular information as well.
   * @return An error flag.
   */
  int CalculateLorentzVectors(std::vector <double> const& parameters) override;
};
}  // namespace KLFitter

#endif  // KLFITTER_LIKELIHOODTOPLEPTONJETS_JETANGLES_H_
