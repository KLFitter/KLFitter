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

#include <iostream>
#include <vector>

#include "KLFitter/LikelihoodTopLeptonJets.h"

class TLorentzVector;

// ---------------------------------------------------------

/**
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter {
class ResolutionBase;

/**
  * \class KLFitter::LikelihoodTopLeptonJets_JetAngles
  * \brief Add brief description here
  *
  * Add a detailed description here.
  */
class LikelihoodTopLeptonJets_JetAngles : public KLFitter::LikelihoodTopLeptonJets {
 public:
  /** \name Constructors and destructors */
  /* @{ */

  /**
    * The default constructor.
    */
  LikelihoodTopLeptonJets_JetAngles();

  /**
    * The (defaulted) destructor.
    */
  ~LikelihoodTopLeptonJets_JetAngles();

  /* @} */
  /** \name Member functions (Set)  */
  /* @{ */

  /**
    * Enumerator for the parameters.
    */
  enum Parameters { parBhadE, parBlepE, parLQ1E, parLQ2E, parLepE, parNuPx, parNuPy, parNuPz, parBhadEta, parBlepEta, parLQ1Eta, parLQ2Eta, parBhadPhi, parBlepPhi, parLQ1Phi, parLQ2Phi, parTopM };

  /* @} */
  /** \name Member functions (misc)  */
  /* @{ */

  /**
    * Define the parameters of the fit.
    */
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
    * 4:  TF_lep
    * 5:  TF_METx
    * 6:  TF_METy
    * 7:  TFeta_bhad
    * 8:  TFeta_blep
    * 9:  TFeta_lq1
    *10:  TFeta_lq2
    *11:  TFphi_bhad
    *12:  TFphi_blep
    *13:  TFphi_lq1
    *14:  TFphi_lq2
    *15:  BW_Whad
    *16:  BW_Wlep
    *17:  BW_Thad
    *18:  BW_Tlep
    */
  std::vector<double> LogLikelihoodComponents(std::vector <double> parameters) override;

  /**
    * Get initial values for the parameters with a dummy of "0.0" for the neutrino pz.
    * The decision on the initial value for the neutrino pz then needs to be done in
    * GetInitialParameters().
    * @return vector of initial values.
    */
  std::vector<double> GetInitialParametersWoNeutrinoPz() override;

  /* @} */

 protected:
  /** \name Member functions (misc)  */
  /* @{ */

  /**
    * Update 4-vectors of model particles.
    * @return An error flag.
    */
  int CalculateLorentzVectors(std::vector <double> const& parameters) override;

  /**
    * Adjust parameter ranges
    */
  int AdjustParameterRanges() override;

  /* @} */
};
}  // namespace KLFitter

#endif  // KLFITTER_LIKELIHOODTOPLEPTONJETS_JETANGLES_H_
