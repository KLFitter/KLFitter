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

#ifndef KLFITTER_LIKELIHOODTOPLEPTONJETS_ANGULAR_H_
#define KLFITTER_LIKELIHOODTOPLEPTONJETS_ANGULAR_H_

#include <iostream>
#include <vector>

#include "KLFitter/LikelihoodTopLeptonJets.h"

class TLorentzVector;

// ---------------------------------------------------------

namespace KLFitter {
class ResolutionBase;

/**
  * This is an extension of LikelihoodTopLeptonJets, which adds angular
  * information to the likelihood.
  */
class LikelihoodTopLeptonJets_Angular : public KLFitter::LikelihoodTopLeptonJets {
 public:
  /// \name Constructors and destructors
  /// @{

  /// The (defaulted) constructor.
  LikelihoodTopLeptonJets_Angular();

  /// The (defaulted) destructor.
  ~LikelihoodTopLeptonJets_Angular();

  /// @}
  /// \name Member functions (misc)
  /// @{

  /**
   * The posterior probability definition, overloaded from BCModel. In addition
   * to what is implemented in LikelihoodTopLeptonJets::LogLikelihood(), this
   * also includes angular information in the calculation of the log likelihood.
   * @param parameters A vector of parameters (double values).
   * @return The logarithm of the prior probability.
   */
  double LogLikelihood(const std::vector <double> & parameters) override;

  /// @}

 protected:
  /**
   * Adjust ranges of the parameters. As opposed to
   * LikelihoodTopLeptonJets::AdjustParameterRanges(), this does \emph not
   * utilize the flag #m_flag_get_par_sigmas_from_TFs, but only uses a fixed
   * sigma to calculate the parameter range and set it with SetParameterRange().
   * If #m_flag_top_mass_fixed is set, the top mass will be fixed to the pole
   * mass.
   * @return An error code.
   */
  int AdjustParameterRanges() override;
};
}  // namespace KLFitter

#endif  // KLFITTER_LIKELIHOODTOPLEPTONJETS_ANGULAR_H_
