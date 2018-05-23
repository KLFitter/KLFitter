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

/**
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter {
class ResolutionBase;

/**
  * \class KLFitter::LikelihoodTopLeptonJets_Angular
  * \brief Add brief description here
  *
  * Add detailed description here.
  */
class LikelihoodTopLeptonJets_Angular : public KLFitter::LikelihoodTopLeptonJets {
 public:
  /** \name Constructors and destructors */
  /* @{ */

  /**
    * The default constructor.
    */
  LikelihoodTopLeptonJets_Angular();

  /**
    * The (defaulted) destructor.
    */
  ~LikelihoodTopLeptonJets_Angular();

  /* @} */
  /** \name Member functions (misc)  */
  /* @{ */

  /**
    * The posterior probability definition, overloaded from BCModel.
    * @param parameters A vector of parameters (double values).
    * @return The logarithm of the prior probability.
    */
  double LogLikelihood(const std::vector <double> & parameters) override;

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

 protected:
  double wlep_fit_e;
  double wlep_fit_px;
  double wlep_fit_py;
  double wlep_fit_pz;

  double whad_fit_e;
  double whad_fit_px;
  double whad_fit_py;
  double whad_fit_pz;
};
}  // namespace KLFitter

#endif  // KLFITTER_LIKELIHOODTOPLEPTONJETS_ANGULAR_H_
