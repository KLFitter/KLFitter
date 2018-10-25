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

#ifndef KLFITTER_RESGAUSS_H_
#define KLFITTER_RESGAUSS_H_

#include "KLFitter/ResolutionBase.h"

// ---------------------------------------------------------

namespace KLFitter {
/**
 * \class KLFitter::ResGauss
 * \brief A class describing a Gaussian resolution.
 *
 * This class offers a simple parameterization of a resolution. The
 * parameterization is a Gaussian with a width of a constant times the
 * square root of the true parameter.
 */
class ResGauss : public ResolutionBase {
 public:
  /// The default constructor.
  explicit ResGauss(const char * filename);

  /**
   * A constructor.
   * @param sigma The width of the Gaussian.
   */
  explicit ResGauss(double sigma);

  /// The (defaulted) destructor.
  ~ResGauss();

  /** \name Member functions (Get)  */
  /** @{ */

  /**
   * Return the width of the TF depending on the measured value of x.
   * Used to adjust the range of the fit parameter that correspond to the TF.
   * @param par Dummy parameter. Only needed to satisfy the interface.
   * @return The width.
   */
  double GetSigma(double par = 0) override;

  /**
   * Return the probability of the true value of x given the
   * measured value, xmeas.
   * @param x The true value of x.
   * @param xmeas The measured value of x.
   * @param good False if problem with TF.
   * @param par Optional additional parameter (not used here).
   * @return The probability.
   */
  double p(double x, double xmeas, bool *good, double par = 0) override;

  /** @} */
  /** \name Member functions (Set)  */
  /** @{ */

  /**
   * Set the width of the Gaussian
   * @param sigma The width of the Gaussian.
   */
  void SetSigma(double sigma) {
    if (sigma < 0)
      sigma = - sigma;
    this -> SetPar(0, sigma);
  }

  /** @} */
};
}  // namespace KLFitter

#endif  // KLFITTER_RESGAUSS_H_
