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

#ifndef KLFITTER_RESSINGLEGAUSSLINEARBASE_H_
#define KLFITTER_RESSINGLEGAUSSLINEARBASE_H_

#include <iostream>
#include <vector>

#include "KLFitter/ResolutionBase.h"

// ---------------------------------------------------------

/**
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter {
/**
  * \class KLFitter::ResDoubleGaussBase
  * \brief A class describing a resolution parameterized with a double Gaussian.
  *
  * This class offers a simple parameterization of a resolution. The
  * parameterization is a double Gaussian with energy dependent
  * parameters.
  */
class ResSingleGaussLinearBase : public ResolutionBase {
 public:
  /** \name Constructors and destructors */
  /* @{ */

  /**
    * The default constructor.
    */
  explicit ResSingleGaussLinearBase(const char * filename);

  /**
    * A constructor.
    * @param parameters The parameters of the parameterization.
    */
  explicit ResSingleGaussLinearBase(std::vector<double> const& parameters);

  /**
    * The (defaulted) destructor.
    */
  virtual ~ResSingleGaussLinearBase();

  /* @} */
  /** \name Member functions (Get)  */
  /* @{ */

  /**
    * Calculate the mean of the Gaussian from the TF parameters and the value of x.
    * @param x The value of x.
    * @return The mean.
    */
  virtual double GetMean(double x) = 0;

  /**
    * Calculate the width of the Gaussian from the TF parameters and the value of x.
    * @param x The value of x.
    * @return The width.
    */
  virtual double GetSigma(double x) override = 0;

  /**
    * Return the probability of the true value of x given the
    * measured value, xmeas.
    * @param x The true value of x.
    * @param xmeas The measured value of x.
    * @param good False if problem with TF.
    * @param par Optional additional parameter (not used here).
    * @return Logarithm of the probability.
    */
  double logp(double x, double xmeas, bool *good, double /*par*/ = 0) override;

  /* @} */

  /**
    * Sanity check for double gaussian parameters p2, p3 and p5 (1st sigma, scale and 2nd sigma).
    * @param sigma1 (the 1st sigma).
    * @param amplitude2 (the scale parameter).
    * @param sigma2 (the 2nd sigma).
    * @return False if problem with TF.
    */
  static bool CheckSingleGaussianSanity(double *sigma) {
    if (*sigma <= 0.) {
      *sigma = 0.00000001;
      return false;
    }

    return true;
  }
};
}  // namespace KLFitter

#endif  // KLFITTER_RESSINGLEGAUSSLINEARBASE_H_
