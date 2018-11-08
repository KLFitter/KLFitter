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

#ifndef KLFITTER_RESCRYSTALBALLBASE_H_
#define KLFITTER_RESCRYSTALBALLBASE_H_

#include <iostream>
#include <vector>

#include "KLFitter/ResolutionBase.h"

// ---------------------------------------------------------

namespace KLFitter {
/**
  * This class offers a simple parameterisation of a resolution. The
  * parametrisation is a CrystalBall function with energy dependent
  * parameters.
  */
class ResCrystalBallBase : public ResolutionBase {
 public:
  ///The default constructor.
  explicit ResCrystalBallBase(const char * filename);

  /**
    * A constructor that takes parameters directly, unlike the default 
    * constructor that takes a path to the file with TFs.
    * @param parameters The parameters of the parameterization.
    */
  explicit ResCrystalBallBase(std::vector<double> const& parameters);

  /**
    * The (defaulted) destructor.
    */
  virtual ~ResCrystalBallBase();

  /**
    * Calculate the alpha parameter from the TF parameters and the value of x.
    * @param x The value of x.
    * @return The alpha.
    */
  virtual double GetAlpha(double x) = 0;

  /**
    * Calculate the n parameter from the TF parameters and the value of x.
    * @param x The value of x.
    * @return The parameter n.
    */
  virtual double GetN(double x) = 0;

  /**
    * Calculate the mean from the TF parameters and the value of x.
    * @param x The value of x.
    * @return The mean.
    */
  virtual double GetMean(double x) = 0;

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

  /**
    * Sanity check for the crystal ball parameters sigma, and n (1st sigma, scale and 2nd sigma).
    * @param sigma.
    * @param n.
    * @return False if problem with TF.
    */
  static bool CheckCrystalBallSanity(double *sigma, double *n) {
    if (*sigma <= 0.) {
      *sigma = 0.000000001;
      return false;
    }
    if (*n <= 1.) {
      *n = 1.00000001;
      return false;
    }

    return true;
  }

  /**
   * A function that calculates CrystalBall
   * See twiki: https://en.wikipedia.org/wiki/Crystal_Ball_function
   * @param x
   * @param Alpha parameter
   * @param n parameter
   * @param sigma parameter
   * @param mean parameter
   * @return CrystalBall value for X
   */
  double CrystalBallFunction(double x, double alpha, double n, double sigma, double mean);

  /**
   * An approximation of the error function needed to calculate crystal ball normalization
   * with precision < 1e-4. Decreases computation time by about 10%.
   * @param x
   * @return Approximate value of the error function for x
   */
  double ApproxError(double x);
};
}  // namespace KLFitter

#endif  // KLFITTER_RESCRYSTALBALLBASE_H_
