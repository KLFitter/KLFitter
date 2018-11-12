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

#include "KLFitter/Resolutions/ResCrystalBallBase.h"

// Needed for CrystalBall
#include "Math/Math.h"

#include <cmath>
#include <iostream>

// ---------------------------------------------------------
KLFitter::ResCrystalBall::ResCrystalBallBase::ResCrystalBallBase(const char * filename) :
  KLFitter::ResolutionBase(8) {
  // read parameters from file
  ReadParameters(filename, 8);
}

// ---------------------------------------------------------
KLFitter::ResCrystalBall::ResCrystalBallBase::ResCrystalBallBase(std::vector<double> const& parameters) :
  KLFitter::ResolutionBase(parameters) {
  // check number of parameters
  if (parameters.size() != 8) {
    std::cout << "KLFitter::ResCrystalBallBase::ResCrystalBallBase(). Number of parameters != 8." << std::endl;
    return;
  }
}

// ---------------------------------------------------------
KLFitter::ResCrystalBall::ResCrystalBallBase::~ResCrystalBallBase() = default;

// ---------------------------------------------------------
double KLFitter::ResCrystalBall::ResCrystalBallBase::logp(double x, double xmeas, bool *good, double /*par*/) {
  static constexpr double overSqrt2 = 1./std::sqrt(2.);
  static constexpr double sqrtPiHalf = std::sqrt(M_PI/2.);

  double alpha = GetAlpha(x);
  double n = GetN(x);
  double sigma = GetSigma(x);
  double mean = GetMean(x);

  // sanity checks for n and sigma
  *good = CheckCrystalBallSanity(&sigma, &n);

  double dx = (x - xmeas) / x;

  // Needed for normalization
  const double C = n/std::fabs(alpha) * 1./(n-1.) * std::exp(-0.5*alpha*alpha);
  const double D = sqrtPiHalf*(1.+ApproxError(std::fabs(alpha)*overSqrt2));

  return (-std::log(sigma*(C+D))) + LogCrystalBallFunction(dx, alpha, n, sigma, mean);
}
  
// ---------------------------------------------------------
double KLFitter::ResCrystalBall::ResCrystalBallBase::LogCrystalBallFunction(double x, double alpha,
    double n, double sigma, double mean) {
	// evaluate the crystal ball function
  if (sigma < 0.) {
    return -9999999;
  }
  double z = (x - mean)/sigma; 
  if (alpha < 0) {
    z = -z;
  }
  double abs_alpha = std::abs(alpha);
  if (z  > - abs_alpha) {
    return (- 0.5 * z * z);
  } else {
    double nDivAlpha = n/abs_alpha;
    double B = nDivAlpha -abs_alpha;
    double arg = nDivAlpha/(B-z);

    return (-0.5*abs_alpha*abs_alpha) + n*std::log(arg);
  }
}

// ---------------------------------------------------------
double KLFitter::ResCrystalBall::ResCrystalBallBase::ApproxError(double x) {
  static constexpr double a1 = 0.278393;
  static constexpr double a2 = 0.230389;
  static constexpr double a3 = 0.000972;
  static constexpr double a4 = 0.078108;

  const double denom = 1.+ a1*x + a2*x*x + a3*x*x*x + a4*x*x*x*x;

  return 1. - (1./denom/denom/denom/denom);
}
