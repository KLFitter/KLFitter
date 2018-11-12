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

#include "KLFitter/ResSingleGaussLinearBase.h"

#include "TMath.h"

#include <cmath>
#include <iostream>

// ---------------------------------------------------------
KLFitter::ResSingleGaussLinearBase::ResSingleGaussLinearBase(const char * filename) : KLFitter::ResolutionBase(8) {
  // read parameters from file
  ReadParameters(filename, 2);
}

// ---------------------------------------------------------
KLFitter::ResSingleGaussLinearBase::ResSingleGaussLinearBase(std::vector<double> const& parameters) : KLFitter::ResolutionBase(parameters) {
  // check number of parameters
  if (parameters.size() != 2) {
    std::cout << "KLFitter::ResSingleGaussBase::ResSingleGaussBase(). Number of parameters != 2." << std::endl;
    return;
  }
}

// ---------------------------------------------------------
KLFitter::ResSingleGaussLinearBase::~ResSingleGaussLinearBase() = default;

// ---------------------------------------------------------
double KLFitter::ResSingleGaussLinearBase::logp(double x, double xmeas, bool *good, double /*par*/) {
  static constexpr double logSqrtTwoPi = 0.5*std::log(2*M_PI);

  double sigma = GetSigma(x);
  // sanity checks for p2, p3 and p5
  *good = CheckSingleGaussianSanity(&sigma);

  return -logSqrtTwoPi - std::log(sigma) - (0.5*(xmeas-x)*(xmeas-x)/sigma/sigma);
}
