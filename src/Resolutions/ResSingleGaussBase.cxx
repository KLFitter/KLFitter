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

#include "KLFitter/Resolutions/ResSingleGaussBase.h"

#include "TMath.h"

#include <cmath>
#include <iostream>

// ---------------------------------------------------------
KLFitter::ResSingleGauss::ResSingleGaussBase::ResSingleGaussBase(const char * filename) : KLFitter::ResolutionBase(8) {
  // read parameters from file
  ReadParameters(filename, 4);
}

// ---------------------------------------------------------
KLFitter::ResSingleGauss::ResSingleGaussBase::ResSingleGaussBase(std::vector<double> const& parameters) : KLFitter::ResolutionBase(parameters) {
  // check number of parameters
  if (parameters.size() != 4) {
    std::cout << "KLFitter::ResSingleGauss::ResSingleGaussBase::ResSingleGaussBase(). Number of parameters != 4." << std::endl;
    return;
  }
}

// ---------------------------------------------------------
KLFitter::ResSingleGauss::ResSingleGaussBase::~ResSingleGaussBase() = default;

// ---------------------------------------------------------
double KLFitter::ResSingleGauss::ResSingleGaussBase::logp(double x, double xmeas, bool *good, double /*par*/) {
  static constexpr double logSquareTwoPi = 0.5*std::log(2*M_PI);

  double mean = GetMean(x);
  double sigma = GetSigma(x);

  // sanity checks for p2, p3 and p5
  *good = CheckSingleGaussianSanity(&sigma);

  double dx = (x - xmeas) / x;

  return -logSquareTwoPi - std::log(sigma) - (0.5*(dx-mean)*(dx-mean)/sigma/sigma); 
}
