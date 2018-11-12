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

#include "KLFitter/ResGaussE.h"

#include <iostream>
#include "TMath.h"

// ---------------------------------------------------------
KLFitter::ResGaussE::ResGaussE(const char * filename) : KLFitter::ResolutionBase(3) {
  // read parameters from file
  ReadParameters(filename, 3);
}
// ---------------------------------------------------------
KLFitter::ResGaussE::ResGaussE(std::vector<double> const& parameters) :KLFitter::ResolutionBase(parameters) {
  // check number of parameters
  if (parameters.size() != 3) {
    std::cout << "KLFitter::ResGaussE::ResGaussE(). Number of parameters != 3." << std::endl;
    return;
  }
}
// ---------------------------------------------------------1
KLFitter::ResGaussE::~ResGaussE() = default;

// ---------------------------------------------------------
double KLFitter::ResGaussE::GetSigma(double par) {
  return sqrt(fParameters[0]*fParameters[0]*par*par + fParameters[1]*fParameters[1]*par + fParameters[2]*fParameters[2]);
}

// ---------------------------------------------------------
double KLFitter::ResGaussE::logp(double x, double xmeas, bool *good, double /*par*/) {
  static constexpr double logSqrtTwoPi = 0.5*std::log(2*M_PI);

  *good = true;
  double sigma = GetSigma(x);
  return -logSqrtTwoPi - std::log(sigma) - (0.5*(xmeas-x)*(xmeas-x)/sigma/sigma);
}
