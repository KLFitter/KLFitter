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

#include "KLFitter/ResGauss.h"

#include <iostream>

#include "TMath.h"

// ---------------------------------------------------------
KLFitter::ResGauss::ResGauss(const char * filename) : KLFitter::ResolutionBase(1) {
  // read parameters from file
  ReadParameters(filename, 1);
}

// ---------------------------------------------------------
KLFitter::ResGauss::ResGauss(double sigma) : KLFitter::ResolutionBase(1) {
  // set parameter
  SetPar(0, sigma);
}

// ---------------------------------------------------------
KLFitter::ResGauss::~ResGauss() = default;

// ---------------------------------------------------------
double KLFitter::ResGauss::GetSigma(double /*par*/) {
  return fParameters[0];
}

// ---------------------------------------------------------
double KLFitter::ResGauss::logp(double x, double xmeas, bool *good, double /*par*/) {
  static constexpr double logSqrtTwoPi = 0.5*std::log(2*M_PI);

  const double sigma = GetSigma(x);
  *good = true;

  return -logSqrtTwoPi - std::log(sigma) - (0.5*(xmeas-x)*(xmeas-x)/sigma/sigma);
}
