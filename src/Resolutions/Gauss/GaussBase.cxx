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

#include "KLFitter/Resolutions/Gauss/GaussBase.h"

#include <iostream>

#include "TMath.h"

// ---------------------------------------------------------
KLFitter::Resolutions::Gauss::GaussBase::GaussBase(const char * filename) : KLFitter::Resolutions::ResolutionBase(1) {
  // read parameters from file
  ReadParameters(filename, 1);
}

// ---------------------------------------------------------
KLFitter::Resolutions::Gauss::GaussBase::GaussBase(std::vector<double> const& parameters) : KLFitter::Resolutions::ResolutionBase(parameters) {
  // check number of parameters
  if (parameters.size() != 2) {
    std::cout << "KLFitter::Resolutions::Gauss::GaussBase::GaussBase(). Number of parameters != 2." << std::endl;
    return;
  }
}

// ---------------------------------------------------------
KLFitter::Resolutions::Gauss::GaussBase::GaussBase(double sigma) : KLFitter::Resolutions::ResolutionBase(1) {
  // set parameter
  SetPar(0, sigma);
}

// ---------------------------------------------------------
KLFitter::Resolutions::Gauss::GaussBase::~GaussBase() = default;

// ---------------------------------------------------------
double KLFitter::Resolutions::Gauss::GaussBase::GetSigma(double /*par*/) {
  return fParameters[0];
}

// ---------------------------------------------------------
double KLFitter::Resolutions::Gauss::GaussBase::logp(double x, double xmeas, bool *good, double /*par*/) {
  static constexpr double logSqrtTwoPi = 0.5*std::log(2*M_PI);

  const double sigma = GetSigma(x);
  *good = true;

  return -logSqrtTwoPi - std::log(sigma) - (0.5*(xmeas-x)*(xmeas-x)/sigma/sigma);
}
