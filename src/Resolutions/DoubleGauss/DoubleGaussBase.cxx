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

#include "KLFitter/Resolutions/DoubleGauss/DoubleGaussBase.h"

#include <cmath>
#include <iostream>

// ---------------------------------------------------------
KLFitter::Resolutions::DoubleGauss::DoubleGaussBase::DoubleGaussBase(const char * filename) : KLFitter::Resolutions::ResolutionBase(10) {
  // read parameters from file
  ReadParameters(filename, 10);
}

// ---------------------------------------------------------
KLFitter::Resolutions::DoubleGauss::DoubleGaussBase::DoubleGaussBase(std::vector<double> const& parameters) : KLFitter::Resolutions::ResolutionBase(parameters) {
  // check number of parameters
  if (parameters.size() != 10) {
    std::cout << "KLFitter::Resolutions::DoubleGauss::DoubleGaussBase::DoubleGaussBase(). Number of parameters != 10." << std::endl;
    return;
  }
}

// ---------------------------------------------------------
KLFitter::Resolutions::DoubleGauss::DoubleGaussBase::~DoubleGaussBase() = default;

// ---------------------------------------------------------
double KLFitter::Resolutions::DoubleGauss::DoubleGaussBase::GetSigma(double par) {
  // Calculate mean width of both gaussians; weight the width of the 2nd one with its amplitude
  double sigma1 = GetSigma1(par);
  double sigma2 = GetSigma2(par);
  double amplitude2 = GetAmplitude2(par);
  double sigma = (sigma1 + amplitude2*sigma2) / (1+amplitude2);

  // sigma estimates the fractional resolution, but we want absolute
  return sigma*par;
}

// ---------------------------------------------------------
double KLFitter::Resolutions::DoubleGauss::DoubleGaussBase::logp(double x, double xmeas, bool *good, double /*par*/) {
  double m1 = GetMean1(x);
  double s1 = GetSigma1(x);
  double a2 = GetAmplitude2(x);
  double m2 = GetMean2(x);
  double s2 = GetSigma2(x);

  // sanity checks for p2, p3 and p5
  *good = CheckDoubleGaussianSanity(&s1, &a2, &s2);

  double dx = (x - xmeas) / x;

  // calculate double-Gaussian
  const double p =  1./sqrt(2.*M_PI) / (s1 + a2 * s2) * (exp(-(dx-m1)*(dx-m1)/(2 * s1*s1)) + a2 * exp(-(dx-m2)*(dx-m2)/(2 * s2 * s2)));
  return std::log(p);
}
