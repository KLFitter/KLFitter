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

#include "KLFitter/Resolutions/DoubleGauss/DoubleGaussE_3.h"

#include <cmath>
#include <iostream>

// ---------------------------------------------------------
KLFitter::Resolutions::DoubleGauss::DoubleGaussE_3::DoubleGaussE_3(const char * filename) : KLFitter::Resolutions::DoubleGauss::DoubleGaussBase(filename) { }

// ---------------------------------------------------------
KLFitter::Resolutions::DoubleGauss::DoubleGaussE_3::DoubleGaussE_3(std::vector<double> const& parameters) : KLFitter::Resolutions::DoubleGauss::DoubleGaussBase(parameters) { }

// ---------------------------------------------------------
KLFitter::Resolutions::DoubleGauss::DoubleGaussE_3::~DoubleGaussE_3() = default;

// ---------------------------------------------------------
double KLFitter::Resolutions::DoubleGauss::DoubleGaussE_3::GetMean1(double x) {
  return fParameters[0] / x + fParameters[1];
}

// ---------------------------------------------------------
double KLFitter::Resolutions::DoubleGauss::DoubleGaussE_3::GetSigma1(double x) {
  return sqrt(fParameters[2]*fParameters[2] / x + fParameters[3]*fParameters[3]);
}

// ---------------------------------------------------------
double KLFitter::Resolutions::DoubleGauss::DoubleGaussE_3::GetAmplitude2(double x) {
  return fParameters[4] + fParameters[5] * x;
}

// ---------------------------------------------------------
double KLFitter::Resolutions::DoubleGauss::DoubleGaussE_3::GetMean2(double x) {
  return fParameters[6] / x + fParameters[7];
}

// ---------------------------------------------------------
double KLFitter::Resolutions::DoubleGauss::DoubleGaussE_3::GetSigma2(double x) {
  return fParameters[8]/ x + fParameters[9];
}
