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

#include "KLFitter/Resolutions/ResDoubleGaussE_1.h"

#include <cmath>
#include <iostream>

// ---------------------------------------------------------
KLFitter::ResDoubleGauss::ResDoubleGaussE_1::ResDoubleGaussE_1(const char * filename) : KLFitter::ResDoubleGauss::ResDoubleGaussBase(filename) { }

// ---------------------------------------------------------
KLFitter::ResDoubleGauss::ResDoubleGaussE_1::ResDoubleGaussE_1(std::vector<double> const& parameters) : KLFitter::ResDoubleGauss::ResDoubleGaussBase(parameters) { }

// ---------------------------------------------------------
KLFitter::ResDoubleGauss::ResDoubleGaussE_1::~ResDoubleGaussE_1() = default;

// ---------------------------------------------------------
double KLFitter::ResDoubleGauss::ResDoubleGaussE_1::GetMean1(double x) {
  return fParameters[0] + fParameters[1] * x;
}

// ---------------------------------------------------------
double KLFitter::ResDoubleGauss::ResDoubleGaussE_1::GetSigma1(double x) {
  return fParameters[2] / sqrt(x) + fParameters[3];
}

// ---------------------------------------------------------
double KLFitter::ResDoubleGauss::ResDoubleGaussE_1::GetAmplitude2(double x) {
  return fParameters[4] + fParameters[5] * x;
}

// ---------------------------------------------------------
double KLFitter::ResDoubleGauss::ResDoubleGaussE_1::GetMean2(double x) {
  return fParameters[6] + fParameters[7] * x;
}

// ---------------------------------------------------------
double KLFitter::ResDoubleGauss::ResDoubleGaussE_1::GetSigma2(double x) {
  return fParameters[8] + fParameters[9] * x;
}
