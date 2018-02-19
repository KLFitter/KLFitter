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

#include "KLFitter/ResDoubleGaussE_4.h"
#include <iostream>

#include <cmath>
// used for lJets/bjets in mc12
// ---------------------------------------------------------
KLFitter::ResDoubleGaussE_4::ResDoubleGaussE_4(const char * filename) : KLFitter::ResDoubleGaussBase(filename) { }

// ---------------------------------------------------------
KLFitter::ResDoubleGaussE_4::ResDoubleGaussE_4(std::vector<double> const& parameters) : KLFitter::ResDoubleGaussBase(parameters) { }

// ---------------------------------------------------------
KLFitter::ResDoubleGaussE_4::~ResDoubleGaussE_4() {
  ;
}

// ---------------------------------------------------------
double KLFitter::ResDoubleGaussE_4::GetMean1(double x) {
  return fParameters[0] + fParameters[1] / x;
}

// ---------------------------------------------------------
double KLFitter::ResDoubleGaussE_4::GetSigma1(double x) {
  return fParameters[2] + fParameters[3] / sqrt(x);
}

// ---------------------------------------------------------
double KLFitter::ResDoubleGaussE_4::GetAmplitude2(double x) {
  return fParameters[4] + fParameters[5] / x;
}

// ---------------------------------------------------------
double KLFitter::ResDoubleGaussE_4::GetMean2(double x) {
  return fParameters[6] + fParameters[7]  / sqrt(x);
}

// ---------------------------------------------------------
double KLFitter::ResDoubleGaussE_4::GetSigma2(double x) {
  return fParameters[8] + fParameters[9] * x;
}

// ---------------------------------------------------------
