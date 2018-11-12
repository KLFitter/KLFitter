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

#include "KLFitter/Resolutions/ResGauss_linearMET.h"

#include <cmath>
#include <iostream>

// ---------------------------------------------------------
KLFitter::ResGauss::ResGauss_linearMET::ResGauss_linearMET(const char * filename) : KLFitter::ResGauss::ResGaussBase(filename) { }

// ---------------------------------------------------------
KLFitter::ResGauss::ResGauss_linearMET::ResGauss_linearMET(std::vector<double> const& parameters) : KLFitter::ResGauss::ResGaussBase(parameters) { }

// ---------------------------------------------------------
KLFitter::ResGauss::ResGauss_linearMET::~ResGauss_linearMET() = default;

// ---------------------------------------------------------
double KLFitter::ResGauss::ResGauss_linearMET::GetSigma(double x) {
  return fParameters[0] + fParameters[1] * x;
}
