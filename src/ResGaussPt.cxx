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

#include "KLFitter/ResGaussPt.h"

#include <iostream>
#include "TMath.h"

// ---------------------------------------------------------
KLFitter::ResGaussPt::ResGaussPt(const char * filename) : KLFitter::ResolutionBase(2) {
  // read parameters from file
  ReadParameters(filename, 2);
}
// ---------------------------------------------------------
KLFitter::ResGaussPt::ResGaussPt(std::vector<double> const& parameters) :KLFitter::ResolutionBase(parameters) {
  // check number of parameters
  if (parameters.size() != 2) {
    std::cout << "KLFitter::ResGaussPt::ResGaussPt(). Number of parameters != 2." << std::endl;
    return;
  }
}
// ---------------------------------------------------------1
KLFitter::ResGaussPt::~ResGaussPt() {
}

// ---------------------------------------------------------
double KLFitter::ResGaussPt::GetSigma(double x) {
  if (x <= 200.0)
    return fParameters[0];
  else
    return fParameters[1];
}

// ---------------------------------------------------------
double KLFitter::ResGaussPt::p(double x, double xmeas, bool *good) {
  *good = true;
   double sigma = GetSigma(x);
  return TMath::Gaus(xmeas, x, sigma, true);
}
