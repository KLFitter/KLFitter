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

#include "KLFitter/ResGauss_MET.h"

#include <iostream>

#include "TMath.h"

// ---------------------------------------------------------
KLFitter::ResGauss_MET::ResGauss_MET(const char * filename) : KLFitter::ResolutionBase(4) {
  // read parameters from file
  ReadParameters(filename, 4);
}

// ---------------------------------------------------------
KLFitter::ResGauss_MET::ResGauss_MET(std::vector<double> const& parameters) :KLFitter::ResolutionBase(parameters) {
  // check number of parameters
  if (parameters.size() != 4) {
    std::cout << "KLFitter::ResGauss_MET::ResGauss_MET(). Number of parameters != 4." << std::endl;
    return;
  }
}

// ---------------------------------------------------------
KLFitter::ResGauss_MET::~ResGauss_MET() = default;

// ---------------------------------------------------------
double KLFitter::ResGauss_MET::GetSigma(double sumet) {
  return fParameters[0]+fParameters[1]/(1+exp(-fParameters[2]*(sumet-fParameters[3])));
}

// ---------------------------------------------------------
double KLFitter::ResGauss_MET::logp(double x, double xmeas, bool *good, double sumet) {
  *good = true;
  // calculate MET TF with 4 parameters (MC10b or later)
  double sigma = GetSigma(sumet);
  return TMath::Gaus(xmeas, x, sigma, true);
}
