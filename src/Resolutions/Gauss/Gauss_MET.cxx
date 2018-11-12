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

#include "KLFitter/Resolutions/Gauss/Gauss_MET.h"

#include <iostream>

#include "TMath.h"

// ---------------------------------------------------------
KLFitter::Resolutions::Gauss::Gauss_MET::Gauss_MET(const char * filename) : KLFitter::Resolutions::Gauss::GaussBase(filename) {}

// ---------------------------------------------------------
KLFitter::Resolutions::Gauss::Gauss_MET::Gauss_MET(std::vector<double> const& parameters) : KLFitter::Resolutions::Gauss::GaussBase(parameters) { }

// ---------------------------------------------------------
KLFitter::Resolutions::Gauss::Gauss_MET::~Gauss_MET() = default;

// ---------------------------------------------------------
double KLFitter::Resolutions::Gauss::Gauss_MET::GetSigma(double sumet) {
  return fParameters[0]+fParameters[1]/(1+exp(-fParameters[2]*(sumet-fParameters[3])));
}
