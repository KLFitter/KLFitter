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

#include "KLFitter/Resolutions/Gauss/GaussPt.h"

#include <iostream>
#include "TMath.h"

// ---------------------------------------------------------
KLFitter::Resolutions::Gauss::GaussPt::GaussPt(const char * filename) : KLFitter::Resolutions::Gauss::GaussBase(filename) { }

// ---------------------------------------------------------
KLFitter::Resolutions::Gauss::GaussPt::GaussPt(std::vector<double> const& parameters) : KLFitter::Resolutions::Gauss::GaussBase(parameters) { }

// ---------------------------------------------------------1
KLFitter::Resolutions::Gauss::GaussPt::~GaussPt() = default;

// ---------------------------------------------------------
double KLFitter::Resolutions::Gauss::GaussPt::GetSigma(double par) {
  if (par <= 200.0)
    return fParameters[0];
  else
    return fParameters[1];
}
