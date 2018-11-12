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

#include "KLFitter/Resolutions/ResCrystalBallJets.h"

#include <cmath>
#include <iostream>

// ---------------------------------------------------------
KLFitter::ResCrystalBall::ResCrystalBallJets::ResCrystalBallJets(const char * filename) : KLFitter::ResCrystalBall::ResCrystalBallBase(filename) { }

// ---------------------------------------------------------
KLFitter::ResCrystalBall::ResCrystalBallJets::ResCrystalBallJets(std::vector<double> const& parameters) : KLFitter::ResCrystalBall::ResCrystalBallBase(parameters) { }

// ---------------------------------------------------------
KLFitter::ResCrystalBall::ResCrystalBallJets::~ResCrystalBallJets() = default;

// ---------------------------------------------------------
double KLFitter::ResCrystalBall::ResCrystalBallJets::GetAlpha(double x) {
  return fParameters[0] + fParameters[1] * x;
}

// ---------------------------------------------------------
double KLFitter::ResCrystalBall::ResCrystalBallJets::GetN(double x) {
  return fParameters[2] + fParameters[3] * x;
}

// ---------------------------------------------------------
double KLFitter::ResCrystalBall::ResCrystalBallJets::GetSigma(double x) {
  return fParameters[4] + fParameters[5] / std::sqrt(x);
}

// ---------------------------------------------------------
double KLFitter::ResCrystalBall::ResCrystalBallJets::GetMean(double x) {
  return fParameters[6] + fParameters[7] / std::sqrt(x);
}
