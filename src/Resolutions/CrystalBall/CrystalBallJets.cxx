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

#include "KLFitter/Resolutions/CrystalBall/CrystalBallJets.h"

#include <cmath>
#include <iostream>

// ---------------------------------------------------------
KLFitter::Resolutions::CrystalBall::CrystalBallJets::CrystalBallJets(const char * filename) : KLFitter::Resolutions::CrystalBall::CrystalBallBase(filename) { }

// ---------------------------------------------------------
KLFitter::Resolutions::CrystalBall::CrystalBallJets::CrystalBallJets(std::vector<double> const& parameters) : KLFitter::Resolutions::CrystalBall::CrystalBallBase(parameters) { }

// ---------------------------------------------------------
KLFitter::Resolutions::CrystalBall::CrystalBallJets::~CrystalBallJets() = default;

// ---------------------------------------------------------
double KLFitter::Resolutions::CrystalBall::CrystalBallJets::GetAlpha(double x) {
  return fParameters[0] + fParameters[1] * x;
}

// ---------------------------------------------------------
double KLFitter::Resolutions::CrystalBall::CrystalBallJets::GetN(double x) {
  return fParameters[2] + fParameters[3] * x;
}

// ---------------------------------------------------------
double KLFitter::Resolutions::CrystalBall::CrystalBallJets::GetSigma(double x) {
  return fParameters[4] + fParameters[5] / std::sqrt(x);
}

// ---------------------------------------------------------
double KLFitter::Resolutions::CrystalBall::CrystalBallJets::GetMean(double x) {
  return fParameters[6] + fParameters[7] / std::sqrt(x);
}
