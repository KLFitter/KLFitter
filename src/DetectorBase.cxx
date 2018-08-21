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

#include "KLFitter/DetectorBase.h"

#include <iostream>

#include "KLFitter/ResolutionBase.h"

namespace KLFitter {
// ---------------------------------------------------------
DetectorBase::DetectorBase(std::string /*folder*/)
  : fResEnergyLightJet(0)
  , fResEnergyBJet(0)
  , fResEnergyGluonJet(0)
  , fResEnergyElectron(0)
  , fResEnergyMuon(0)
  , fResEnergyPhoton(0)
  , fResMissingET(0) {
}

// ---------------------------------------------------------
DetectorBase::~DetectorBase() = default;

// ---------------------------------------------------------
int DetectorBase::SetResEnergyBJet(ResolutionBase * res) {
  // set resolution
  fResEnergyBJet = res;

  // no error
  return 1;
}

// ---------------------------------------------------------
int DetectorBase::SetResEnergyLightJet(ResolutionBase * res) {
  // set resolution
  fResEnergyLightJet = res;

  // no error
  return 1;
}

// ---------------------------------------------------------
int DetectorBase::SetResEnergyGluonJet(ResolutionBase * res) {
  // set resolution
  fResEnergyGluonJet = res;

  // no error
  return 1;
}

// ---------------------------------------------------------
int DetectorBase::SetResEnergyElectron(ResolutionBase * res) {
  // set resolution
  fResEnergyElectron = res;

  // no error
  return 1;
}

// ---------------------------------------------------------
int DetectorBase::SetResEnergyMuon(ResolutionBase * res) {
  // set resolution
  fResEnergyMuon = res;

  // no error
  return 1;
}

// ---------------------------------------------------------
int DetectorBase::SetResEnergyPhoton(ResolutionBase * res) {
  // set resolution
  fResEnergyPhoton = res;

  // no error
  return 1;
}

// ---------------------------------------------------------
int DetectorBase::SetResMissingET(ResolutionBase * res) {
  // set resolution
  fResMissingET = res;

  // no error
  return 1;
}

// ---------------------------------------------------------
int DetectorBase::Status() {
  if (!fResEnergyLightJet) {
    std::cout << "DetectorBase::Status(). Energy resolution of light jets not defined." << std::endl;
    return 0;
  }

  if (!fResEnergyBJet) {
    std::cout << "DetectorBase::Status(). Energy resolution of b jets not defined." << std::endl;
    return 0;
  }

  // if (!fResEnergyGluonJet) {
  //   std::cout << "DetectorBase::Status(). Energy resolution of gluon jets not defined." << std::endl;
  //   return 0;
  // }

  if (!fResEnergyElectron) {
    std::cout << "DetectorBase::Status(). Energy resolution of electrons not defined." << std::endl;
    return 0;
  }

  if (!fResEnergyMuon) {
    std::cout << "DetectorBase::Status(). Energy resolution of muons not defined." << std::endl;
    return 0;
  }

  // if (!fResEnergyPhoton) {
  //   std::cout << "DetectorBase::Status(). Energy resolution of photons not defined." << std::endl;
  //   return 0;
  // }

  if (!fResMissingET) {
    std::cout << "DetectorBase::Status(). Missing ET resolution not defined." << std::endl;
    return 0;
  }

  // no error
  return 1;
}
}  // namespace KLFitter
