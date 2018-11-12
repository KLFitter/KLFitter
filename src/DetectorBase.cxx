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

#include "KLFitter/Resolutions/ResolutionBase.h"

namespace KLFitter {
// ---------------------------------------------------------
DetectorBase::DetectorBase(std::string /*folder*/) {}

// ---------------------------------------------------------
DetectorBase::~DetectorBase() = default;

// ---------------------------------------------------------
int DetectorBase::Status() {
  try {
    if (res_type_requested.find(ResolutionType::EnergyLightJet) != res_type_requested.end()) {
      if (!ResEnergyLightJet(0)->Status()) {
        ResolutionParametersUnavailable("ResEnergyLightJet");
      }
    }
    if (res_type_requested.find(ResolutionType::EnergyBJet) != res_type_requested.end()) {
      if (!ResEnergyBJet(0)->Status()) {
        ResolutionParametersUnavailable("ResEnergyBJet");
      }
    }
    if (res_type_requested.find(ResolutionType::EnergyGluonJet) != res_type_requested.end()) {
      if (!ResEnergyGluonJet(0)->Status()) {
        ResolutionParametersUnavailable("ResEnergyGluonJet");
      }
    }
    if (res_type_requested.find(ResolutionType::EnergyElectron) != res_type_requested.end()) {
      if (!ResEnergyElectron(0)->Status()) {
        ResolutionParametersUnavailable("ResEnergyElectron");
      }
    }
    if (res_type_requested.find(ResolutionType::EnergyMuon) != res_type_requested.end()) {
      if (!ResEnergyMuon(0)->Status()) {
        ResolutionParametersUnavailable("ResEnergyMuon");
      }
    }
    if (res_type_requested.find(ResolutionType::EnergyPhoton) != res_type_requested.end()) {
      if (!ResEnergyPhoton(0)->Status()) {
        ResolutionParametersUnavailable("ResEnergyPhoton");
      }
    }
    if (res_type_requested.find(ResolutionType::MissingET) != res_type_requested.end()) {
      if (!ResMissingET()->Status()) {
        ResolutionParametersUnavailable("ResMissingET");
      }
    }
    if (res_type_requested.find(ResolutionType::EtaLightJet) != res_type_requested.end()) {
      if (!ResEtaLightJet(0)->Status()) {
        ResolutionParametersUnavailable("ResEtaLightJet");
      }
    }
    if (res_type_requested.find(ResolutionType::EtaBJet) != res_type_requested.end()) {
      if (!ResEtaBJet(0)->Status()) {
        ResolutionParametersUnavailable("ResEtaBJet");
      }
    }
    if (res_type_requested.find(ResolutionType::PhiLightJet) != res_type_requested.end()) {
      if (!ResPhiLightJet(0)->Status()) {
        ResolutionParametersUnavailable("ResPhiLightJet");
      }
    }
    if (res_type_requested.find(ResolutionType::PhiBJet) != res_type_requested.end()) {
      if (!ResPhiBJet(0)->Status()) {
        ResolutionParametersUnavailable("ResPhiBJet");
      }
    }
  } catch (std::invalid_argument& excp) {
    // Catch exception and exit gracefully with an error code.
    std::cerr << "KLFitter::DetectorBase:";
    std::cerr << " Caught std::invalid_argument: " << excp.what() << "\n";
    return 0;
  }

  // no error
  return 1;
}

// ---------------------------------------------------------
void DetectorBase::RequestResolutionType(const ResolutionType& type) {
  res_type_requested.insert(type);
}

// ---------------------------------------------------------
void DetectorBase::ResolutionParametersUnavailable(const std::string& type) {
  std::cerr << "KLFitter::DetectorBase: Parametrization of \"" << type;
  std::cerr << "\" needed by chosen likelihood, but parameters could not be read.";
  std::cerr << " Maybe the used transfer-function parameter set does not provide them?\n";
  throw std::invalid_argument(type);
}

// ---------------------------------------------------------
ResolutionBase* DetectorBase::ResolutionUndefined(const std::string& type) {
  std::cerr << "KLFitter::DetectorBase: Resolution object of type \"" << type;
  std::cerr << "\" requested by likelihood, but not supported by the detector class.\n";
  throw std::invalid_argument(type + " is undefined");
  return nullptr;
}
}  // namespace KLFitter
