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

#include "KLFitter/DetectorSnowmass.h"

#include <cmath>
#include <iostream>

#include "KLFitter/ResolutionBase.h"
#include "KLFitter/ResGauss.h"
#include "KLFitter/ResGaussE.h"
#include "KLFitter/ResGaussPt.h"
#include "KLFitter/ResGauss_MET.h"
#include "TString.h"

namespace KLFitter {
// ---------------------------------------------------------
DetectorSnowmass::DetectorSnowmass(std::string folder) : DetectorBase() {
  std::cout << "Using TFs from SnowMass ..." << std::endl;
  m_res_jet_energy_eta1      = std::unique_ptr<ResolutionBase>(new ResGaussE{Form("%s/par_energy_jets_eta1.txt", folder.c_str())});
  m_res_jet_energy_eta2      = std::unique_ptr<ResolutionBase>(new ResGaussE{Form("%s/par_energy_jets_eta2.txt", folder.c_str())});
  m_res_jet_energy_eta3      = std::unique_ptr<ResolutionBase>(new ResGaussE{Form("%s/par_energy_jets_eta3.txt", folder.c_str())});

  m_res_electron_energy_eta1 = std::unique_ptr<ResolutionBase>(new ResGaussE{Form("%s/par_energy_electrons_eta1.txt", folder.c_str())});
  m_res_electron_energy_eta2 = std::unique_ptr<ResolutionBase>(new ResGaussE{Form("%s/par_energy_electrons_eta2.txt", folder.c_str())});

  m_res_muon_momentum_eta1   = std::unique_ptr<ResolutionBase>(new ResGaussPt{Form("%s/par_pt_muons_eta1.txt", folder.c_str())});
  m_res_muon_momentum_eta2   = std::unique_ptr<ResolutionBase>(new ResGaussPt{Form("%s/par_pt_muons_eta2.txt", folder.c_str())});

  m_res_missing_ET      = std::unique_ptr<ResolutionBase>(new ResGauss_MET{Form("%s/par_misset.txt", folder.c_str())});
}

// ---------------------------------------------------------
DetectorSnowmass::~DetectorSnowmass() = default;

// ---------------------------------------------------------
ResolutionBase* DetectorSnowmass::ResEnergyLightJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    return m_res_jet_energy_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    return m_res_jet_energy_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    return m_res_jet_energy_eta3.get();
  } else {
    std::cout << "DetectorSnowmass::ResEnergyLightJet(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorSnowmass::ResEnergyBJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    return m_res_jet_energy_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    return m_res_jet_energy_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    return m_res_jet_energy_eta3.get();
  } else {
    std::cout << "DetectorSnowmass::ResEnergyBJet(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorSnowmass::ResEnergyElectron(double eta) {
  if (fabs(eta) < m_electron_eta_bin_1) {
    return m_res_electron_energy_eta1.get();
  } else if (fabs(eta) < m_electron_eta_bin_2) {
    return m_res_electron_energy_eta2.get();
  } else {
    std::cout << "DetectorSnowmass::ResEnergyElectron(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorSnowmass::ResEnergyMuon(double eta) {
  if (fabs(eta) < m_muon_eta_bin_1) {
    return m_res_muon_momentum_eta1.get();
  } else if (fabs(eta) < m_muon_eta_bin_2) {
    return m_res_muon_momentum_eta2.get();
  } else {
    std::cout << "DetectorSnowmass::ResEnergyMuon(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorSnowmass::ResMissingET() {
  return m_res_missing_ET.get();
}
}  // namespace KLFitter
