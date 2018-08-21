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

#include "KLFitter/DetectorAtlas_7TeV.h"

#include <stdlib.h>

#include <cmath>
#include <iostream>

#include "KLFitter/ResDoubleGaussE_1.h"
#include "KLFitter/ResDoubleGaussE_2.h"
#include "KLFitter/ResDoubleGaussE_3.h"
#include "KLFitter/ResDoubleGaussPt.h"
#include "KLFitter/ResGauss.h"
#include "KLFitter/ResGauss_MET.h"
#include "KLFitter/ResolutionBase.h"
#include "TString.h"

namespace KLFitter {
// ---------------------------------------------------------
DetectorAtlas_7TeV::DetectorAtlas_7TeV(std::string folder) : DetectorBase() {
  // check: powheg sample with 7TeV? Must use 8!
  if (strstr(folder.c_str(), "mc11c_powheg")) {
    std::cout << "ERROR! Don't use PowHeg TFs with the 7TeV Detector class!!! Exiting..." << std::endl;
    exit(1);
  }
  // check: MC11b? New parametrization!
  if ((strstr(folder.c_str(), "mc11b")) || (strstr(folder.c_str(), "mc11c"))) {
    std::cout << "Using TF from MC11b or later..." << std::endl;
    m_res_energy_bjet_eta1     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_bJets_eta1.txt", folder.c_str())});
    m_res_energy_bjet_eta2     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_bJets_eta2.txt", folder.c_str())});
    m_res_energy_bjet_eta3     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_bJets_eta3.txt", folder.c_str())});
    m_res_energy_bjet_eta4     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_bJets_eta4.txt", folder.c_str())});
    m_res_energy_bjet_eta5     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_bJets_eta5.txt", folder.c_str())});
  } else {
    std::cout << "Using TF from MC11a or earlier..." << std::endl;
    m_res_energy_bjet_eta1     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_2{Form("%s/par_energy_bJets_eta1.txt", folder.c_str())});
    m_res_energy_bjet_eta2     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_2{Form("%s/par_energy_bJets_eta2.txt", folder.c_str())});
    m_res_energy_bjet_eta3     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_2{Form("%s/par_energy_bJets_eta3.txt", folder.c_str())});
    m_res_energy_bjet_eta4     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_2{Form("%s/par_energy_bJets_eta4.txt", folder.c_str())});
    m_res_energy_bjet_eta5     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_2{Form("%s/par_energy_bJets_eta5.txt", folder.c_str())});
  }

  // Remaining energy resolutions
  m_res_energy_light_jet_eta1 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_lJets_eta1.txt", folder.c_str())});
  m_res_energy_light_jet_eta2 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_lJets_eta2.txt", folder.c_str())});
  m_res_energy_light_jet_eta3 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_lJets_eta3.txt", folder.c_str())});
  m_res_energy_light_jet_eta4 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_lJets_eta4.txt", folder.c_str())});
  m_res_energy_light_jet_eta5 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_lJets_eta5.txt", folder.c_str())});

  m_res_energy_gluon_jet_eta1 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_gluon_eta1.txt", folder.c_str())});
  m_res_energy_gluon_jet_eta2 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_gluon_eta2.txt", folder.c_str())});
  m_res_energy_gluon_jet_eta3 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_gluon_eta3.txt", folder.c_str())});
  m_res_energy_gluon_jet_eta4 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_gluon_eta4.txt", folder.c_str())});

  m_res_energy_electron_eta1 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_Electrons_eta1.txt", folder.c_str())});
  m_res_energy_electron_eta2 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_Electrons_eta2.txt", folder.c_str())});
  m_res_energy_electron_eta3 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_Electrons_eta3.txt", folder.c_str())});
  m_res_energy_electron_eta4 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_Electrons_eta4.txt", folder.c_str())});

  m_res_energy_muon_eta1     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussPt{Form("%s/par_energy_Muons_eta1.txt", folder.c_str())});
  m_res_energy_muon_eta2     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussPt{Form("%s/par_energy_Muons_eta2.txt", folder.c_str())});
  m_res_energy_muon_eta3     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussPt{Form("%s/par_energy_Muons_eta3.txt", folder.c_str())});

  m_res_energy_photon_eta1   = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_energy_photon_eta1.txt", folder.c_str())});
  m_res_energy_photon_eta2   = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_energy_photon_eta2.txt", folder.c_str())});
  m_res_energy_photon_eta3   = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_energy_photon_eta3.txt", folder.c_str())});
  m_res_energy_photon_eta4   = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_energy_photon_eta4.txt", folder.c_str())});

  // eta resolution
  m_res_eta_light_jet_eta1 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_lJets_eta1.txt", folder.c_str())});
  m_res_eta_light_jet_eta2 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_lJets_eta2.txt", folder.c_str())});
  m_res_eta_light_jet_eta3 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_lJets_eta3.txt", folder.c_str())});
  m_res_eta_light_jet_eta4 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_lJets_eta4.txt", folder.c_str())});

  m_res_eta_bjet_eta1 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_bJets_eta1.txt", folder.c_str())});
  m_res_eta_bjet_eta2 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_bJets_eta2.txt", folder.c_str())});
  m_res_eta_bjet_eta3 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_bJets_eta3.txt", folder.c_str())});
  m_res_eta_bjet_eta4 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_bJets_eta4.txt", folder.c_str())});

  // phi resolution
  m_res_phi_light_jet_eta1 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_lJets_eta1.txt", folder.c_str())});
  m_res_phi_light_jet_eta2 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_lJets_eta2.txt", folder.c_str())});
  m_res_phi_light_jet_eta3 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_lJets_eta3.txt", folder.c_str())});
  m_res_phi_light_jet_eta4 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_lJets_eta4.txt", folder.c_str())});

  m_res_phi_bjet_eta1 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_bJets_eta1.txt", folder.c_str())});
  m_res_phi_bjet_eta2 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_bJets_eta2.txt", folder.c_str())});
  m_res_phi_bjet_eta3 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_bJets_eta3.txt", folder.c_str())});
  m_res_phi_bjet_eta4 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_bJets_eta4.txt", folder.c_str())});

  m_res_missing_ET = std::unique_ptr<ResolutionBase>(new ResGauss_MET{Form("%s/par_misset.txt", folder.c_str())});
}

// ---------------------------------------------------------
DetectorAtlas_7TeV::~DetectorAtlas_7TeV() = default;

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_7TeV::ResEnergyLightJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    return m_res_energy_light_jet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    return m_res_energy_light_jet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    return m_res_energy_light_jet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    return m_res_energy_light_jet_eta4.get();
  } else if (fabs(eta) <= m_jet_eta_bin_5) {
    return m_res_energy_light_jet_eta5.get();
  } else {
    std::cout << "DetectorAtlas_7TeV::ResEnergyLightJet(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_7TeV::ResEnergyBJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    return m_res_energy_bjet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    return m_res_energy_bjet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    return m_res_energy_bjet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    return m_res_energy_bjet_eta4.get();
  } else if (fabs(eta) <= m_jet_eta_bin_5) {
    return m_res_energy_bjet_eta5.get();
  } else {
    std::cout << "DetectorAtlas_7TeV::ResEnergyBJet(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_7TeV::ResEnergyGluonJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    return m_res_energy_gluon_jet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    return m_res_energy_gluon_jet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    return m_res_energy_gluon_jet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    return m_res_energy_gluon_jet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_7TeV::ResEnergyGluonJet(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_7TeV::ResEnergyElectron(double eta) {
  if (fabs(eta) < m_electron_eta_bin_1) {
    return m_res_energy_electron_eta1.get();
  } else if (fabs(eta) < m_electron_eta_bin_2) {
    return m_res_energy_electron_eta2.get();
  } else if (fabs(eta) < m_electron_eta_bin_3) {
    std::cout << "DetectorAtlas_7TeV::ResEnergyElectron(). Electron in crack region" << std::endl;
    return nullptr;
  } else if (fabs(eta) <= m_electron_eta_bin_4) {
    return m_res_energy_electron_eta4.get();
  } else {
    std::cout << "DetectorAtlas_7TeV::ResEnergyElectron(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_7TeV::ResEnergyMuon(double eta) {
  if (fabs(eta) < m_muon_eta_bin_1) {
    return m_res_energy_muon_eta1.get();
  } else if (fabs(eta) < m_muon_eta_bin_2) {
    return m_res_energy_muon_eta2.get();
  } else if (fabs(eta) < m_muon_eta_bin_3) {
    return m_res_energy_muon_eta3.get();
  } else {
    std::cout << "DetectorAtlas_7TeV::ResEnergyMuon(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_7TeV::ResEnergyPhoton(double eta) {
  if (fabs(eta) < m_photon_eta_bin_1) {
    return m_res_energy_photon_eta1.get();
  } else if (fabs(eta) < m_photon_eta_bin_2) {
    return m_res_energy_photon_eta2.get();
  } else if (fabs(eta) < m_photon_eta_bin_3) {
    return m_res_energy_photon_eta3.get();
  } else if (fabs(eta) <= m_photon_eta_bin_4) {
    return m_res_energy_photon_eta4.get();
  } else {
    std::cout << "DetectorAtlas_7TeV::ResEnergyPhoton(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_7TeV::ResEtaLightJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    return m_res_eta_light_jet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    return m_res_eta_light_jet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    return m_res_eta_light_jet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    return m_res_eta_light_jet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_7TeV::ResEtaLightJet(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_7TeV::ResEtaBJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    return m_res_eta_bjet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    return m_res_eta_bjet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    return m_res_eta_bjet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    return m_res_eta_bjet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_7TeV::ResEtaBJet(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_7TeV::ResPhiLightJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    return m_res_phi_light_jet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    return m_res_phi_light_jet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    return m_res_phi_light_jet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    return m_res_phi_light_jet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_7TeV::ResPhiLightJet(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_7TeV::ResPhiBJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    return m_res_phi_bjet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    return m_res_phi_bjet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    return m_res_phi_bjet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    return m_res_phi_light_jet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_7TeV::ResPhiBJet(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_7TeV::ResMissingET() {
  return m_res_missing_ET.get();
}
}  // namespace KLFitter
