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

#include "KLFitter/DetectorAtlas_8TeV.h"

#include <cmath>
#include <iostream>

#include "KLFitter/ResolutionBase.h"
#include "KLFitter/ResDoubleGaussE_1.h"
#include "KLFitter/ResDoubleGaussE_4.h"
#include "KLFitter/ResDoubleGaussE_5.h"
#include "KLFitter/ResDoubleGaussPt.h"
#include "KLFitter/ResGauss.h"
#include "KLFitter/ResGauss_MET.h"
#include "TString.h"

namespace KLFitter {
// ---------------------------------------------------------
DetectorAtlas_8TeV::DetectorAtlas_8TeV(std::string folder) : DetectorBase() {
  std::cout << "Using TF from MC12 ..." << std::endl;
  // energy resolution
  m_res_energy_light_jet_eta1 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_4{Form("%s/par_energy_lJets_eta1.txt", folder.c_str())});
  m_res_energy_light_jet_eta2 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_4{Form("%s/par_energy_lJets_eta2.txt", folder.c_str())});
  m_res_energy_light_jet_eta3 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_4{Form("%s/par_energy_lJets_eta3.txt", folder.c_str())});
  m_res_energy_light_jet_eta4 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_4{Form("%s/par_energy_lJets_eta4.txt", folder.c_str())});

  m_res_energy_bjet_eta1     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_4{Form("%s/par_energy_bJets_eta1.txt", folder.c_str())});
  m_res_energy_bjet_eta2     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_4{Form("%s/par_energy_bJets_eta2.txt", folder.c_str())});
  m_res_energy_bjet_eta3     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_4{Form("%s/par_energy_bJets_eta3.txt", folder.c_str())});
  m_res_energy_bjet_eta4     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_4{Form("%s/par_energy_bJets_eta4.txt", folder.c_str())});

  // m_res_energy_gluon_jet_eta1 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_gluon_eta1.txt", folder.c_str())});
  // m_res_energy_gluon_jet_eta2 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_gluon_eta2.txt", folder.c_str())});
  // m_res_energy_gluon_jet_eta3 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_gluon_eta3.txt", folder.c_str())});
  // m_res_energy_gluon_jet_eta4 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_1{Form("%s/par_energy_gluon_eta4.txt", folder.c_str())});

  m_res_energy_electron_eta1 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_5{Form("%s/par_energy_Electrons_eta1.txt", folder.c_str())});
  m_res_energy_electron_eta2 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_5{Form("%s/par_energy_Electrons_eta2.txt", folder.c_str())});
  m_res_energy_electron_eta3 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_5{Form("%s/par_energy_Electrons_eta3.txt", folder.c_str())});
  m_res_energy_electron_eta4 = std::unique_ptr<ResolutionBase>(new ResDoubleGaussE_5{Form("%s/par_energy_Electrons_eta4.txt", folder.c_str())});

  m_res_energy_muon_eta1     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussPt{Form("%s/par_energy_Muons_eta1.txt", folder.c_str())});
  m_res_energy_muon_eta2     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussPt{Form("%s/par_energy_Muons_eta2.txt", folder.c_str())});
  m_res_energy_muon_eta3     = std::unique_ptr<ResolutionBase>(new ResDoubleGaussPt{Form("%s/par_energy_Muons_eta3.txt", folder.c_str())});

  // m_res_energy_photon_eta1   = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_energy_photon_eta1.txt", folder.c_str())});
  // m_res_energy_photon_eta2   = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_energy_photon_eta2.txt", folder.c_str())});
  // m_res_energy_photon_eta3   = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_energy_photon_eta3.txt", folder.c_str())});
  // m_res_energy_photon_eta4   = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_energy_photon_eta4.txt", folder.c_str())});

  // eta resolution
  // m_res_eta_light_jet_eta1 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_lJets_eta1.txt", folder.c_str())});
  // m_res_eta_light_jet_eta2 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_lJets_eta2.txt", folder.c_str())});
  // m_res_eta_light_jet_eta3 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_lJets_eta3.txt", folder.c_str())});
  // m_res_eta_light_jet_eta4 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_lJets_eta4.txt", folder.c_str())});

  // m_res_eta_bjet_eta1 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_bJets_eta1.txt", folder.c_str())});
  // m_res_eta_bjet_eta2 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_bJets_eta2.txt", folder.c_str())});
  // m_res_eta_bjet_eta3 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_bJets_eta3.txt", folder.c_str())});
  // m_res_eta_bjet_eta4 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_bJets_eta4.txt", folder.c_str())});

  // phi resolution
  // m_res_phi_light_jet_eta1 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_lJets_eta1.txt", folder.c_str())});
  // m_res_phi_light_jet_eta2 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_lJets_eta2.txt", folder.c_str())});
  // m_res_phi_light_jet_eta3 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_lJets_eta3.txt", folder.c_str())});
  // m_res_phi_light_jet_eta4 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_lJets_eta4.txt", folder.c_str())});

  // m_res_phi_bjet_eta1 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_bJets_eta1.txt", folder.c_str())});
  // m_res_phi_bjet_eta2 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_bJets_eta2.txt", folder.c_str())});
  // m_res_phi_bjet_eta3 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_bJets_eta3.txt", folder.c_str())});
  // m_res_phi_bjet_eta4 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_bJets_eta4.txt", folder.c_str())});

  m_res_missing_ET = std::unique_ptr<ResolutionBase>(new ResGauss_MET{Form("%s/par_misset.txt", folder.c_str())});

  // default settings
  fResEnergyLightJet = m_res_energy_light_jet_eta1.get();
  fResEnergyBJet     = m_res_energy_bjet_eta1.get();
  fResEnergyGluonJet = m_res_energy_gluon_jet_eta1.get();
  fResEnergyElectron = m_res_energy_electron_eta1.get();
  fResEnergyMuon     = m_res_energy_muon_eta1.get();
  fResEnergyPhoton   = m_res_energy_photon_eta1.get();
  fResEtaLightJet    = m_res_eta_light_jet_eta1.get();
  fResEtaBJet        = m_res_eta_bjet_eta1.get();
  fResPhiLightJet    = m_res_phi_light_jet_eta1.get();
  fResPhiBJet        = m_res_phi_bjet_eta1.get();
  fResMissingET      = m_res_missing_ET.get();

  // Set eta binning for different objects starting with eta = 0
  m_jet_eta_bin_1 = 0.8;
  m_jet_eta_bin_2 = 1.37;
  m_jet_eta_bin_3 = 1.52;
  m_jet_eta_bin_4 = 2.50001;
  //  m_jet_eta_bin_5 = 4.5;

  m_electron_eta_bin_1 = 0.8;
  m_electron_eta_bin_2 = 1.37;
  m_electron_eta_bin_3 = 1.52;
  m_electron_eta_bin_4 = 2.50001;
}

// ---------------------------------------------------------
DetectorAtlas_8TeV::~DetectorAtlas_8TeV() = default;

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_8TeV::ResEnergyLightJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    fResEnergyLightJet = m_res_energy_light_jet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    fResEnergyLightJet = m_res_energy_light_jet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    fResEnergyLightJet = m_res_energy_light_jet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    fResEnergyLightJet = m_res_energy_light_jet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_8TeV::ResEnergyLightJet(). Eta range exceeded." << std::endl;
    return 0;
  }

  return fResEnergyLightJet;
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_8TeV::ResEnergyBJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    fResEnergyBJet = m_res_energy_bjet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    fResEnergyBJet = m_res_energy_bjet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    fResEnergyBJet = m_res_energy_bjet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    fResEnergyBJet = m_res_energy_bjet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_8TeV::ResEnergyBJet(). Eta range exceeded." << std::endl;
    return 0;
  }

  return fResEnergyBJet;
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_8TeV::ResEnergyGluonJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    fResEnergyGluonJet = m_res_energy_gluon_jet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    fResEnergyGluonJet = m_res_energy_gluon_jet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    fResEnergyGluonJet = m_res_energy_gluon_jet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    fResEnergyGluonJet = m_res_energy_gluon_jet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_8TeV::ResEnergyGluonJet(). Eta range exceeded." << std::endl;
    return 0;
  }

  return fResEnergyGluonJet;
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_8TeV::ResEnergyElectron(double eta) {
  if (fabs(eta) < m_electron_eta_bin_1) {
    fResEnergyElectron = m_res_energy_electron_eta1.get();
  } else if (fabs(eta) < m_electron_eta_bin_2) {
    fResEnergyElectron = m_res_energy_electron_eta2.get();
  } else if (fabs(eta) < m_electron_eta_bin_3) {
    std::cout << "DetectorAtlas_8TeV::ResEnergyElectron(). Electron in crack region" << std::endl;
    return 0;
  } else if (fabs(eta) <= m_electron_eta_bin_4) {
    fResEnergyElectron = m_res_energy_electron_eta4.get();
  } else {
    std::cout << "DetectorAtlas_8TeV::ResEnergyElectron(). Eta range exceeded." << std::endl;
    return 0;
  }

  return fResEnergyElectron;
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_8TeV::ResEnergyMuon(double eta) {
  if (fabs(eta) < m_muon_eta_bin_1) {
    fResEnergyMuon = m_res_energy_muon_eta1.get();
  } else if (fabs(eta) < m_muon_eta_bin_2) {
    fResEnergyMuon = m_res_energy_muon_eta2.get();
  } else if (fabs(eta) < m_muon_eta_bin_3) {
    fResEnergyMuon = m_res_energy_muon_eta3.get();
  } else {
    std::cout << "DetectorAtlas_8TeV::ResEnergyMuon(). Eta range exceeded." << std::endl;
    return 0;
  }

  return fResEnergyMuon;
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_8TeV::ResEnergyPhoton(double eta) {
  if (fabs(eta) < m_photon_eta_bin_1) {
    fResEnergyPhoton = m_res_energy_photon_eta1.get();
  } else if (fabs(eta) < m_photon_eta_bin_2) {
    fResEnergyPhoton = m_res_energy_photon_eta2.get();
  } else if (fabs(eta) < m_photon_eta_bin_3) {
    fResEnergyPhoton = m_res_energy_photon_eta3.get();
  } else if (fabs(eta) <= m_photon_eta_bin_4) {
    fResEnergyPhoton = m_res_energy_photon_eta4.get();
  } else {
    std::cout << "DetectorAtlas_8TeV::ResEnergyPhoton(). Eta range exceeded." << std::endl;
    return 0;
  }

  return fResEnergyPhoton;
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_8TeV::ResEtaLightJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    fResEtaLightJet = m_res_eta_light_jet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    fResEtaLightJet = m_res_eta_light_jet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    fResEtaLightJet = m_res_eta_light_jet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    fResEtaLightJet = m_res_eta_light_jet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_8TeV::ResEtaLightJet(). Eta range exceeded." << std::endl;
    return 0;
  }

  return fResEtaLightJet;
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_8TeV::ResEtaBJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    fResEtaBJet = m_res_eta_bjet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    fResEtaBJet = m_res_eta_bjet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    fResEtaBJet = m_res_eta_bjet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    fResEtaBJet = m_res_eta_bjet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_8TeV::ResEtaBJet(). Eta range exceeded." << std::endl;
    return 0;
  }

  return fResEtaBJet;
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_8TeV::ResPhiLightJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    fResPhiLightJet = m_res_phi_light_jet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    fResPhiLightJet = m_res_phi_light_jet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    fResPhiLightJet = m_res_phi_light_jet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    fResPhiLightJet = m_res_phi_light_jet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_8TeV::ResPhiLightJet(). Eta range exceeded." << std::endl;
    return 0;
  }

  return fResPhiLightJet;
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_8TeV::ResPhiBJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    fResPhiBJet = m_res_phi_bjet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    fResPhiBJet = m_res_phi_bjet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    fResPhiBJet = m_res_phi_bjet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    fResPhiLightJet = m_res_phi_light_jet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_8TeV::ResPhiBJet(). Eta range exceeded." << std::endl;
    return 0;
  }

  return fResPhiBJet;
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_8TeV::ResMissingET() {
  return fResMissingET;
}
}  // namespace KLFitter
