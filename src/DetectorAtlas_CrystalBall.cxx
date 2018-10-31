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

#include "KLFitter/DetectorAtlas_CrystalBall.h"

#include <cmath>
#include <iostream>

#include "KLFitter/ResolutionBase.h"
#include "KLFitter/ResSingleGaussE.h"
#include "KLFitter/ResSingleGaussMET.h"
#include "KLFitter/ResSingleGaussPt.h"
#include "KLFitter/ResCrystalBallJets.h"
#include "KLFitter/ResGauss.h"
#include "TString.h"

namespace KLFitter {
// ---------------------------------------------------------
DetectorAtlas_CrystalBall::DetectorAtlas_CrystalBall(std::string folder) : DetectorBase() {
  std::cout << "Using TF with crystal ball for jets ..." << std::endl;
  // energy resolution
  m_res_energy_light_jet_eta1 = std::unique_ptr<ResolutionBase>(new ResCrystalBallJets{Form("%s/par_energy_lJets_eta1.txt", folder.c_str())});
  m_res_energy_light_jet_eta2 = std::unique_ptr<ResolutionBase>(new ResCrystalBallJets{Form("%s/par_energy_lJets_eta2.txt", folder.c_str())});
  m_res_energy_light_jet_eta3 = std::unique_ptr<ResolutionBase>(new ResCrystalBallJets{Form("%s/par_energy_lJets_eta3.txt", folder.c_str())});
  m_res_energy_light_jet_eta4 = std::unique_ptr<ResolutionBase>(new ResCrystalBallJets{Form("%s/par_energy_lJets_eta4.txt", folder.c_str())});

  m_res_energy_bjet_eta1     = std::unique_ptr<ResolutionBase>(new ResCrystalBallJets{Form("%s/par_energy_bJets_eta1.txt", folder.c_str())});
  m_res_energy_bjet_eta2     = std::unique_ptr<ResolutionBase>(new ResCrystalBallJets{Form("%s/par_energy_bJets_eta2.txt", folder.c_str())});
  m_res_energy_bjet_eta3     = std::unique_ptr<ResolutionBase>(new ResCrystalBallJets{Form("%s/par_energy_bJets_eta3.txt", folder.c_str())});
  m_res_energy_bjet_eta4     = std::unique_ptr<ResolutionBase>(new ResCrystalBallJets{Form("%s/par_energy_bJets_eta4.txt", folder.c_str())});

  m_res_energy_gluon_jet_eta1 = std::unique_ptr<ResolutionBase>(new ResCrystalBallJets{Form("%s/par_energy_gluon_eta1.txt", folder.c_str())});
  m_res_energy_gluon_jet_eta2 = std::unique_ptr<ResolutionBase>(new ResCrystalBallJets{Form("%s/par_energy_gluon_eta2.txt", folder.c_str())});
  m_res_energy_gluon_jet_eta3 = std::unique_ptr<ResolutionBase>(new ResCrystalBallJets{Form("%s/par_energy_gluon_eta3.txt", folder.c_str())});
  m_res_energy_gluon_jet_eta4 = std::unique_ptr<ResolutionBase>(new ResCrystalBallJets{Form("%s/par_energy_gluon_eta4.txt", folder.c_str())});

  m_res_energy_electron_eta1 = std::unique_ptr<ResolutionBase>(new ResSingleGaussE{Form("%s/par_energy_Electrons_eta1.txt", folder.c_str())});
  m_res_energy_electron_eta2 = std::unique_ptr<ResolutionBase>(new ResSingleGaussE{Form("%s/par_energy_Electrons_eta2.txt", folder.c_str())});
  m_res_energy_electron_eta3 = std::unique_ptr<ResolutionBase>(new ResSingleGaussE{Form("%s/par_energy_Electrons_eta3.txt", folder.c_str())});
  m_res_energy_electron_eta4 = std::unique_ptr<ResolutionBase>(new ResSingleGaussE{Form("%s/par_energy_Electrons_eta4.txt", folder.c_str())});

  m_res_energy_muon_eta1     = std::unique_ptr<ResolutionBase>(new ResSingleGaussPt{Form("%s/par_energy_Muons_eta1.txt", folder.c_str())});
  m_res_energy_muon_eta2     = std::unique_ptr<ResolutionBase>(new ResSingleGaussPt{Form("%s/par_energy_Muons_eta2.txt", folder.c_str())});
  m_res_energy_muon_eta3     = std::unique_ptr<ResolutionBase>(new ResSingleGaussPt{Form("%s/par_energy_Muons_eta3.txt", folder.c_str())});

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

  m_res_missing_ET = std::unique_ptr<ResolutionBase>(new ResSingleGaussMET{Form("%s/par_misset.txt", folder.c_str())});
}

// ---------------------------------------------------------
DetectorAtlas_CrystalBall::~DetectorAtlas_CrystalBall() = default;

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_CrystalBall::ResEnergyLightJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    return m_res_energy_light_jet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    return m_res_energy_light_jet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    return m_res_energy_light_jet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    return m_res_energy_light_jet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_CrystalBall::ResEnergyLightJet(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_CrystalBall::ResEnergyBJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    return m_res_energy_bjet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    return m_res_energy_bjet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    return m_res_energy_bjet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    return m_res_energy_bjet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_CrystalBall::ResEnergyBJet(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_CrystalBall::ResEnergyGluonJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    return m_res_energy_gluon_jet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    return m_res_energy_gluon_jet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    return m_res_energy_gluon_jet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    return m_res_energy_gluon_jet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_CrystalBall::ResEnergyGluonJet(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_CrystalBall::ResEnergyElectron(double eta) {
  if (fabs(eta) < m_electron_eta_bin_1) {
    return m_res_energy_electron_eta1.get();
  } else if (fabs(eta) < m_electron_eta_bin_2) {
    return m_res_energy_electron_eta2.get();
  } else if (fabs(eta) < m_electron_eta_bin_3) {
    std::cout << "DetectorAtlas_CrystalBall::ResEnergyElectron(). Electron in crack region" << std::endl;
    return nullptr;
  } else if (fabs(eta) <= m_electron_eta_bin_4) {
    return m_res_energy_electron_eta4.get();
  } else {
    std::cout << "DetectorAtlas_CrystalBall::ResEnergyElectron(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_CrystalBall::ResEnergyMuon(double eta) {
  if (fabs(eta) < m_muon_eta_bin_1) {
    return m_res_energy_muon_eta1.get();
  } else if (fabs(eta) < m_muon_eta_bin_2) {
    return m_res_energy_muon_eta2.get();
  } else if (fabs(eta) < m_muon_eta_bin_3) {
    return m_res_energy_muon_eta3.get();
  } else {
    std::cout << "DetectorAtlas_CrystalBall::ResEnergyMuon(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_CrystalBall::ResEnergyPhoton(double eta) {
  if (fabs(eta) < m_photon_eta_bin_1) {
    return m_res_energy_photon_eta1.get();
  } else if (fabs(eta) < m_photon_eta_bin_2) {
    return m_res_energy_photon_eta2.get();
  } else if (fabs(eta) < m_photon_eta_bin_3) {
    return m_res_energy_photon_eta3.get();
  } else if (fabs(eta) <= m_photon_eta_bin_4) {
    return m_res_energy_photon_eta4.get();
  } else {
    std::cout << "DetectorAtlas_CrystalBall::ResEnergyPhoton(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_CrystalBall::ResEtaLightJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    return m_res_eta_light_jet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    return m_res_eta_light_jet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    return m_res_eta_light_jet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    return m_res_eta_light_jet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_CrystalBall::ResEtaLightJet(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_CrystalBall::ResEtaBJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    return m_res_eta_bjet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    return m_res_eta_bjet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    return m_res_eta_bjet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    return m_res_eta_bjet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_CrystalBall::ResEtaBJet(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_CrystalBall::ResPhiLightJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    return m_res_phi_light_jet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    return m_res_phi_light_jet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    return m_res_phi_light_jet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    return m_res_phi_light_jet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_CrystalBall::ResPhiLightJet(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_CrystalBall::ResPhiBJet(double eta) {
  if (fabs(eta) < m_jet_eta_bin_1) {
    return m_res_phi_bjet_eta1.get();
  } else if (fabs(eta) < m_jet_eta_bin_2) {
    return m_res_phi_bjet_eta2.get();
  } else if (fabs(eta) < m_jet_eta_bin_3) {
    return m_res_phi_bjet_eta3.get();
  } else if (fabs(eta) <= m_jet_eta_bin_4) {
    return m_res_phi_light_jet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_CrystalBall::ResPhiBJet(). Eta range exceeded." << std::endl;
    return nullptr;
  }
}

// ---------------------------------------------------------
ResolutionBase* DetectorAtlas_CrystalBall::ResMissingET() {
  return m_res_missing_ET.get();
}
}  // namespace KLFitter
