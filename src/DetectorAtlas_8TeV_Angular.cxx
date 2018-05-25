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

#include "KLFitter/DetectorAtlas_8TeV_Angular.h"

#include <iostream>

#include "KLFitter/ResGauss.h"
#include "KLFitter/ResolutionBase.h"
#include "TString.h"

namespace KLFitter {
// ---------------------------------------------------------
DetectorAtlas_8TeV_Angular::DetectorAtlas_8TeV_Angular(std::string folder) : DetectorAtlas_8TeV::DetectorAtlas_8TeV(folder) {
  // eta resolution
  fResEtaLightJet_eta1 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_lJets_eta1.txt", folder.c_str())});
  fResEtaLightJet_eta2 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_lJets_eta2.txt", folder.c_str())});
  fResEtaLightJet_eta3 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_lJets_eta3.txt", folder.c_str())});
  fResEtaLightJet_eta4 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_lJets_eta4.txt", folder.c_str())});

  fResEtaBJet_eta1 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_bJets_eta1.txt", folder.c_str())});
  fResEtaBJet_eta2 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_bJets_eta2.txt", folder.c_str())});
  fResEtaBJet_eta3 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_bJets_eta3.txt", folder.c_str())});
  fResEtaBJet_eta4 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_eta_bJets_eta4.txt", folder.c_str())});

  // phi resolution
  fResPhiLightJet_eta1 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_lJets_eta1.txt", folder.c_str())});
  fResPhiLightJet_eta2 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_lJets_eta2.txt", folder.c_str())});
  fResPhiLightJet_eta3 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_lJets_eta3.txt", folder.c_str())});
  fResPhiLightJet_eta4 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_lJets_eta4.txt", folder.c_str())});

  fResPhiBJet_eta1 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_bJets_eta1.txt", folder.c_str())});
  fResPhiBJet_eta2 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_bJets_eta2.txt", folder.c_str())});
  fResPhiBJet_eta3 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_bJets_eta3.txt", folder.c_str())});
  fResPhiBJet_eta4 = std::unique_ptr<ResolutionBase>(new ResGauss{Form("%s/par_phi_bJets_eta4.txt", folder.c_str())});

  // default settings
  fResEtaLightJet    = fResEtaLightJet_eta1.get();
  fResEtaBJet        = fResEtaBJet_eta1.get();
  fResPhiLightJet    = fResPhiLightJet_eta1.get();
  fResPhiBJet        = fResPhiBJet_eta1.get();
}

// ---------------------------------------------------------
DetectorAtlas_8TeV_Angular::~DetectorAtlas_8TeV_Angular() = default;

// ---------------------------------------------------------
ResolutionBase * DetectorAtlas_8TeV_Angular::ResEtaLightJet(double eta) {
  if (fabs(eta) < fJetEtaBin_1) {
    fResEtaLightJet = fResEtaLightJet_eta1.get();
  } else if (fabs(eta) < fJetEtaBin_2) {
    fResEtaLightJet = fResEtaLightJet_eta2.get();
  } else if (fabs(eta) < fJetEtaBin_3) {
    fResEtaLightJet = fResEtaLightJet_eta3.get();
  } else if (fabs(eta) <= fJetEtaBin_4) {
    fResEtaLightJet = fResEtaLightJet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_8TeV_Angular::ResEtaLightJet(). Eta range exceeded." << std::endl;
    return 0;
  }

  return fResEtaLightJet;
}

// ---------------------------------------------------------
ResolutionBase * DetectorAtlas_8TeV_Angular::ResEtaBJet(double eta) {
  if (fabs(eta) < fJetEtaBin_1) {
    fResEtaBJet = fResEtaBJet_eta1.get();
  } else if (fabs(eta) < fJetEtaBin_2) {
    fResEtaBJet = fResEtaBJet_eta2.get();
  } else if (fabs(eta) < fJetEtaBin_3) {
    fResEtaBJet = fResEtaBJet_eta3.get();
  } else if (fabs(eta) <= fJetEtaBin_4) {
    fResEtaBJet = fResEtaBJet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_8TeV_Angular::ResEtaBJet(). Eta range exceeded." << std::endl;
    return 0;
  }

  return fResEtaBJet;
}

// ---------------------------------------------------------
ResolutionBase * DetectorAtlas_8TeV_Angular::ResPhiLightJet(double eta) {
  if (fabs(eta) < fJetEtaBin_1) {
    fResPhiLightJet = fResPhiLightJet_eta1.get();
  } else if (fabs(eta) < fJetEtaBin_2) {
    fResPhiLightJet = fResPhiLightJet_eta2.get();
  } else if (fabs(eta) < fJetEtaBin_3) {
    fResPhiLightJet = fResPhiLightJet_eta3.get();
  } else if (fabs(eta) <= fJetEtaBin_4) {
    fResPhiLightJet = fResPhiLightJet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_8TeV_Angular::ResPhiLightJet(). Eta range exceeded." << std::endl;
    return 0;
  }

  return fResPhiLightJet;
}

// ---------------------------------------------------------
ResolutionBase * DetectorAtlas_8TeV_Angular::ResPhiBJet(double eta) {
  if (fabs(eta) < fJetEtaBin_1) {
    fResPhiBJet = fResPhiBJet_eta1.get();
  } else if (fabs(eta) < fJetEtaBin_2) {
    fResPhiBJet = fResPhiBJet_eta2.get();
  } else if (fabs(eta) < fJetEtaBin_3) {
    fResPhiBJet = fResPhiBJet_eta3.get();
  } else if (fabs(eta) <= fJetEtaBin_4) {
    fResPhiLightJet = fResPhiLightJet_eta4.get();
  } else {
    std::cout << "DetectorAtlas_8TeV_Angular::ResPhiBJet(). Eta range exceeded." << std::endl;
    return 0;
  }

  return fResPhiBJet;
}
}   // namespace KLFitter
