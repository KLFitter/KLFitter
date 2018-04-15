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
// --------------------------------------------------------- 
KLFitter::DetectorSnowmass::DetectorSnowmass(std::string folder) : DetectorBase() 
{
  std::cout << "Using TFs from SnowMass ..." << std::endl;
  fResEnergyJet_eta1      = new KLFitter::ResGaussE(Form("%s/par_energy_jets_eta1.txt", folder.c_str()));
  fResEnergyElectron_eta1 = new KLFitter::ResGaussE(Form("%s/par_energy_electrons_eta1.txt", folder.c_str()));
  fResMomentumMuon_eta1   = new KLFitter::ResGaussPt(Form("%s/par_pt_muons_eta1.txt", folder.c_str()));
  //fResMomentumMuon_eta1   = new KLFitter::ResGaussE(Form("%s/par_energy_electrons_eta1.txt", folder.c_str()));
  
  fResEnergyJet_eta2      = new KLFitter::ResGaussE(Form("%s/par_energy_jets_eta2.txt", folder.c_str()));
  fResEnergyElectron_eta2 = new KLFitter::ResGaussE(Form("%s/par_energy_electrons_eta2.txt", folder.c_str()));
  fResMomentumMuon_eta2   = new KLFitter::ResGaussPt(Form("%s/par_pt_muons_eta2.txt", folder.c_str()));
  //fResMomentumMuon_eta2   = new KLFitter::ResGaussE(Form("%s/par_energy_electrons_eta2.txt", folder.c_str()));
  
  fResEnergyJet_eta3      = new KLFitter::ResGaussE(Form("%s/par_energy_jets_eta3.txt", folder.c_str()));
  
  // missing et resolution in x and y
  fResMissingET      = new KLFitter::ResGauss_MET(Form("%s/par_misset.txt", folder.c_str()));
  
  // default settings
  fResEnergyLightJet = fResEnergyJet_eta1;
  fResEnergyBJet     = fResEnergyJet_eta1;
  fResEnergyGluonJet = fResEnergyJet_eta1;
  fResEnergyElectron = fResEnergyElectron_eta1;
  fResEnergyPhoton   = fResEnergyElectron_eta1;
  fResEnergyMuon     = fResMomentumMuon_eta1;

  // Set eta binning for different objects starting with eta = 0
  fJetEtaBin_1 = 1.7;
  fJetEtaBin_2 = 3.2;
  fJetEtaBin_3 = 4.9;

  fElectronEtaBin_1 = 3.0;
  fElectronEtaBin_2 = 5.0;

  fMuonEtaBin_1 = 1.5;
  fMuonEtaBin_2 = 2.5;

}
// --------------------------------------------------------- 
KLFitter::DetectorSnowmass::~DetectorSnowmass() 
{
  if (fResEnergyJet_eta1)
    delete fResEnergyJet_eta1;
  
  if (fResEnergyElectron_eta1)
    delete fResEnergyElectron_eta1;

  if (fResMomentumMuon_eta1)
    delete fResMomentumMuon_eta1;
  
  if (fResEnergyJet_eta2)
    delete fResEnergyJet_eta2;
  
  if (fResEnergyElectron_eta2)
    delete fResEnergyElectron_eta2;

  if (fResMomentumMuon_eta2)
    delete fResMomentumMuon_eta2;
  
  if (fResEnergyJet_eta3)
    delete fResEnergyJet_eta3; 

}
// ---------------------------------------------------------
KLFitter::ResolutionBase * KLFitter::DetectorSnowmass::ResEnergyLightJet(double eta) {
  if (fabs(eta) < fJetEtaBin_1) {
    fResEnergyLightJet = fResEnergyJet_eta1;
  } else if (fabs(eta) < fJetEtaBin_2) {
    fResEnergyLightJet = fResEnergyJet_eta2;
  } else if (fabs(eta) < fJetEtaBin_3) {
    fResEnergyLightJet = fResEnergyJet_eta3;
  } else {
    std::cout << "KLFitter::DetectorSnowmass::ResEnergyLightJet(). Eta range exceeded." << std::endl;
    return 0;
  }

  return fResEnergyLightJet;
}
// ---------------------------------------------------------
KLFitter::ResolutionBase * KLFitter::DetectorSnowmass::ResEnergyBJet(double eta) {
  if (fabs(eta) < fJetEtaBin_1) {
    fResEnergyBJet = fResEnergyJet_eta1;
  } else if (fabs(eta) < fJetEtaBin_2) {
    fResEnergyBJet = fResEnergyJet_eta2;
  } else if (fabs(eta) < fJetEtaBin_3) {
    fResEnergyBJet = fResEnergyJet_eta3;
  } else {
    std::cout << "KLFitter::DetectorSnowmass::ResEnergyBJet(). Eta range exceeded." << std::endl;
    return 0;
  }

  return fResEnergyBJet;
}
// ---------------------------------------------------------
KLFitter::ResolutionBase * KLFitter::DetectorSnowmass::ResEnergyElectron(double eta) {
  if (fabs(eta) < fElectronEtaBin_1) {
    fResEnergyElectron = fResEnergyElectron_eta1;
  } else if (fabs(eta) < fElectronEtaBin_2) {
    fResEnergyElectron = fResEnergyElectron_eta2;
  } else {
    std::cout << "KLFitter::DetectorSnowmass::ResEnergyElectron(). Eta range exceeded." << std::endl;
    return 0;
  }
  return fResEnergyElectron;
}
// ---------------------------------------------------------
KLFitter::ResolutionBase * KLFitter::DetectorSnowmass::ResEnergyMuon(double eta) {
  if (fabs(eta) < fMuonEtaBin_1) {
    fResEnergyMuon = fResMomentumMuon_eta1;
  } else if (fabs(eta) < fMuonEtaBin_2) {
    fResEnergyMuon = fResMomentumMuon_eta2;
  } else {
    std::cout << "KLFitter::DetectorSnowmass::ResEnergyMuon(). Eta range exceeded." << std::endl;
    return 0;
  }

  return fResEnergyMuon;
}
// ---------------------------------------------------------
KLFitter::ResolutionBase * KLFitter::DetectorSnowmass::ResMissingET() {
  return fResMissingET;
}
// ---------------------------------------------------------
