#include "DetectorAtlas_10TeV.h" 
#include "ResolutionBase.h" 
#include "ResDoubleGaussPt.h" 
#include "ResDoubleGaussE_1.h"
#include "ResDoubleGaussE_2.h" 
#include "ResDoubleGaussE_3.h" 
#include "ResGauss.h"

#include <cmath>
#include <iostream>

// --------------------------------------------------------- 

KLFitter::DetectorAtlas_10TeV::DetectorAtlas_10TeV() : DetectorBase() 
{
  // energy resolution 
  fResEnergyLightJet_eta1 = new KLFitter::ResDoubleGaussE_3("../../transferfunctions/ttbar/par_energy_light_eta1_10TeV.txt"); 
  fResEnergyBJet_eta1     = new KLFitter::ResDoubleGaussE_3("../../transferfunctions/ttbar/par_energy_b_eta1_10TeV.txt"); 
  fResEnergyGluonJet_eta1 = new KLFitter::ResDoubleGaussE_3("../../transferfunctions/ttbar/par_energy_gluon_eta1_10TeV.txt"); 
  fResEnergyElectron_eta1 = new KLFitter::ResDoubleGaussE_3("../../transferfunctions/ttbar/par_energy_electron_eta1_10TeV.txt"); 
  fResEnergyMuon_eta1     = new KLFitter::ResGauss("../../transferfunctions/ttbar/par_energy_muon_eta1_10TeV.txt"); 
  fResEnergyPhoton_eta1   = new KLFitter::ResGauss("../../transferfunctions/ttbar/par_energy_photon_eta1_10TeV.txt"); 

  fResEnergyLightJet_eta2 = new KLFitter::ResDoubleGaussE_3("../../transferfunctions/ttbar/par_energy_light_eta2_10TeV.txt"); 
  fResEnergyBJet_eta2     = new KLFitter::ResDoubleGaussE_2("../../transferfunctions/ttbar/par_energy_b_eta2_10TeV.txt"); 
  fResEnergyGluonJet_eta2 = new KLFitter::ResDoubleGaussE_3("../../transferfunctions/ttbar/par_energy_gluon_eta2_10TeV.txt"); 
  fResEnergyElectron_eta2 = new KLFitter::ResDoubleGaussE_3("../../transferfunctions/ttbar/par_energy_electron_eta2_10TeV.txt"); 
  fResEnergyMuon_eta2     = new KLFitter::ResGauss("../../transferfunctions/ttbar/par_energy_muon_eta2_10TeV.txt"); 
  fResEnergyPhoton_eta2   = new KLFitter::ResGauss("../../transferfunctions/ttbar/par_energy_photon_eta2_10TeV.txt"); 

  fResEnergyLightJet_eta3 = new KLFitter::ResDoubleGaussE_3("../../transferfunctions/ttbar/par_energy_light_eta3_10TeV.txt"); 
  fResEnergyBJet_eta3     = new KLFitter::ResDoubleGaussE_2("../../transferfunctions/ttbar/par_energy_b_eta3_10TeV.txt"); 
  fResEnergyGluonJet_eta3 = new KLFitter::ResDoubleGaussE_3("../../transferfunctions/ttbar/par_energy_gluon_eta3_10TeV.txt"); 
  fResEnergyElectron_eta3 = new KLFitter::ResDoubleGaussE_3("../../transferfunctions/ttbar/par_energy_electron_eta3_10TeV.txt"); 
  fResEnergyMuon_eta3     = new KLFitter::ResGauss("../../transferfunctions/ttbar/par_energy_muon_eta3_10TeV.txt"); 
  fResEnergyPhoton_eta3   = new KLFitter::ResGauss("../../transferfunctions/ttbar/par_energy_photon_eta3_10TeV.txt"); 

  // eta resolution 
  fResEtaLightJet_eta1 = new KLFitter::ResGauss("../../transferfunctions/ttbar/par_eta_light_eta1_10TeV.txt"); 
  fResEtaLightJet_eta2 = new KLFitter::ResGauss("../../transferfunctions/ttbar/par_eta_light_eta2_10TeV.txt"); 
  fResEtaLightJet_eta3 = new KLFitter::ResGauss("../../transferfunctions/ttbar/par_eta_light_eta3_10TeV.txt"); 

  fResEtaBJet_eta1 = new KLFitter::ResGauss("../../transferfunctions/ttbar/par_eta_b_eta1_10TeV.txt"); 
  fResEtaBJet_eta2 = new KLFitter::ResGauss("../../transferfunctions/ttbar/par_eta_b_eta2_10TeV.txt"); 
  fResEtaBJet_eta3 = new KLFitter::ResGauss("../../transferfunctions/ttbar/par_eta_b_eta3_10TeV.txt"); 

  // phi resolution 
  fResPhiLightJet_eta1 = new KLFitter::ResGauss("../../transferfunctions/ttbar/par_phi_light_eta1_10TeV.txt"); 
  fResPhiLightJet_eta2 = new KLFitter::ResGauss("../../transferfunctions/ttbar/par_phi_light_eta2_10TeV.txt"); 
  fResPhiLightJet_eta3 = new KLFitter::ResGauss("../../transferfunctions/ttbar/par_phi_light_eta3_10TeV.txt"); 

  fResPhiBJet_eta1 = new KLFitter::ResGauss("../../transferfunctions/ttbar/par_phi_b_eta1_10TeV.txt"); 
  fResPhiBJet_eta2 = new KLFitter::ResGauss("../../transferfunctions/ttbar/par_phi_b_eta2_10TeV.txt"); 
  fResPhiBJet_eta3 = new KLFitter::ResGauss("../../transferfunctions/ttbar/par_phi_b_eta3_10TeV.txt"); 

  // missing et resolution in x and y 
  fResMissingET      = new KLFitter::ResGauss("../../transferfunctions/ttbar/par_misset_10TeV.txt"); 

  // default settings 
  fResEnergyLightJet = fResEnergyLightJet_eta1; 
  fResEnergyBJet     = fResEnergyBJet_eta1; 
  fResEnergyGluonJet = fResEnergyGluonJet_eta1; 
  fResEnergyElectron = fResEnergyElectron_eta1; 
  fResEnergyMuon     = fResEnergyMuon_eta1; 
  fResEnergyPhoton   = fResEnergyPhoton_eta1; 
  fResEtaLightJet    = fResEtaLightJet_eta1; 
  fResEtaBJet        = fResEtaBJet_eta1; 
  fResPhiLightJet    = fResPhiLightJet_eta1; 
  fResPhiBJet        = fResPhiBJet_eta1; 
}

// --------------------------------------------------------- 

KLFitter::DetectorAtlas_10TeV::~DetectorAtlas_10TeV() 
{
  if (fResEnergyLightJet_eta1) 
    delete fResEnergyLightJet_eta1; 
        
  if (fResEnergyLightJet_eta2)
    delete fResEnergyLightJet_eta2; 

  if (fResEnergyLightJet_eta3)
    delete fResEnergyLightJet_eta3; 

  if (fResEnergyBJet_eta1)
    delete fResEnergyBJet_eta1; 

  if (fResEnergyBJet_eta2)
    delete fResEnergyBJet_eta2; 

  if (fResEnergyBJet_eta3)
    delete fResEnergyBJet_eta3; 

  if (fResEnergyGluonJet_eta1)
    delete fResEnergyGluonJet_eta1; 

  if (fResEnergyGluonJet_eta2)
    delete fResEnergyGluonJet_eta2; 

  if (fResEnergyGluonJet_eta3)
    delete fResEnergyGluonJet_eta3; 

  if (fResEnergyElectron_eta1)
    delete fResEnergyElectron_eta1; 

  if (fResEnergyElectron_eta2)
    delete fResEnergyElectron_eta2; 

  if (fResEnergyElectron_eta3)
    delete fResEnergyElectron_eta3; 

  if (fResEnergyMuon_eta1)
    delete fResEnergyMuon_eta1; 

  if (fResEnergyMuon_eta2)
    delete fResEnergyMuon_eta2; 

  if (fResEnergyMuon_eta3)
    delete fResEnergyMuon_eta3; 

  if (fResEnergyPhoton_eta1)
    delete fResEnergyPhoton_eta1; 

  if (fResEnergyPhoton_eta2)
    delete fResEnergyPhoton_eta2; 

  if (fResEnergyPhoton_eta3)
    delete fResEnergyPhoton_eta3; 

  if (fResEtaLightJet_eta1) 
    delete fResEtaLightJet_eta1; 

  if (fResEtaLightJet_eta2) 
    delete fResEtaLightJet_eta2; 

  if (fResEtaLightJet_eta3) 
    delete fResEtaLightJet_eta3; 
        
  if (fResEtaBJet_eta1) 
    delete fResEtaBJet_eta1; 

  if (fResEtaBJet_eta2) 
    delete fResEtaBJet_eta2; 

  if (fResEtaBJet_eta3) 
    delete fResEtaBJet_eta3; 
        
  if (fResPhiLightJet_eta1) 
    delete fResPhiLightJet_eta1; 

  if (fResPhiLightJet_eta2) 
    delete fResPhiLightJet_eta2; 

  if (fResPhiLightJet_eta3) 
    delete fResPhiLightJet_eta3; 

  if (fResPhiBJet_eta1) 
    delete fResPhiBJet_eta1; 

  if (fResPhiBJet_eta2) 
    delete fResPhiBJet_eta2; 

  if (fResPhiBJet_eta3) 
    delete fResPhiBJet_eta3; 
        
  if (fResMissingET)
    delete fResMissingET; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_10TeV::ResEnergyLightJet(double eta)
{
  if (fabs(eta) < 1.0) 
    fResEnergyLightJet = fResEnergyLightJet_eta1; 
  else if (fabs(eta) < 1.7) 
    fResEnergyLightJet = fResEnergyLightJet_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResEnergyLightJet = fResEnergyLightJet_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_10TeV::ResEnergyLightJet(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEnergyLightJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_10TeV::ResEnergyBJet(double eta)
{
  if (fabs(eta) < 1.0) 
    fResEnergyBJet = fResEnergyBJet_eta1; 
  else if (fabs(eta) < 1.7) 
    fResEnergyBJet = fResEnergyBJet_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResEnergyBJet = fResEnergyBJet_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_10TeV::ResEnergyBJet(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEnergyBJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_10TeV::ResEnergyGluonJet(double eta)
{
  if (fabs(eta) < 1.0) 
    fResEnergyGluonJet = fResEnergyGluonJet_eta1; 
  else if (fabs(eta) < 1.7) 
    fResEnergyGluonJet = fResEnergyGluonJet_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResEnergyGluonJet = fResEnergyGluonJet_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_10TeV::ResEnergyGluonJet(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEnergyGluonJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_10TeV::ResEnergyElectron(double eta)
{
  if (fabs(eta) < 1.0) 
    fResEnergyElectron = fResEnergyElectron_eta1; 
  else if (fabs(eta) < 1.7) 
    fResEnergyElectron = fResEnergyElectron_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResEnergyElectron = fResEnergyElectron_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_10TeV::ResEnergyElectron(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEnergyElectron; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_10TeV::ResEnergyMuon(double eta)
{
  if (fabs(eta) < 1.0) 
    fResEnergyMuon = fResEnergyMuon_eta1; 
  else if (fabs(eta) < 1.7) 
    fResEnergyMuon = fResEnergyMuon_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResEnergyMuon = fResEnergyMuon_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_10TeV::ResEnergyMuon(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEnergyMuon; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_10TeV::ResEnergyPhoton(double eta)
{
  if (fabs(eta) < 1.0) 
    fResEnergyPhoton = fResEnergyPhoton_eta1; 
  else if (fabs(eta) < 1.7) 
    fResEnergyPhoton = fResEnergyPhoton_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResEnergyPhoton = fResEnergyPhoton_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_10TeV::ResEnergyPhoton(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEnergyPhoton; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_10TeV::ResEtaLightJet(double eta)
{
  if (fabs(eta) < 1.0) 
    fResEtaLightJet = fResEtaLightJet_eta1; 
  else if (fabs(eta) < 1.7) 
    fResEtaLightJet = fResEtaLightJet_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResEtaLightJet = fResEtaLightJet_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_10TeV::ResEtaLightJet(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEtaLightJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_10TeV::ResEtaBJet(double eta)
{
  if (fabs(eta) < 1.0) 
    fResEtaBJet = fResEtaBJet_eta1; 
  else if (fabs(eta) < 1.7) 
    fResEtaBJet = fResEtaBJet_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResEtaBJet = fResEtaBJet_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_10TeV::ResEtaBJet(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEtaBJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_10TeV::ResPhiLightJet(double eta)
{
  if (fabs(eta) < 1.0) 
    fResPhiLightJet = fResPhiLightJet_eta1; 
  else if (fabs(eta) < 1.7) 
    fResPhiLightJet = fResPhiLightJet_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResPhiLightJet = fResPhiLightJet_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_10TeV::ResPhiLightJet(). Phi range exceeded." << std::endl; 
      return 0; 
    }

  return fResPhiLightJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_10TeV::ResPhiBJet(double eta)
{
  if (fabs(eta) < 1.0) 
    fResPhiBJet = fResPhiBJet_eta1; 
  else if (fabs(eta) < 1.7) 
    fResPhiBJet = fResPhiBJet_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResPhiBJet = fResPhiBJet_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_10TeV::ResPhiBJet(). Phi range exceeded." << std::endl; 
      return 0; 
    }

  return fResPhiBJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_10TeV::ResMissingET()
{
  return fResMissingET; 
}

// --------------------------------------------------------- 
