#include "DetectorAtlas.h" 
#include "ResolutionBase.h" 
#include "ResDoubleGaussE.h" 
#include "ResGauss.h"

#include <cmath>
#include <iostream>

// --------------------------------------------------------- 

KLFitter::DetectorAtlas::DetectorAtlas() : DetectorBase() 
{
  // energy resolution 
  fResEnergyLightJet_eta1 = new KLFitter::ResDoubleGaussE("par_energy_light_eta1.txt"); 
  fResEnergyBJet_eta1     = new KLFitter::ResDoubleGaussE("par_energy_b_eta1.txt"); 
  fResEnergyGluonJet_eta1 = new KLFitter::ResDoubleGaussE("par_energy_gluon_eta1.txt"); 
  fResEnergyElectron_eta1 = new KLFitter::ResDoubleGaussE("par_energy_electron_eta1.txt"); 
  //    fResEnergyMuon_eta1     = new KLFitter::ResDoubleGaussPt("par_energy_muon_eta1.txt"); 
  fResEnergyMuon_eta1     = new KLFitter::ResGauss("par_energy_muon_eta1.txt"); 
  fResEnergyPhoton_eta1   = new KLFitter::ResGauss("par_energy_photon_eta1.txt"); 

  fResEnergyLightJet_eta2 = new KLFitter::ResDoubleGaussE("par_energy_light_eta2.txt"); 
  fResEnergyBJet_eta2     = new KLFitter::ResDoubleGaussE("par_energy_b_eta2.txt"); 
  fResEnergyGluonJet_eta2 = new KLFitter::ResDoubleGaussE("par_energy_gluon_eta2.txt"); 
  fResEnergyElectron_eta2 = new KLFitter::ResDoubleGaussE("par_energy_electron_eta2.txt"); 
  //    fResEnergyMuon_eta2     = new KLFitter::ResDoubleGaussPt("par_energy_muon_eta2.txt"); 
  fResEnergyMuon_eta2     = new KLFitter::ResGauss("par_energy_muon_eta2.txt"); 
  fResEnergyPhoton_eta2   = new KLFitter::ResGauss("par_energy_photon_eta2.txt"); 

  fResEnergyLightJet_eta3 = new KLFitter::ResDoubleGaussE("par_energy_light_eta3.txt"); 
  fResEnergyBJet_eta3     = new KLFitter::ResDoubleGaussE("par_energy_b_eta3.txt"); 
  fResEnergyGluonJet_eta3 = new KLFitter::ResDoubleGaussE("par_energy_gluon_eta3.txt"); 
  fResEnergyElectron_eta3 = new KLFitter::ResDoubleGaussE("par_energy_electron_eta3.txt"); 
  //    fResEnergyMuon_eta3     = new KLFitter::ResDoubleGaussPt("par_energy_muon_eta3.txt"); 
  fResEnergyMuon_eta3     = new KLFitter::ResGauss("par_energy_muon_eta3.txt"); 
  fResEnergyPhoton_eta3   = new KLFitter::ResGauss("par_energy_photon_eta3.txt"); 

  // eta resolution 
  fResEtaLightJet_eta1 = new KLFitter::ResGauss("par_eta_light_eta1.txt"); 
  fResEtaLightJet_eta2 = new KLFitter::ResGauss("par_eta_light_eta2.txt"); 
  fResEtaLightJet_eta3 = new KLFitter::ResGauss("par_eta_light_eta3.txt"); 

  fResEtaBJet_eta1 = new KLFitter::ResGauss("par_eta_b_eta1.txt"); 
  fResEtaBJet_eta2 = new KLFitter::ResGauss("par_eta_b_eta2.txt"); 
  fResEtaBJet_eta3 = new KLFitter::ResGauss("par_eta_b_eta3.txt"); 

  // phi resolution 
  fResPhiLightJet_eta1 = new KLFitter::ResGauss("par_phi_light_eta1.txt"); 
  fResPhiLightJet_eta2 = new KLFitter::ResGauss("par_phi_light_eta2.txt"); 
  fResPhiLightJet_eta3 = new KLFitter::ResGauss("par_phi_light_eta3.txt"); 

  fResPhiBJet_eta1 = new KLFitter::ResGauss("par_phi_b_eta1.txt"); 
  fResPhiBJet_eta2 = new KLFitter::ResGauss("par_phi_b_eta2.txt"); 
  fResPhiBJet_eta3 = new KLFitter::ResGauss("par_phi_b_eta3.txt"); 

  // missing et resolution in x and y 
  fResMissingET      = new KLFitter::ResGauss("par_misset.txt"); 

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

KLFitter::DetectorAtlas::~DetectorAtlas() 
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
KLFitter::ResolutionBase * KLFitter::DetectorAtlas::ResEnergyLightJet(double eta)
{
  if (fabs(eta) < 1.0) 
    fResEnergyLightJet = fResEnergyLightJet_eta1; 
  else if (fabs(eta) < 1.7) 
    fResEnergyLightJet = fResEnergyLightJet_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResEnergyLightJet = fResEnergyLightJet_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorATLAS::ResEnergyLightJet(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEnergyLightJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas::ResEnergyBJet(double eta)
{
  if (fabs(eta) < 1.0) 
    fResEnergyBJet = fResEnergyBJet_eta1; 
  else if (fabs(eta) < 1.7) 
    fResEnergyBJet = fResEnergyBJet_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResEnergyBJet = fResEnergyBJet_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorATLAS::ResEnergyBJet(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEnergyBJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas::ResEnergyGluonJet(double eta)
{
  if (fabs(eta) < 1.0) 
    fResEnergyGluonJet = fResEnergyGluonJet_eta1; 
  else if (fabs(eta) < 1.7) 
    fResEnergyGluonJet = fResEnergyGluonJet_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResEnergyGluonJet = fResEnergyGluonJet_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorATLAS::ResEnergyGluonJet(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEnergyGluonJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas::ResEnergyElectron(double eta)
{
  if (fabs(eta) < 1.0) 
    fResEnergyElectron = fResEnergyElectron_eta1; 
  else if (fabs(eta) < 1.7) 
    fResEnergyElectron = fResEnergyElectron_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResEnergyElectron = fResEnergyElectron_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorATLAS::ResEnergyElectron(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEnergyElectron; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas::ResEnergyMuon(double eta)
{
  if (fabs(eta) < 1.0) 
    fResEnergyMuon = fResEnergyMuon_eta1; 
  else if (fabs(eta) < 1.7) 
    fResEnergyMuon = fResEnergyMuon_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResEnergyMuon = fResEnergyMuon_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorATLAS::ResEnergyMuon(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEnergyMuon; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas::ResEnergyPhoton(double eta)
{
  if (fabs(eta) < 1.0) 
    fResEnergyPhoton = fResEnergyPhoton_eta1; 
  else if (fabs(eta) < 1.7) 
    fResEnergyPhoton = fResEnergyPhoton_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResEnergyPhoton = fResEnergyPhoton_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorATLAS::ResEnergyPhoton(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEnergyPhoton; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas::ResEtaLightJet(double eta)
{
  if (fabs(eta) < 1.0) 
    fResEtaLightJet = fResEtaLightJet_eta1; 
  else if (fabs(eta) < 1.7) 
    fResEtaLightJet = fResEtaLightJet_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResEtaLightJet = fResEtaLightJet_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorATLAS::ResEtaLightJet(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEtaLightJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas::ResEtaBJet(double eta)
{
  if (fabs(eta) < 1.0) 
    fResEtaBJet = fResEtaBJet_eta1; 
  else if (fabs(eta) < 1.7) 
    fResEtaBJet = fResEtaBJet_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResEtaBJet = fResEtaBJet_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorATLAS::ResEtaBJet(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEtaBJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas::ResPhiLightJet(double eta)
{
  if (fabs(eta) < 1.0) 
    fResPhiLightJet = fResPhiLightJet_eta1; 
  else if (fabs(eta) < 1.7) 
    fResPhiLightJet = fResPhiLightJet_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResPhiLightJet = fResPhiLightJet_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorATLAS::ResPhiLightJet(). Phi range exceeded." << std::endl; 
      return 0; 
    }

  return fResPhiLightJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas::ResPhiBJet(double eta)
{
  if (fabs(eta) < 1.0) 
    fResPhiBJet = fResPhiBJet_eta1; 
  else if (fabs(eta) < 1.7) 
    fResPhiBJet = fResPhiBJet_eta2; 
  else if (fabs(eta) <= 2.5) 
    fResPhiBJet = fResPhiBJet_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorATLAS::ResPhiBJet(). Phi range exceeded." << std::endl; 
      return 0; 
    }

  return fResPhiBJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas::ResMissingET()
{
  return fResMissingET; 
}

// --------------------------------------------------------- 
