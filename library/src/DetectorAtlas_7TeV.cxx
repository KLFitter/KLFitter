#include "DetectorAtlas_7TeV.h" 
#include "ResolutionBase.h" 
#include "ResDoubleGaussPt.h" 
#include "ResDoubleGaussE_1.h"
#include "ResDoubleGaussE_2.h"
#include "ResDoubleGaussE_3.h"  
#include "ResGauss.h"
#include "ResGauss_MET.h"

#include "TString.h"

#include <cmath>
#include <iostream>
#include <stdlib.h>

// --------------------------------------------------------- 

KLFitter::DetectorAtlas_7TeV::DetectorAtlas_7TeV(std::string folder) : DetectorBase() 
{
  //check: powheg sample with 7TeV? Must use 8!
  if (strstr(folder.c_str(), "mc11c_powheg")) {
std::cout<<"ERROR! Don't use PowHeg TFs with the 7TeV Detector class!!! Exiting..."<<std::endl;
  exit(1);
  }
  //check: MC11b? New parametrization!
  if ((strstr(folder.c_str(), "mc11b"))||(strstr(folder.c_str(), "mc11c"))){
std::cout<<"Using TF from MC11b or later..."<<std::endl;
  // energy resolution 
  fResEnergyLightJet_eta1 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_lJets_eta1.txt", folder.c_str())); 
  fResEnergyBJet_eta1     = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_bJets_eta1.txt", folder.c_str())); 
  fResEnergyGluonJet_eta1 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_gluon_eta1.txt", folder.c_str())); 
  fResEnergyElectron_eta1 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_Electrons_eta1.txt", folder.c_str())); 
  fResEnergyMuon_eta1     = new KLFitter::ResDoubleGaussPt(Form("%s/par_energy_Muons_eta1.txt", folder.c_str())); 
  fResEnergyPhoton_eta1   = new KLFitter::ResGauss(Form("%s/par_energy_photon_eta1.txt", folder.c_str())); 

  fResEnergyLightJet_eta2 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_lJets_eta2.txt", folder.c_str())); 
  fResEnergyBJet_eta2     = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_bJets_eta2.txt", folder.c_str())); 
  fResEnergyGluonJet_eta2 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_gluon_eta2.txt", folder.c_str())); 
  fResEnergyElectron_eta2 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_Electrons_eta2.txt", folder.c_str())); 
  fResEnergyMuon_eta2     = new KLFitter::ResDoubleGaussPt(Form("%s/par_energy_Muons_eta2.txt", folder.c_str())); 
  fResEnergyPhoton_eta2   = new KLFitter::ResGauss(Form("%s/par_energy_photon_eta2.txt", folder.c_str())); 

  fResEnergyLightJet_eta3 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_lJets_eta3.txt", folder.c_str())); 
  fResEnergyBJet_eta3     = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_bJets_eta3.txt", folder.c_str())); 
  fResEnergyGluonJet_eta3 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_gluon_eta3.txt", folder.c_str())); 
  fResEnergyElectron_eta3 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_Electrons_eta3.txt", folder.c_str())); 
  fResEnergyMuon_eta3     = new KLFitter::ResDoubleGaussPt(Form("%s/par_energy_Muons_eta3.txt", folder.c_str())); 
  fResEnergyPhoton_eta3   = new KLFitter::ResGauss(Form("%s/par_energy_photon_eta3.txt", folder.c_str())); 

  fResEnergyLightJet_eta4 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_lJets_eta4.txt", folder.c_str())); 
  fResEnergyBJet_eta4     = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_bJets_eta4.txt", folder.c_str())); 
  fResEnergyGluonJet_eta4 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_gluon_eta4.txt", folder.c_str())); 
  fResEnergyElectron_eta4 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_Electrons_eta4.txt", folder.c_str())); 
  fResEnergyPhoton_eta4   = new KLFitter::ResGauss(Form("%s/par_energy_photon_eta4.txt", folder.c_str()));

  fResEnergyLightJet_eta5 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_lJets_eta5.txt", folder.c_str())); 
  fResEnergyBJet_eta5     = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_bJets_eta5.txt", folder.c_str()));
  } else  {
std::cout<<"Using TF from MC11a or earlier..."<<std::endl;
  // energy resolution 
  fResEnergyLightJet_eta1 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_lJets_eta1.txt", folder.c_str())); 
  fResEnergyBJet_eta1     = new KLFitter::ResDoubleGaussE_2(Form("%s/par_energy_bJets_eta1.txt", folder.c_str())); 
  fResEnergyGluonJet_eta1 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_gluon_eta1.txt", folder.c_str())); 
  fResEnergyElectron_eta1 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_Electrons_eta1.txt", folder.c_str())); 
  fResEnergyMuon_eta1     = new KLFitter::ResDoubleGaussPt(Form("%s/par_energy_Muons_eta1.txt", folder.c_str())); 
  fResEnergyPhoton_eta1   = new KLFitter::ResGauss(Form("%s/par_energy_photon_eta1.txt", folder.c_str())); 

  fResEnergyLightJet_eta2 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_lJets_eta2.txt", folder.c_str())); 
  fResEnergyBJet_eta2     = new KLFitter::ResDoubleGaussE_2(Form("%s/par_energy_bJets_eta2.txt", folder.c_str())); 
  fResEnergyGluonJet_eta2 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_gluon_eta2.txt", folder.c_str())); 
  fResEnergyElectron_eta2 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_Electrons_eta2.txt", folder.c_str())); 
  fResEnergyMuon_eta2     = new KLFitter::ResDoubleGaussPt(Form("%s/par_energy_Muons_eta2.txt", folder.c_str())); 
  fResEnergyPhoton_eta2   = new KLFitter::ResGauss(Form("%s/par_energy_photon_eta2.txt", folder.c_str())); 

  fResEnergyLightJet_eta3 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_lJets_eta3.txt", folder.c_str())); 
  fResEnergyBJet_eta3     = new KLFitter::ResDoubleGaussE_2(Form("%s/par_energy_bJets_eta3.txt", folder.c_str())); 
  fResEnergyGluonJet_eta3 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_gluon_eta3.txt", folder.c_str())); 
  fResEnergyElectron_eta3 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_Electrons_eta3.txt", folder.c_str())); 
  fResEnergyMuon_eta3     = new KLFitter::ResDoubleGaussPt(Form("%s/par_energy_Muons_eta3.txt", folder.c_str())); 
  fResEnergyPhoton_eta3   = new KLFitter::ResGauss(Form("%s/par_energy_photon_eta3.txt", folder.c_str())); 

  fResEnergyLightJet_eta4 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_lJets_eta4.txt", folder.c_str())); 
  fResEnergyBJet_eta4     = new KLFitter::ResDoubleGaussE_2(Form("%s/par_energy_bJets_eta4.txt", folder.c_str())); 
  fResEnergyGluonJet_eta4 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_gluon_eta4.txt", folder.c_str())); 
  fResEnergyElectron_eta4 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_Electrons_eta4.txt", folder.c_str())); 
  fResEnergyPhoton_eta4   = new KLFitter::ResGauss(Form("%s/par_energy_photon_eta4.txt", folder.c_str()));

  fResEnergyLightJet_eta5 = new KLFitter::ResDoubleGaussE_1(Form("%s/par_energy_lJets_eta5.txt", folder.c_str())); 
  fResEnergyBJet_eta5     = new KLFitter::ResDoubleGaussE_2(Form("%s/par_energy_bJets_eta5.txt", folder.c_str()));
  }  

  // eta resolution 
  fResEtaLightJet_eta1 = new KLFitter::ResGauss(Form("%s/par_eta_lJets_eta1.txt", folder.c_str())); 
  fResEtaLightJet_eta2 = new KLFitter::ResGauss(Form("%s/par_eta_lJets_eta2.txt", folder.c_str())); 
  fResEtaLightJet_eta3 = new KLFitter::ResGauss(Form("%s/par_eta_lJets_eta3.txt", folder.c_str()));
  fResEtaLightJet_eta4 = new KLFitter::ResGauss(Form("%s/par_eta_lJets_eta4.txt", folder.c_str()));  

  fResEtaBJet_eta1 = new KLFitter::ResGauss(Form("%s/par_eta_bJets_eta1.txt", folder.c_str())); 
  fResEtaBJet_eta2 = new KLFitter::ResGauss(Form("%s/par_eta_bJets_eta2.txt", folder.c_str())); 
  fResEtaBJet_eta3 = new KLFitter::ResGauss(Form("%s/par_eta_bJets_eta3.txt", folder.c_str()));
  fResEtaBJet_eta4 = new KLFitter::ResGauss(Form("%s/par_eta_bJets_eta4.txt", folder.c_str()));  

  // phi resolution 
  fResPhiLightJet_eta1 = new KLFitter::ResGauss(Form("%s/par_phi_lJets_eta1.txt", folder.c_str())); 
  fResPhiLightJet_eta2 = new KLFitter::ResGauss(Form("%s/par_phi_lJets_eta2.txt", folder.c_str())); 
  fResPhiLightJet_eta3 = new KLFitter::ResGauss(Form("%s/par_phi_lJets_eta3.txt", folder.c_str()));
  fResPhiLightJet_eta4 = new KLFitter::ResGauss(Form("%s/par_phi_lJets_eta4.txt", folder.c_str()));  

  fResPhiBJet_eta1 = new KLFitter::ResGauss(Form("%s/par_phi_bJets_eta1.txt", folder.c_str())); 
  fResPhiBJet_eta2 = new KLFitter::ResGauss(Form("%s/par_phi_bJets_eta2.txt", folder.c_str())); 
  fResPhiBJet_eta3 = new KLFitter::ResGauss(Form("%s/par_phi_bJets_eta3.txt", folder.c_str())); 
  fResPhiBJet_eta4 = new KLFitter::ResGauss(Form("%s/par_phi_bJets_eta4.txt", folder.c_str()));

  // missing et resolution in x and y 
  fResMissingET      = new KLFitter::ResGauss_MET(Form("%s/par_misset.txt", folder.c_str())); 

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

  // Set eta binning for different objects starting with eta=0
  fJetEtaBin_1 = 0.8;
  fJetEtaBin_2 = 1.37; 
  fJetEtaBin_3 = 1.52;
  fJetEtaBin_4 = 2.5;
  fJetEtaBin_5 = 4.5;

  fElectronEtaBin_1 = 0.8;
  fElectronEtaBin_2 = 1.37;
  fElectronEtaBin_3 = 1.52;
  fElectronEtaBin_4 = 2.5;

  fMuonEtaBin_1 = 1.11;
  fMuonEtaBin_2 = 1.25;
  fMuonEtaBin_3 = 2.5;

  fPhotonEtaBin_1 = 1.11;
  fPhotonEtaBin_2 = 1.25;
  fPhotonEtaBin_3 = 2.5;
  fPhotonEtaBin_4 = 3.0;
    
}

// --------------------------------------------------------- 

KLFitter::DetectorAtlas_7TeV::~DetectorAtlas_7TeV() 
{
  if (fResEnergyLightJet_eta1) 
    delete fResEnergyLightJet_eta1; 
        
  if (fResEnergyLightJet_eta2)
    delete fResEnergyLightJet_eta2; 

  if (fResEnergyLightJet_eta3)
    delete fResEnergyLightJet_eta3;

  if (fResEnergyLightJet_eta4)
    delete fResEnergyLightJet_eta4;

  if (fResEnergyLightJet_eta5)
    delete fResEnergyLightJet_eta5;  

  if (fResEnergyBJet_eta1)
    delete fResEnergyBJet_eta1; 

  if (fResEnergyBJet_eta2)
    delete fResEnergyBJet_eta2; 

  if (fResEnergyBJet_eta3)
    delete fResEnergyBJet_eta3;

  if (fResEnergyBJet_eta4)
    delete fResEnergyBJet_eta4;

  if (fResEnergyBJet_eta5)
    delete fResEnergyBJet_eta5;    

  if (fResEnergyGluonJet_eta1)
    delete fResEnergyGluonJet_eta1; 

  if (fResEnergyGluonJet_eta2)
    delete fResEnergyGluonJet_eta2; 

  if (fResEnergyGluonJet_eta3)
    delete fResEnergyGluonJet_eta3;

 if (fResEnergyGluonJet_eta4)
    delete fResEnergyGluonJet_eta4; 

  if (fResEnergyElectron_eta1)
    delete fResEnergyElectron_eta1; 

  if (fResEnergyElectron_eta2)
    delete fResEnergyElectron_eta2; 

  if (fResEnergyElectron_eta3)
    delete fResEnergyElectron_eta3;

  if (fResEnergyElectron_eta4)
    delete fResEnergyElectron_eta4; 

  if (fResEnergyPhoton_eta1)
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

  if (fResEnergyPhoton_eta4)
    delete fResEnergyPhoton_eta4; 

  if (fResEtaLightJet_eta1) 
    delete fResEtaLightJet_eta1; 

  if (fResEtaLightJet_eta2) 
    delete fResEtaLightJet_eta2; 

  if (fResEtaLightJet_eta3) 
    delete fResEtaLightJet_eta3;

  if (fResEtaLightJet_eta4) 
    delete fResEtaLightJet_eta4; 
        
  if (fResEtaBJet_eta1) 
    delete fResEtaBJet_eta1; 

  if (fResEtaBJet_eta2) 
    delete fResEtaBJet_eta2; 

  if (fResEtaBJet_eta3) 
    delete fResEtaBJet_eta3;

  if (fResEtaBJet_eta4) 
    delete fResEtaBJet_eta4; 
        
  if (fResPhiLightJet_eta1) 
    delete fResPhiLightJet_eta1; 

  if (fResPhiLightJet_eta2) 
    delete fResPhiLightJet_eta2; 

  if (fResPhiLightJet_eta3) 
    delete fResPhiLightJet_eta3;

  if (fResPhiLightJet_eta4) 
    delete fResPhiLightJet_eta4;  

  if (fResPhiBJet_eta1) 
    delete fResPhiBJet_eta1; 

  if (fResPhiBJet_eta2) 
    delete fResPhiBJet_eta2; 

  if (fResPhiBJet_eta3) 
    delete fResPhiBJet_eta3;

  if (fResPhiBJet_eta4) 
    delete fResPhiBJet_eta4;  
        
  if (fResMissingET)
    delete fResMissingET; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_7TeV::ResEnergyLightJet(double eta)
{
  if (fabs(eta) < fJetEtaBin_1) 
    fResEnergyLightJet = fResEnergyLightJet_eta1; 
  else if (fabs(eta) < fJetEtaBin_2) 
    fResEnergyLightJet = fResEnergyLightJet_eta2; 
  else if (fabs(eta) < fJetEtaBin_3) 
    fResEnergyLightJet = fResEnergyLightJet_eta3;
  else if (fabs(eta) <= fJetEtaBin_4) 
    fResEnergyLightJet = fResEnergyLightJet_eta4;
  else if (fabs(eta) <= fJetEtaBin_5) 
    fResEnergyLightJet = fResEnergyLightJet_eta5;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_7TeV::ResEnergyLightJet(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEnergyLightJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_7TeV::ResEnergyBJet(double eta)
{
  if (fabs(eta) < fJetEtaBin_1) 
    fResEnergyBJet = fResEnergyBJet_eta1; 
  else if (fabs(eta) < fJetEtaBin_2) 
    fResEnergyBJet = fResEnergyBJet_eta2; 
  else if (fabs(eta) < fJetEtaBin_3) 
    fResEnergyBJet = fResEnergyBJet_eta3;
  else if (fabs(eta) <= fJetEtaBin_4) 
    fResEnergyBJet = fResEnergyBJet_eta4;
  else if (fabs(eta) <= fJetEtaBin_5) 
    fResEnergyBJet = fResEnergyBJet_eta5;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_7TeV::ResEnergyBJet(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEnergyBJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_7TeV::ResEnergyGluonJet(double eta)
{
  if (fabs(eta) < fJetEtaBin_1) 
    fResEnergyGluonJet = fResEnergyGluonJet_eta1; 
  else if (fabs(eta) < fJetEtaBin_2) 
    fResEnergyGluonJet = fResEnergyGluonJet_eta2; 
  else if (fabs(eta) < fJetEtaBin_3) 
    fResEnergyGluonJet = fResEnergyGluonJet_eta3;
  else if (fabs(eta) <= fJetEtaBin_4) 
    fResEnergyGluonJet = fResEnergyGluonJet_eta4;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_7TeV::ResEnergyGluonJet(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEnergyGluonJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_7TeV::ResEnergyElectron(double eta)
{
  if (fabs(eta) < fElectronEtaBin_1) 
    fResEnergyElectron = fResEnergyElectron_eta1; 
  else if (fabs(eta) < fElectronEtaBin_2) 
    fResEnergyElectron = fResEnergyElectron_eta2; 
  else if (fabs(eta) < fElectronEtaBin_3){ 
    std::cout << "KLFitter::DetectorAtlas_7TeV::ResEnergyElectron(). Electron in crack region" << std::endl;
    return 0;
  }
  else if (fabs(eta) <= fElectronEtaBin_4) 
    fResEnergyElectron = fResEnergyElectron_eta4;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_7TeV::ResEnergyElectron(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEnergyElectron; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_7TeV::ResEnergyMuon(double eta)
{
  if (fabs(eta) < fMuonEtaBin_1) 
    fResEnergyMuon = fResEnergyMuon_eta1; 
  else if (fabs(eta) < fMuonEtaBin_2) 
    fResEnergyMuon = fResEnergyMuon_eta2; 
  else if (fabs(eta) < fMuonEtaBin_3) 
    fResEnergyMuon = fResEnergyMuon_eta3;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_7TeV::ResEnergyMuon(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEnergyMuon; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_7TeV::ResEnergyPhoton(double eta)
{
  if (fabs(eta) < fPhotonEtaBin_1) 
    fResEnergyPhoton = fResEnergyPhoton_eta1; 
  else if (fabs(eta) < fPhotonEtaBin_2) 
    fResEnergyPhoton = fResEnergyPhoton_eta2; 
  else if (fabs(eta) < fPhotonEtaBin_3) 
    fResEnergyPhoton = fResEnergyPhoton_eta3;
  else if (fabs(eta) <= fPhotonEtaBin_4) 
    fResEnergyPhoton = fResEnergyPhoton_eta4;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_7TeV::ResEnergyPhoton(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEnergyPhoton; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_7TeV::ResEtaLightJet(double eta)
{
  if (fabs(eta) < fJetEtaBin_1) 
    fResEtaLightJet = fResEtaLightJet_eta1; 
  else if (fabs(eta) < fJetEtaBin_2) 
    fResEtaLightJet = fResEtaLightJet_eta2; 
  else if (fabs(eta) < fJetEtaBin_3) 
    fResEtaLightJet = fResEtaLightJet_eta3;
  else if (fabs(eta) <= fJetEtaBin_4) 
    fResEtaLightJet = fResEtaLightJet_eta4;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_7TeV::ResEtaLightJet(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEtaLightJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_7TeV::ResEtaBJet(double eta)
{
  if (fabs(eta) < fJetEtaBin_1) 
    fResEtaBJet = fResEtaBJet_eta1; 
  else if (fabs(eta) < fJetEtaBin_2) 
    fResEtaBJet = fResEtaBJet_eta2; 
  else if (fabs(eta) < fJetEtaBin_3) 
    fResEtaBJet = fResEtaBJet_eta3;
  else if (fabs(eta) <= fJetEtaBin_4) 
    fResEtaBJet = fResEtaBJet_eta4;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_7TeV::ResEtaBJet(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResEtaBJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_7TeV::ResPhiLightJet(double eta)
{
  if (fabs(eta) < fJetEtaBin_1) 
    fResPhiLightJet = fResPhiLightJet_eta1; 
  else if (fabs(eta) < fJetEtaBin_2) 
    fResPhiLightJet = fResPhiLightJet_eta2; 
  else if (fabs(eta) < fJetEtaBin_3) 
    fResPhiLightJet = fResPhiLightJet_eta3;
  else if (fabs(eta) <= fJetEtaBin_4) 
    fResPhiLightJet = fResPhiLightJet_eta4;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_7TeV::ResPhiLightJet(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResPhiLightJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_7TeV::ResPhiBJet(double eta)
{
  if (fabs(eta) < fJetEtaBin_1) 
    fResPhiBJet = fResPhiBJet_eta1; 
  else if (fabs(eta) < fJetEtaBin_2) 
    fResPhiBJet = fResPhiBJet_eta2; 
  else if (fabs(eta) < fJetEtaBin_3) 
    fResPhiBJet = fResPhiBJet_eta3;
  else if (fabs(eta) <= fJetEtaBin_4) 
    fResPhiLightJet = fResPhiLightJet_eta4;
  else 
    {
      std::cout << "KLFitter::DetectorAtlas_7TeV::ResPhiBJet(). Eta range exceeded." << std::endl; 
      return 0; 
    }

  return fResPhiBJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::DetectorAtlas_7TeV::ResMissingET()
{
  return fResMissingET; 
}

// --------------------------------------------------------- 
