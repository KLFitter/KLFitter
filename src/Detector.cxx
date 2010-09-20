#include "Detector.h" 
#include "ResolutionBase.h" 
#include "ResGauss.h"
#include "ResDoubleGaussE.h"
#include "ResDoubleGaussPt.h"

#include <cmath>
#include <iostream>

// --------------------------------------------------------- 

KLFitter::Detector::Detector(std::vector<double> MaxValuesForEtaRegions,
                             std::string eResTypeLightJet,
                             std::string eResTypeBJet,
                             std::string eResTypeGluonJet,
                             std::string eResTypeElectron,
                             std::string eResTypeMuon,
                             std::string eResTypePhoton,
                             std::vector<double> eResParLightJet,
                             std::vector<double> eResParBJet,
                             std::vector<double> eResParGluonJet,
                             std::vector<double> eResParElectron,
                             std::vector<double> eResParMuon,
                             std::vector<double> eResParPhoton) : DetectorBase() {

  fMaxValuesForEtaRegions = MaxValuesForEtaRegions;
  fNEtaRegions = fMaxValuesForEtaRegions.size();

  fResEnergyLightJet_etaRegions = new std::vector<KLFitter::ResolutionBase*>;
  fResEnergyBJet_etaRegions     = new std::vector<KLFitter::ResolutionBase*>;
  fResEnergyGluonJet_etaRegions = new std::vector<KLFitter::ResolutionBase*>;
  fResEnergyElectron_etaRegions = new std::vector<KLFitter::ResolutionBase*>;
  fResEnergyMuon_etaRegions     = new std::vector<KLFitter::ResolutionBase*>;
  fResEnergyPhoton_etaRegions     = new std::vector<KLFitter::ResolutionBase*>;

  for (unsigned int iEtaRegion = 0; iEtaRegion < fNEtaRegions; iEtaRegion++) {
    if (eResTypeLightJet=="Gauss")
      fResEnergyLightJet_etaRegions->push_back(new KLFitter::ResGauss(eResParLightJet.at(iEtaRegion)));
    else if (eResTypeLightJet=="DoubleGauss")
      fResEnergyLightJet_etaRegions->push_back(new KLFitter::ResDoubleGaussE(getDoubleGaussParameters(eResParLightJet, iEtaRegion)));
    else
      std::cout << "energy resolution type for light jets not recognised" << std::endl;

    if (eResTypeBJet=="Gauss")
      fResEnergyBJet_etaRegions->push_back(new KLFitter::ResGauss(eResParBJet.at(0)));
    else if (eResTypeBJet=="DoubleGauss")
      fResEnergyBJet_etaRegions->push_back(new KLFitter::ResDoubleGaussE(getDoubleGaussParameters(eResParBJet, iEtaRegion)));
    else
      std::cout << "energy resolution type for b jets not recognised" << std::endl;

    if (eResTypeGluonJet=="Gauss")
      fResEnergyGluonJet_etaRegions->push_back(new KLFitter::ResGauss(eResParGluonJet.at(0)));
    else if (eResTypeGluonJet=="DoubleGauss")
      fResEnergyGluonJet_etaRegions->push_back(new KLFitter::ResDoubleGaussE(getDoubleGaussParameters(eResParGluonJet, iEtaRegion)));
    else
      std::cout << "energy resolution type for gluon jets not recognised" << std::endl;

    if (eResTypeElectron=="Gauss")
      fResEnergyElectron_etaRegions->push_back(new KLFitter::ResGauss(eResParElectron.at(0)));
    else if (eResTypeElectron=="DoubleGauss")
      fResEnergyElectron_etaRegions->push_back(new KLFitter::ResDoubleGaussE(getDoubleGaussParameters(eResParElectron, iEtaRegion)));
    else
      std::cout << "energy resolution type for electrons not recognised" << std::endl;

    if (eResTypeMuon=="Gauss")
      fResEnergyMuon_etaRegions->push_back(new KLFitter::ResGauss(eResParMuon.at(0)));
    else if (eResTypeMuon=="DoubleGauss")
      fResEnergyMuon_etaRegions->push_back(new KLFitter::ResDoubleGaussPt(getDoubleGaussParameters(eResParMuon, iEtaRegion)));
    else
      std::cout << "energy resolution type for muons not recognised" << std::endl;

    if (eResTypePhoton=="Gauss")
      fResEnergyPhoton_etaRegions->push_back(new KLFitter::ResGauss(eResParPhoton.at(0)));
    else if (eResTypePhoton=="DoubleGauss")
      fResEnergyPhoton_etaRegions->push_back(new KLFitter::ResDoubleGaussE(getDoubleGaussParameters(eResParPhoton, iEtaRegion)));
    else
      std::cout << "energy resolution type for photons not recognised" << std::endl;

  }

  // just to make DetectorBase::Status() happy
  fResEnergyLightJet = new KLFitter::ResGauss(1.0);
  fResEnergyBJet     = new KLFitter::ResGauss(1.0);
  fResEnergyGluonJet = new KLFitter::ResGauss(1.0);
  fResEnergyElectron = new KLFitter::ResGauss(1.0);
  fResEnergyMuon     = new KLFitter::ResGauss(1.0);
  fResEnergyPhoton     = new KLFitter::ResGauss(1.0);

  fResMissingET = new KLFitter::ResGauss(1.0); 
}

// --------------------------------------------------------- 

KLFitter::Detector::Detector(std::vector<double> MaxValuesForEtaRegions,
                             std::string eResTypeLightJet,
                             std::string eResTypeBJet,
                             std::string eResTypeGluonJet,
                             std::string eResTypeElectron,
                             std::string eResTypeMuon,
                             std::string eResTypePhoton,
                             std::string etaResTypeLightJet,
                             std::string etaResTypeBJet,
                             std::string phiResTypeLightJet,
                             std::string phiResTypeBJet,
                             std::string metResType,
                             std::vector<double> eResParLightJet,
                             std::vector<double> eResParBJet,
                             std::vector<double> eResParGluonJet,
                             std::vector<double> eResParElectron,
                             std::vector<double> eResParMuon,
                             std::vector<double> eResParPhoton,
                             std::vector<double> etaResParLightJet,
                             std::vector<double> etaResParBJet,
                             std::vector<double> phiResParLightJet,
                             std::vector<double> phiResParBJet,
                             std::vector<double> metResPar) : DetectorBase() {

  fMaxValuesForEtaRegions = MaxValuesForEtaRegions;
  fNEtaRegions = fMaxValuesForEtaRegions.size();

  fResEnergyLightJet_etaRegions = new std::vector<KLFitter::ResolutionBase*>;
  fResEnergyBJet_etaRegions     = new std::vector<KLFitter::ResolutionBase*>;
  fResEnergyGluonJet_etaRegions = new std::vector<KLFitter::ResolutionBase*>;
  fResEnergyElectron_etaRegions = new std::vector<KLFitter::ResolutionBase*>;
  fResEnergyMuon_etaRegions     = new std::vector<KLFitter::ResolutionBase*>;
  fResEnergyPhoton_etaRegions   = new std::vector<KLFitter::ResolutionBase*>;
  fResEtaLightJet_etaRegions    = new std::vector<KLFitter::ResolutionBase*>;
  fResEtaBJet_etaRegions        = new std::vector<KLFitter::ResolutionBase*>;
  fResPhiLightJet_etaRegions    = new std::vector<KLFitter::ResolutionBase*>;
  fResPhiBJet_etaRegions        = new std::vector<KLFitter::ResolutionBase*>;
  fResMissingET_etaRegions      = new std::vector<KLFitter::ResolutionBase*>;

  for (unsigned int iEtaRegion = 0; iEtaRegion < fNEtaRegions; iEtaRegion++) {
    if (eResTypeLightJet=="Gauss")
      fResEnergyLightJet_etaRegions->push_back(new KLFitter::ResGauss(eResParLightJet.at(iEtaRegion)));
    else if (eResTypeLightJet=="DoubleGauss")
      fResEnergyLightJet_etaRegions->push_back(new KLFitter::ResDoubleGaussE(getDoubleGaussParameters(eResParLightJet, iEtaRegion)));
    else
      std::cout << "energy resolution type for light jets not recognised" << std::endl;

    if (eResTypeBJet=="Gauss")
      fResEnergyBJet_etaRegions->push_back(new KLFitter::ResGauss(eResParBJet.at(0)));
    else if (eResTypeBJet=="DoubleGauss")
      fResEnergyBJet_etaRegions->push_back(new KLFitter::ResDoubleGaussE(getDoubleGaussParameters(eResParBJet, iEtaRegion)));
    else
      std::cout << "energy resolution type for b jets not recognised" << std::endl;

    if (eResTypeGluonJet=="Gauss")
      fResEnergyGluonJet_etaRegions->push_back(new KLFitter::ResGauss(eResParGluonJet.at(0)));
    else if (eResTypeGluonJet=="DoubleGauss")
      fResEnergyGluonJet_etaRegions->push_back(new KLFitter::ResDoubleGaussE(getDoubleGaussParameters(eResParGluonJet, iEtaRegion)));
    else
      std::cout << "energy resolution type for gluon jets not recognised" << std::endl;

    if (eResTypeElectron=="Gauss")
      fResEnergyElectron_etaRegions->push_back(new KLFitter::ResGauss(eResParElectron.at(0)));
    else if (eResTypeElectron=="DoubleGauss")
      fResEnergyElectron_etaRegions->push_back(new KLFitter::ResDoubleGaussE(getDoubleGaussParameters(eResParElectron, iEtaRegion)));
    else
      std::cout << "energy resolution type for electrons not recognised" << std::endl;

    if (eResTypeMuon=="Gauss")
      fResEnergyMuon_etaRegions->push_back(new KLFitter::ResGauss(eResParMuon.at(0)));
    else if (eResTypeMuon=="DoubleGauss")
      fResEnergyMuon_etaRegions->push_back(new KLFitter::ResDoubleGaussPt(getDoubleGaussParameters(eResParMuon, iEtaRegion)));
    else
      std::cout << "energy resolution type for muons not recognised" << std::endl;

    if (eResTypePhoton=="Gauss")
      fResEnergyPhoton_etaRegions->push_back(new KLFitter::ResGauss(eResParPhoton.at(0)));
    else if (eResTypePhoton=="DoubleGauss")
      fResEnergyPhoton_etaRegions->push_back(new KLFitter::ResDoubleGaussE(getDoubleGaussParameters(eResParPhoton, iEtaRegion)));
    else
      std::cout << "energy resolution type for photons not recognised" << std::endl;

    if (etaResTypeLightJet=="Gauss")
      fResEtaLightJet_etaRegions->push_back(new KLFitter::ResGauss(etaResParLightJet.at(iEtaRegion)));
    else if (etaResTypeLightJet=="DoubleGauss")
      fResEtaLightJet_etaRegions->push_back(new KLFitter::ResDoubleGaussE(getDoubleGaussParameters(etaResParLightJet, iEtaRegion)));
    else
      std::cout << "eta resolution type for light jets not recognised" << std::endl;

    if (etaResTypeBJet=="Gauss")
      fResEtaBJet_etaRegions->push_back(new KLFitter::ResGauss(etaResParBJet.at(0)));
    else if (etaResTypeBJet=="DoubleGauss")
      fResEtaBJet_etaRegions->push_back(new KLFitter::ResDoubleGaussE(getDoubleGaussParameters(etaResParBJet, iEtaRegion)));
    else
      std::cout << "eta resolution type for b jets not recognised" << std::endl;

    if (phiResTypeLightJet=="Gauss")
      fResPhiLightJet_etaRegions->push_back(new KLFitter::ResGauss(phiResParLightJet.at(iEtaRegion)));
    else if (phiResTypeLightJet=="DoubleGauss")
      fResPhiLightJet_etaRegions->push_back(new KLFitter::ResDoubleGaussE(getDoubleGaussParameters(phiResParLightJet, iEtaRegion)));
    else
      std::cout << "phi resolution type for light jets not recognised" << std::endl;

    if (phiResTypeBJet=="Gauss")
      fResPhiBJet_etaRegions->push_back(new KLFitter::ResGauss(phiResParBJet.at(0)));
    else if (phiResTypeBJet=="DoubleGauss")
      fResPhiBJet_etaRegions->push_back(new KLFitter::ResDoubleGaussE(getDoubleGaussParameters(phiResParBJet, iEtaRegion)));
    else
      std::cout << "phi resolution type for b jets not recognised" << std::endl;

    if (metResType=="Gauss")
      fResMissingET_etaRegions->push_back(new KLFitter::ResGauss(metResPar.at(0)));
    else if (metResType=="DoubleGauss")
      fResMissingET_etaRegions->push_back(new KLFitter::ResDoubleGaussE(getDoubleGaussParameters(metResPar, iEtaRegion)));
    else
      std::cout << "resolution type for missing ET not recognised" << std::endl;
  }

  // just to make DetectorBase::Status() happy
  fResEnergyLightJet = new KLFitter::ResGauss(1.0);
  fResEnergyBJet     = new KLFitter::ResGauss(1.0);
  fResEnergyGluonJet = new KLFitter::ResGauss(1.0);
  fResEnergyElectron = new KLFitter::ResGauss(1.0);
  fResEnergyMuon     = new KLFitter::ResGauss(1.0);
  fResEnergyPhoton   = new KLFitter::ResGauss(1.0);
  fResEtaLightJet    = new KLFitter::ResGauss(1.0);
  fResEtaBJet        = new KLFitter::ResGauss(1.0);
  fResPhiLightJet    = new KLFitter::ResGauss(1.0);
  fResPhiBJet        = new KLFitter::ResGauss(1.0);
  fResMissingET = new KLFitter::ResGauss(1.0); 
}

// --------------------------------------------------------- 

KLFitter::Detector::~Detector() {

  if (fResEnergyLightJet_etaRegions) { 
    for (unsigned int i = 0; i < fResEnergyLightJet_etaRegions->size(); i++)
      delete fResEnergyLightJet_etaRegions->at(i);
    delete fResEnergyLightJet_etaRegions; 
  }

  if (fResEnergyBJet_etaRegions) { 
    for (unsigned int i = 0; i < fResEnergyBJet_etaRegions->size(); i++)
      delete fResEnergyBJet_etaRegions->at(i);
    delete fResEnergyBJet_etaRegions; 
  }

  if (fResEnergyGluonJet_etaRegions) { 
    for (unsigned int i = 0; i < fResEnergyGluonJet_etaRegions->size(); i++)
      delete fResEnergyGluonJet_etaRegions->at(i);
    delete fResEnergyGluonJet_etaRegions; 
  }

  if (fResEnergyElectron_etaRegions) { 
    for (unsigned int i = 0; i < fResEnergyElectron_etaRegions->size(); i++)
      delete fResEnergyElectron_etaRegions->at(i);
    delete fResEnergyElectron_etaRegions; 
  }

  if (fResEnergyMuon_etaRegions) { 
    for (unsigned int i = 0; i < fResEnergyMuon_etaRegions->size(); i++)
      delete fResEnergyMuon_etaRegions->at(i);
    delete fResEnergyMuon_etaRegions; 
  }

  if (fResEnergyPhoton_etaRegions) { 
    for (unsigned int i = 0; i < fResEnergyPhoton_etaRegions->size(); i++)
      delete fResEnergyPhoton_etaRegions->at(i);
    delete fResEnergyPhoton_etaRegions; 
  }

  if (fResEtaLightJet_etaRegions) { 
    for (unsigned int i = 0; i < fResEtaLightJet_etaRegions->size(); i++)
      delete fResEtaLightJet_etaRegions->at(i);
    delete fResEtaLightJet_etaRegions; 
  }

  if (fResEtaBJet_etaRegions) { 
    for (unsigned int i = 0; i < fResEtaBJet_etaRegions->size(); i++)
      delete fResEtaBJet_etaRegions->at(i);
    delete fResEtaBJet_etaRegions; 
  }

  if (fResPhiLightJet_etaRegions) { 
    for (unsigned int i = 0; i < fResPhiLightJet_etaRegions->size(); i++)
      delete fResPhiLightJet_etaRegions->at(i);
    delete fResPhiLightJet_etaRegions; 
  }

  if (fResPhiBJet_etaRegions) { 
    for (unsigned int i = 0; i < fResPhiBJet_etaRegions->size(); i++)
      delete fResPhiBJet_etaRegions->at(i);
    delete fResPhiBJet_etaRegions; 
  }

  if (fResMissingET_etaRegions) {
    for (unsigned int i = 0; i < fResMissingET_etaRegions->size(); i++)
      delete fResMissingET_etaRegions->at(i);
    delete fResMissingET_etaRegions; 
  }

}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::Detector::ResEnergyLightJet(double eta) {
  bool foundEtaRegion = false;
  for (unsigned int iEtaRegion = 0; iEtaRegion < fNEtaRegions; iEtaRegion++) {
    if (fabs(eta) <= fMaxValuesForEtaRegions.at(iEtaRegion) && !foundEtaRegion) {
      fResEnergyLightJet = fResEnergyLightJet_etaRegions->at(iEtaRegion);
      foundEtaRegion = true;
    }
  }
  if (!foundEtaRegion) {
    std::cout << "KLFitter::Detector::ResEnergyLightJet(). Eta range exceeded." << std::endl; 
    return 0; 
  }

  return fResEnergyLightJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::Detector::ResEnergyBJet(double eta) {
  bool foundEtaRegion = false;
  for (unsigned int iEtaRegion = 0; iEtaRegion < fNEtaRegions; iEtaRegion++) {
    if (fabs(eta) <= fMaxValuesForEtaRegions.at(iEtaRegion) && !foundEtaRegion) {
      fResEnergyBJet = fResEnergyBJet_etaRegions->at(iEtaRegion);
      foundEtaRegion = true;
    }
  }
  if (!foundEtaRegion) {
    std::cout << "KLFitter::Detector::ResEnergyBJet(). Eta range exceeded." << std::endl; 
    return 0; 
  }

  return fResEnergyBJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::Detector::ResEnergyGluonJet(double eta) {
  bool foundEtaRegion = false;
  for (unsigned int iEtaRegion = 0; iEtaRegion < fNEtaRegions; iEtaRegion++) {
    if (fabs(eta) <= fMaxValuesForEtaRegions.at(iEtaRegion) && !foundEtaRegion) {
      fResEnergyGluonJet = fResEnergyGluonJet_etaRegions->at(iEtaRegion);
      foundEtaRegion = true;
    }
  }
  if (!foundEtaRegion) {
    std::cout << "KLFitter::Detector::ResEnergyGluonJet(). Eta range exceeded." << std::endl; 
    return 0; 
  }

  return fResEnergyGluonJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::Detector::ResEnergyElectron(double eta) {
  bool foundEtaRegion = false;
  for (unsigned int iEtaRegion = 0; iEtaRegion < fNEtaRegions; iEtaRegion++) {
    if (fabs(eta) <= fMaxValuesForEtaRegions.at(iEtaRegion) && !foundEtaRegion) {
      fResEnergyElectron = fResEnergyElectron_etaRegions->at(iEtaRegion);
      foundEtaRegion = true;
    }
  }
  if (!foundEtaRegion) {
    std::cout << "KLFitter::Detector::ResEnergyElectron(). Eta range exceeded." << std::endl; 
    return 0; 
  }

  return fResEnergyElectron; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::Detector::ResEnergyMuon(double eta) {
  bool foundEtaRegion = false;
  for (unsigned int iEtaRegion = 0; iEtaRegion < fNEtaRegions; iEtaRegion++) {
    if (fabs(eta) <= fMaxValuesForEtaRegions.at(iEtaRegion) && !foundEtaRegion) {
      fResEnergyMuon = fResEnergyMuon_etaRegions->at(iEtaRegion);
      foundEtaRegion = true;
    }
  }
  if (!foundEtaRegion) {
    std::cout << "KLFitter::Detector::ResEnergyMuon(). Eta range exceeded." << std::endl; 
    return 0; 
  }

  return fResEnergyMuon; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::Detector::ResEnergyPhoton(double eta) {
  bool foundEtaRegion = false;
  for (unsigned int iEtaRegion = 0; iEtaRegion < fNEtaRegions; iEtaRegion++) {
    if (fabs(eta) <= fMaxValuesForEtaRegions.at(iEtaRegion) && !foundEtaRegion) {
      fResEnergyPhoton = fResEnergyPhoton_etaRegions->at(iEtaRegion);
      foundEtaRegion = true;
    }
  }
  if (!foundEtaRegion) {
    std::cout << "KLFitter::Detector::ResEnergyPhoton(). Eta range exceeded." << std::endl; 
    return 0; 
  }

  return fResEnergyPhoton; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::Detector::ResEtaLightJet(double eta) {
  bool foundEtaRegion = false;
  for (unsigned int iEtaRegion = 0; iEtaRegion < fNEtaRegions; iEtaRegion++) {
    if (fabs(eta) <= fMaxValuesForEtaRegions.at(iEtaRegion) && !foundEtaRegion) {
      fResEtaLightJet = fResEtaLightJet_etaRegions->at(iEtaRegion);
      foundEtaRegion = true;
    }
  }
  if (!foundEtaRegion) {
    std::cout << "KLFitter::Detector::ResEtaLightJet(). Eta range exceeded." << std::endl; 
    return 0; 
  }

  return fResEtaLightJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::Detector::ResEtaBJet(double eta) {
  bool foundEtaRegion = false;
  for (unsigned int iEtaRegion = 0; iEtaRegion < fNEtaRegions; iEtaRegion++) {
    if (fabs(eta) <= fMaxValuesForEtaRegions.at(iEtaRegion) && !foundEtaRegion) {
      fResEtaBJet = fResEtaBJet_etaRegions->at(iEtaRegion);
      foundEtaRegion = true;
    }
  }
  if (!foundEtaRegion) {
    std::cout << "KLFitter::Detector::ResEtaBJet(). Eta range exceeded." << std::endl; 
    return 0; 
  }

  return fResEtaBJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::Detector::ResPhiLightJet(double eta) {
  bool foundEtaRegion = false;
  for (unsigned int iEtaRegion = 0; iEtaRegion < fNEtaRegions; iEtaRegion++) {
    if (fabs(eta) <= fMaxValuesForEtaRegions.at(iEtaRegion) && !foundEtaRegion) {
      fResPhiLightJet = fResPhiLightJet_etaRegions->at(iEtaRegion);
      foundEtaRegion = true;
    }
  }
  if (!foundEtaRegion) {
    std::cout << "KLFitter::Detector::ResPhiLightJet(). Eta range exceeded." << std::endl; 
    return 0; 
  }

  return fResPhiLightJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::Detector::ResPhiBJet(double eta) {
  bool foundEtaRegion = false;
  for (unsigned int iEtaRegion = 0; iEtaRegion < fNEtaRegions; iEtaRegion++) {
    if (fabs(eta) <= fMaxValuesForEtaRegions.at(iEtaRegion) && !foundEtaRegion) {
      fResPhiBJet = fResPhiBJet_etaRegions->at(iEtaRegion);
      foundEtaRegion = true;
    }
  }
  if (!foundEtaRegion) {
    std::cout << "KLFitter::Detector::ResPhiBJet(). Eta range exceeded." << std::endl; 
    return 0; 
  }

  return fResPhiBJet; 
}

// --------------------------------------------------------- 
KLFitter::ResolutionBase * KLFitter::Detector::ResMissingET() {
  fResMissingET = fResMissingET_etaRegions->at(0);

  return fResMissingET; 
}

// --------------------------------------------------------- 
std::vector<double> KLFitter::Detector::getDoubleGaussParameters(std::vector<double> parameters, int iEtaRegion) {
  std::vector<double> currentParameters;

  unsigned indexFirst = iEtaRegion*10;
  unsigned indexLast  = (iEtaRegion+1)*10;
  for (unsigned int index = indexFirst; index < indexLast; index++)
    currentParameters.push_back(parameters.at(index));

  return currentParameters;
}
