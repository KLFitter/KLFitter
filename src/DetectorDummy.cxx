#include "DetectorDummy.h" 
#include "ResolutionBase.h"
#include "ResGauss.h"

// --------------------------------------------------------- 

KLFitter::DetectorDummy::DetectorDummy(std::string folder) : DetectorBase() 
{
  fResEnergyLightJet = new KLFitter::ResGauss(1.0); 
  fResEnergyBJet     = new KLFitter::ResGauss(1.0); 
  fResEnergyGluonJet = new KLFitter::ResGauss(1.0); 
  fResEnergyElectron = new KLFitter::ResGauss(0.1);
  fResEnergyMuon     = new KLFitter::ResGauss(0.1); 
  fResEnergyPhoton   = new KLFitter::ResGauss(0.1); 
  fResMissingET      = new KLFitter::ResGauss(1.0); 
  fResEtaLightJet    = new KLFitter::ResGauss(0.01); 
  fResEtaBJet        = new KLFitter::ResGauss(0.01); 
  fResPhiLightJet    = new KLFitter::ResGauss(0.01); 
  fResPhiBJet        = new KLFitter::ResGauss(0.01); 
}

// --------------------------------------------------------- 

KLFitter::DetectorDummy::~DetectorDummy() 
{
}

// --------------------------------------------------------- 

