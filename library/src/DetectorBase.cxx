#include "DetectorBase.h" 
#include "ResolutionBase.h"
#include <iostream> 

// --------------------------------------------------------- 
KLFitter::DetectorBase::DetectorBase(std::string KLFITTER_UNUSED(folder)) :
  fResEnergyLightJet(0), 
  fResEnergyBJet(0), 
  fResEnergyGluonJet(0), 
  fResEnergyElectron(0), 
  fResEnergyMuon(0), 
  fResEnergyPhoton(0), 
  fResMissingET(0)
{
}

// --------------------------------------------------------- 
KLFitter::DetectorBase::~DetectorBase()
{
}

// --------------------------------------------------------- 
int KLFitter::DetectorBase::SetResEnergyBJet(KLFitter::ResolutionBase * res)
{
  // set resolution 
  fResEnergyBJet = res; 
        
  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::DetectorBase::SetResEnergyLightJet(KLFitter::ResolutionBase * res)
{
  // set resolution 
  fResEnergyLightJet = res; 
        
  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::DetectorBase::SetResEnergyGluonJet(KLFitter::ResolutionBase * res)
{
  // set resolution 
  fResEnergyGluonJet = res; 
        
  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::DetectorBase::SetResEnergyElectron(KLFitter::ResolutionBase * res)
{
  // set resolution 
  fResEnergyElectron = res; 
        
  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::DetectorBase::SetResEnergyMuon(KLFitter::ResolutionBase * res)
{
  // set resolution 
  fResEnergyMuon = res; 
        
  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::DetectorBase::SetResEnergyPhoton(KLFitter::ResolutionBase * res)
{
  // set resolution 
  fResEnergyPhoton = res; 
        
  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::DetectorBase::SetResMissingET(KLFitter::ResolutionBase * res)
{
  // set resolution 
  fResMissingET = res; 
        
  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::DetectorBase::Status()
{
  if (!fResEnergyLightJet)
    {
      std::cout << "KLFitter::DetectorBase::Status(). Energy resolution of light jets not defined." << std::endl;
      return 0; 
    }

  if (!fResEnergyBJet)
    {
      std::cout << "KLFitter::DetectorBase::Status(). Energy resolution of b jets not defined." << std::endl;
      return 0; 
    }

  if (!fResEnergyGluonJet)
    {
      std::cout << "KLFitter::DetectorBase::Status(). Energy resolution of gluon jets not defined." << std::endl;
      return 0; 
    }

  if (!fResEnergyElectron)
    {
      std::cout << "KLFitter::DetectorBase::Status(). Energy resolution of electrons not defined." << std::endl;
      return 0; 
    }

  if (!fResEnergyMuon)
    {
      std::cout << "KLFitter::DetectorBase::Status(). Energy resolution of muons not defined." << std::endl;
      return 0; 
    }

  if (!fResEnergyPhoton)
    {
      std::cout << "KLFitter::DetectorBase::Status(). Energy resolution of photons not defined." << std::endl;
      return 0; 
    }

  if (!fResMissingET)
    {
      std::cout << "KLFitter::DetectorBase::Status(). Missing ET resolution not defined." << std::endl;
      return 0; 
    }

  // no error 
  return 1; 
}

// --------------------------------------------------------- 

