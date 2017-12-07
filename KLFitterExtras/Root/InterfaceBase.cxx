#include "KLFitterExtras/InterfaceBase.h" 
#include "KLFitter/Particles.h" 

// --------------------------------------------------------- 
KLFitter::InterfaceBase::InterfaceBase() :
  fParticles(0),
  fParticlesTruth(0),  
  fWeight(0)
{
}

// --------------------------------------------------------- 
KLFitter::InterfaceBase::~InterfaceBase()
{
  if (fParticles) 
    delete fParticles; 

  if (fParticlesTruth)
    delete fParticlesTruth; 
}

// --------------------------------------------------------- 

