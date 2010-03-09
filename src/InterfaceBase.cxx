#include "InterfaceBase.h" 
#include <iostream> 

// --------------------------------------------------------- 
KLFitter::InterfaceBase::InterfaceBase()
{
	fParticles = 0; 
	fParticlesTruth = 0; 
	fLorentzVectorContainer = 0; 
}

// --------------------------------------------------------- 
KLFitter::InterfaceBase::~InterfaceBase()
{
	delete fParticles; 
	delete fParticlesTruth; 
}

// --------------------------------------------------------- 

