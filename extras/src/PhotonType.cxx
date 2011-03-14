#include "PhotonType.h" 
#include "Particles.h"

#include <iostream>

// --------------------------------------------------------- 
KLFitter::PhotonType::PhotonType()
{
  deltaTopMass = 8.0;
  deltaWMass   = 8.0;
  isNotClassified  = false;
  isRadTopProd     = false;
  isHadTopRadDecay = false;
  isLepTopRadDecay = false;
  isHadWRadDecay   = false;
  isLepWRadDecay   = false;
  fTruthParticles = 0x0;
  fPhysicsConstants = 0x0
    ;}

// --------------------------------------------------------- 
int KLFitter::PhotonType::Classify()
{
  // clear info
  isNotClassified  = false;
  isRadTopProd     = false;
  isHadTopRadDecay = false;
  isLepTopRadDecay = false;
  isHadWRadDecay   = false;
  isLepWRadDecay   = false;

  // make composite particles
  MakeCompositeParticles();

  // call the classification methods
  if (ClassifyRadTopProd())
    isRadTopProd = true;

  if (ClassifyHadTopRadDecay())
    isHadTopRadDecay = true;

  if (ClassifyLepTopRadDecay())
    isLepTopRadDecay = true;

  if (ClassifyHadWRadDecay())
    isHadWRadDecay = true;

  if (ClassifyLepWRadDecay())
    isLepWRadDecay = true;

  isNotClassified = !isRadTopProd && !isHadTopRadDecay && !isLepTopRadDecay && !isHadWRadDecay && !isLepWRadDecay;

  // no error
  return 1;
}

// --------------------------------------------------------- 
bool KLFitter::PhotonType::ClassifyRadTopProd() {
  if (eqTopM(htop.M(), deltaTopMass) &&
      btTopM(htop_g.M(), deltaTopMass) &&
      eqTopM(ltop.M(), deltaTopMass) &&
      btTopM(ltop_g.M(), deltaTopMass) &&
      eqWM(hadW.M(), deltaWMass) &&
      btWM(hadW_g.M(), deltaWMass) &&
      eqWM(lepW.M(), deltaWMass) &&
      btWM(lepW_g.M(), deltaWMass))
    return true;

  return false;
}

// --------------------------------------------------------- 
bool KLFitter::PhotonType::ClassifyHadTopRadDecay() {
  if (ltTopM(htop.M(), deltaTopMass) &&
      eqTopM(htop_g.M(), deltaTopMass) &&
      eqTopM(ltop.M(), deltaTopMass) &&
      btTopM(ltop_g.M(), deltaTopMass) &&
      eqWM(hadW.M(), deltaWMass) &&
      btWM(hadW_g.M(), deltaWMass) &&
      eqWM(lepW.M(), deltaWMass) &&
      btWM(lepW_g.M(), deltaWMass))
    return true;

  return false;
}

// --------------------------------------------------------- 
bool KLFitter::PhotonType::ClassifyLepTopRadDecay() {
  if (eqTopM(htop.M(), deltaTopMass) &&
      btTopM(htop_g.M(), deltaTopMass) &&
      ltTopM(ltop.M(), deltaTopMass) &&
      eqTopM(ltop_g.M(), deltaTopMass) &&
      eqWM(hadW.M(), deltaWMass) &&
      btWM(hadW_g.M(), deltaWMass) &&
      eqWM(lepW.M(), deltaWMass) &&
      btWM(lepW_g.M(), deltaWMass))
    return true;

  return false;
}

// --------------------------------------------------------- 
bool KLFitter::PhotonType::ClassifyHadWRadDecay() {
  if (ltTopM(htop.M(), deltaTopMass) &&
      eqTopM(htop_g.M(), deltaTopMass) &&
      eqTopM(ltop.M(), deltaTopMass) &&
      btTopM(ltop_g.M(), deltaTopMass) &&
      ltWM(hadW.M(), deltaWMass) &&
      eqWM(hadW_g.M(), deltaWMass) &&
      eqWM(lepW.M(), deltaWMass) &&
      btWM(lepW_g.M(), deltaWMass))
    return true;

  return false;
}

// --------------------------------------------------------- 
bool KLFitter::PhotonType::ClassifyLepWRadDecay() {
  if (eqTopM(htop.M(), deltaTopMass) &&
      btTopM(htop_g.M(), deltaTopMass) &&
      ltTopM(ltop.M(), deltaTopMass) &&
      eqTopM(ltop_g.M(), deltaTopMass) &&
      eqWM(hadW.M(), deltaWMass) &&
      btWM(hadW_g.M(), deltaWMass) &&
      ltWM(lepW.M(), deltaWMass) &&
      eqWM(lepW_g.M(), deltaWMass))
    return true;

  return false;
}

// --------------------------------------------------------- 
void KLFitter::PhotonType::MakeCompositeParticles() {
  KLFitter::Particles * p = *fTruthParticles;

  hadW   = *(p->Parton(2)) + *(p->Parton(3));
  hadW_g = hadW + *(p->Photon(0));

  lepW   = *(p->Electron(0)) + *(p->Neutrino(0));
  lepW_g = lepW + *(p->Photon(0));

  htop   = *(p->Parton(0)) + hadW;
  htop_g = *(p->Parton(0)) + hadW_g;

  ltop   = *(p->Parton(1)) + lepW;
  ltop_g = *(p->Parton(1)) + lepW_g;
}
