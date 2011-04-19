#include "LikelihoodTTGamma_HadTopRadDecay.h" 
#include "LikelihoodTopLeptonJets.h" 

// --------------------------------------------------------- 
KLFitter::LikelihoodTTGamma_HadTopRadDecay::LikelihoodTTGamma_HadTopRadDecay()
{
  // calls base class constructor
}

// --------------------------------------------------------- 
KLFitter::LikelihoodTTGamma_HadTopRadDecay::~LikelihoodTTGamma_HadTopRadDecay()
{
  // calls base class destructor
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTGamma_HadTopRadDecay::CalculateLorentzVectors(std::vector <double> parameters)
{
  KLFitter::LikelihoodTTGamma::CalculateLorentzVectors(parameters);

  // hadronic top 
  *(fParticlesModel->Parton(4)) += *(fParticlesModel->Photon(0));

  // no error 
  return 1; 
}
