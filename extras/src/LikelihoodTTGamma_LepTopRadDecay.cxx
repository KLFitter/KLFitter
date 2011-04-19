#include "LikelihoodTTGamma_LepTopRadDecay.h" 
#include "LikelihoodTopLeptonJets.h" 

// --------------------------------------------------------- 
KLFitter::LikelihoodTTGamma_LepTopRadDecay::LikelihoodTTGamma_LepTopRadDecay()
{
  // calls base class constructor
}

// --------------------------------------------------------- 
KLFitter::LikelihoodTTGamma_LepTopRadDecay::~LikelihoodTTGamma_LepTopRadDecay()
{
  // calls base class destructor
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTGamma_LepTopRadDecay::CalculateLorentzVectors(std::vector <double> parameters)
{
  KLFitter::LikelihoodTTGamma::CalculateLorentzVectors(parameters);

  // leptonic top 
  *(fParticlesModel->Parton(5)) += *(fParticlesModel->Photon(0));

  // no error 
  return 1; 
}
