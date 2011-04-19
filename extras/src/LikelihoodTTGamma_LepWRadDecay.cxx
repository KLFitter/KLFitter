#include "LikelihoodTTGamma_LepWRadDecay.h" 

// --------------------------------------------------------- 
KLFitter::LikelihoodTTGamma_LepWRadDecay::LikelihoodTTGamma_LepWRadDecay()
{
  // calls base class constructor
}

// --------------------------------------------------------- 
KLFitter::LikelihoodTTGamma_LepWRadDecay::~LikelihoodTTGamma_LepWRadDecay()
{
  // calls base class destructor
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTGamma_LepWRadDecay::CalculateLorentzVectors(std::vector <double> parameters)
{
  KLFitter::LikelihoodTTGamma::CalculateLorentzVectors(parameters);

  // the leptonic W
  *(fParticlesModel->Boson(1)) += *(fParticlesModel->Photon(0));

  // leptonic top 
  *(fParticlesModel->Parton(5)) += *(fParticlesModel->Photon(0));

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
std::vector<double> KLFitter::LikelihoodTTGamma_LepWRadDecay::GetNeutrinoPzSolutions()
{
  return this->CalculateNeutrinoPzSolutions((*fParticlesPermuted)->Photon(0));
}
