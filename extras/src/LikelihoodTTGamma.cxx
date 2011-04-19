#include "LikelihoodTTGamma.h" 
#include "DetectorBase.h"
#include "ResolutionBase.h"

// --------------------------------------------------------- 
KLFitter::LikelihoodTTGamma::LikelihoodTTGamma() : KLFitter::LikelihoodTopLeptonJets::LikelihoodTopLeptonJets()
{
  // define model particles 
  this->DefineModelParticles(); 

  // define parameters 
  this->DefineParameters(); 

  // reset flag from base class constructor
  fFlagTopMassFixed = true; 
}

// --------------------------------------------------------- 
KLFitter::LikelihoodTTGamma::~LikelihoodTTGamma()
{
  // calls base class destructor
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTGamma::DefineModelParticles()
{
  // call base class constructor
  KLFitter::LikelihoodTopLeptonJets::DefineModelParticles();

  // add the photon
  TLorentzVector * lv = 0; 
  fParticlesModel->AddParticle(lv = new TLorentzVector(),
                               KLFitter::Particles::kPhoton,
                               "photon");
  // no error 
  return 1;
}

// --------------------------------------------------------- 
void KLFitter::LikelihoodTTGamma::DefineParameters()
{
  // check if base class constructor has already been called
  if (this->NParameters() != 17)
    KLFitter::LikelihoodTopLeptonJets::DefineParameters();

  // add photon energy
  this->AddParameter("energy photon", 0.0, 1000.0); // par 17 (the 18th parameter - counting starts well at 0)
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTGamma::AdjustParameterRanges()
{
  // call base class method
  KLFitter::LikelihoodTopLeptonJets::AdjustParameterRanges();

  // add parameter 17: energy of photon
  double nsigmas_photon = 2.0; 
  double E = (*fParticlesPermuted)->Photon(0)->E(); 
  double Emin = TMath::Max(0.001, E - nsigmas_photon * sqrt(E)); 
  double Emax  = E + nsigmas_photon * sqrt(E); 
  this->SetParameterRange(17, Emin, Emax); 

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTGamma::CalculateLorentzVectors(std::vector <double> parameters)
{
  KLFitter::LikelihoodTopLeptonJets::CalculateLorentzVectors(parameters);

  // photon
  TLorentzVector * vect = (*fParticlesPermuted)->Photon(0); 
  double E = parameters.at(17); 
  double px = vect->Px(); 
  double py = vect->Py(); 
  double pz = vect->Pz(); 
  double scale = E / vect->E(); 
  fParticlesModel->Photon(0)->SetPxPyPzE(scale * px, scale * py, scale * pz, E); 

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodTTGamma::LogLikelihood(std::vector <double> parameters)
{
  // calculate 4-vectors 
  this->CalculateLorentzVectors(parameters); 

  // get base likelihood from base class
  // this works also with one additional parameter, because the first
  // parameters are well preserved
  double logprob = KLFitter::LikelihoodTopLeptonJets::LogLikelihood(parameters);

  // add photon transfer function
  bool TFgood(true);
  logprob += log((*fDetector)->ResEnergyPhoton((*fParticlesPermuted)->Photon(0)->Eta())->p(parameters.at(17), (*fParticlesPermuted)->Photon(0)->E(), TFgood) );

  // return log of likelihood 
  return logprob; 
}

// --------------------------------------------------------- 
std::vector<double> KLFitter::LikelihoodTTGamma::GetInitialParametersWoNeutrinoPz()
{
  std::vector<double> values(KLFitter::LikelihoodTopLeptonJets::GetInitialParametersWoNeutrinoPz());

  // add initial value for the photon energy
  values.push_back((*fParticlesPermuted)->Photon(0)->E());

  // return the vector
  return values;
}
