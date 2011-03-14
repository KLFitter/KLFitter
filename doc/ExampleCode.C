////////////////////////
// include these headers
////////////////////////
#include "Fitter.h"
#include "DetectorAtlas_7TeV.h"
#include "LikelihoodTopLeptonJets.h"
#include "PhysicsConstants.h"
#include "Particles.h"
#include "Permutations.h"
#include "TLorentzVector.h"

///////////////////////////
// before the event loop do
///////////////////////////

// create an instance of the fitter 
KLFitter::Fitter * myFitter = new KLFitter::Fitter(); 

// create an instance of the detector, which holds the information on the resolutions (transfer functions);
// it takes as an argument the folder which contains the parameter files for the transfer functions
KLFitter::DetectorBase * myDetector = new KLFitter::DetectorAtlas_7TeV("<path to transferfunctions/ttbar>"); 

// tell the fitter which detector to use
if (!myFitter->SetDetector(myDetector))
  return 0; 

// create an instance of the likelihood for ttbar->l+jets channel and customize it according to your needs
//
// DISCLAIMER: THE CUSTOMIZATION GIVEN HERE IS JUST AN EXAMPLE
//             YOU NEED TO ADAPT IT TO WHAT YOU NEED FOR YOUR ANALYSIS
//
KLFitter::LikelihoodTopLeptonJets * myLikelihood = new KLFitter::LikelihoodTopLeptonJets(); 
// set the lepton type for the hypothesis to be tested by the fitter
// kElectron or kMuon
myLikelihood->SetLeptonType(Fitter::LikelihoodTopLeptonJets::kElectron);
// if true: the top mass is constrained to the Breit-Wigner distribution around a fixed top mass value
myLikelihood->SetFlagTopMassFixed(false);
// set the central value for the fixed top mass for the case SetFlagTopMassFixed == true
myLikelihood->PhysicsConstants()->SetMassTop(172.5);
// if true: use a b-tagging veto to constrain the number of permutations
myLikelihood->SetFlagBTagging(false);
// set the cut value for the b-tagging veto
myLikelihood->SetCutBTag(5.85);
// if true: use the jet mass instead of the quark mass
myLikelihood->SetFlagUseJetMass(false);

// tell the fitter which likelihood to use
myFitter->SetLikelihood(myLikelihood);

////////////////////////
// in the event loop do
////////////////////////

// create an instance of the particles class filled with the particles to be fitted;
// here, you need to make sure that
// - the particles are in the range allowed by the transfer functions (eta and pt)
// - the energies and momenta are in GeV
// - be aware that *all* particles you're adding are considered in the fit
//   (many particles lead to many permutations to be considered and hence a long
//   running time and not necessarily good fitting results due to the many available
//   permutations)
// the arguments taken py AddParticle() are
// - TLorentzVector of the physics 4-momentum
// - detector eta for the evaluation of the transfer functions (for muons: just use the physics eta)
// - type of particle
// - an optional name of the particle
// - for jets: an optional b-tagging weight
KLFitter::Particles * myParticles = new KLFitter::Particles();
// add all jets like this (|eta| may not exceed 2.5):
TLorentzVector * vJet = new TLorentzVector((*Jet_Px)[iJet], (*Jet_Py)[iJet], (*Jet_Pz)[iJet], (*Jet_E)[iJet]));
myParticles->AddParticle(vJet, (*Jet_DetEta)[iJet], KLFitter::Particles::kParton, "", (*Jet_SV0_Weight)[iJet]);
// add all electrons like this (|eta| may not exceed 2.5):
TLorentzVector * vElectron = new TLorentzVector((*Electron_Px)[iElectron], (*Electron_Py)[iElectron], (*Electron_Pz)[iElectron], (*Electron_E)[iElectron]));
myParticles->AddParticle(vElectron, (*Electron_DetEta)[iElectron], KLFitter::Particles::kElectron);
// add all muons like this (|eta| may not exceed 2.5):
TLorentzVector * vMuon = new TLorentzVector((*Muon_Px)[iMuon], (*Muon_Py)[iMuon], (*Muon_Pz)[iMuon], (*Muon_E)[iMuon]);
myParticles->AddParticle(vMuon, (*Muon_Eta)[iMuon], KLFitter::Particles::kMuon);

// check that there are at least 4 jets
if (myParticles->NPartons() < 4)
  continue;
// for the electron hypothesis check that there is exactly one electron and no muon;
// (for the muon hypothesis require one muon and no electron)
if (myParticles->NElectrons() != 1 && myParticles->NMuons() != 0)
  continue;

// add the particles to the fitter
if (!myFitter->SetParticles(myParticles))
  return 0;       

// add the MET x and y components to the fitter
if (!myFitter->SetET_miss_XY(MET_Etx, MET_Ety))
  return 0;

// loop over all permutations
for (int iPerm = 0, nPerms(myFitter->Permutations()->NPermutations()); iPerm < nPerms; ++iPerm) {
  // perform the fit
  myFitter->Fit(iPerm); 

  // get the output from the fitter:
  // - the model particles
  KLFitter::Particles * myModelParticles = myFitter->Likelihood()->ParticlesModel();
  // a bit word with potential problems of the convergence of the fit
  // - the bit masks are defined in the Fitter class:
  //   * MinuitDidNotConvergeMask                 - Minuit fit did not converge
  //   * FitAbortedDueToNaNMask                   - fit was aborted due to a not-a-number value during the fit (typically remove it)
  //   * AtLeastOneFitParameterAtItsLimitMask     - fit converged, but at least one parameter is at its allowed limit (typically keep it)
  //   * InvalidTransferFunctionAtConvergenceMask - invalid use of the transfer functions at the convergence point (typically remove it)
  unsigned int ConvergenceStatusBitWord = myFitter->ConvergenceStatus();
  // in order to check if a certain had a specific problem, just check the bit mask as in this example:
  bool MinuitDidNotConverge = (ConvergenceStatusBitWord & Fitter->MinuitDidNotConvergeMask) == 0;
  // etc.
  // get the log(likelihood) value from the fit
  double LogLikelihood = myFitter->Likelihood()->LogLikelihood(myFitter->Likelihood()->GetBestFitParameters());
  // get the event probability from the fit
  double EventProbability = exp(myFitter->Likelihood()->LogEventProbability());
  // get the fit parameters
  std::vector<double> Parameters = myFitter->Likelihood()->GetBestFitParameters();
  std::vector<double> ParameterErrors = myFitter->Likelihood()->GetBestFitParameterErrors();
 }
}
