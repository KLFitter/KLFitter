#include "LikelihoodTTHLeptonJets.h" 
#include "ResolutionBase.h"
#include "Particles.h"
#include "Permutations.h"
#include "PhysicsConstants.h"
#include "DetectorBase.h"

#include <iostream> 
#include <algorithm> 

#include <BAT/BCMath.h> 

// --------------------------------------------------------- 
KLFitter::LikelihoodTTHLeptonJets::LikelihoodTTHLeptonJets() : KLFitter::LikelihoodBase::LikelihoodBase()
                                                             , fFlagTopMassFixed(false)
							     , fFlagHiggsMassFixed(false)
                                                             , fFlagUseJetMass(false)
                                                             , ETmiss_x(0.)
                                                             , ETmiss_y(0.)
                                                             , SumET(0.)
                                                             , fTypeLepton(kElectron)
                                                             , fTFgood(true)
{


  // define model particles 
  this->DefineModelParticles(); 

  // define parameters 
  // this->DefineParameters(); 
}

// --------------------------------------------------------- 
KLFitter::LikelihoodTTHLeptonJets::~LikelihoodTTHLeptonJets()
{
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTHLeptonJets::SetET_miss_XY_SumET(double etx, double ety, double sumet)
{
  // set missing ET x and y component and the SumET
  ETmiss_x = etx;
  ETmiss_y = ety;
  SumET = sumet;

  // no error
  return 1;
}

// --------------------------------------------------------- 
void KLFitter::LikelihoodTTHLeptonJets::SetLeptonType(LeptonType leptontype)
{
  if (leptontype != kElectron && leptontype != kMuon)
    {
      std::cout << "KLFitter::SetLeptonTyp(). Warning: lepton type not defined. Set electron as lepton type." << std::endl; 
      fTypeLepton = kElectron;
    }
  else
    fTypeLepton = leptontype; 

  // define model particles 
  DefineModelParticles(); 
}

// --------------------------------------------------------- 
void KLFitter::LikelihoodTTHLeptonJets::SetLeptonType(int leptontype)
{
  if (leptontype != 1 && leptontype != 2)
    {
      std::cout << "KLFitter::SetLeptonTyp(). Warning: lepton type not defined. Set electron as lepton type." << std::endl; 
      leptontype = 1;
    }

  if (leptontype == 1)
    SetLeptonType(kElectron);
  else if (leptontype == 2)
    SetLeptonType(kMuon);
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTHLeptonJets::DefineModelParticles()
{
  // check if model particles and lorentz vector container exist and delete
  if (fParticlesModel) {
    delete fParticlesModel; 
    fParticlesModel = 0;
  }

  // create the particles of the model 
  fParticlesModel = new KLFitter::Particles(); 

  // add model particles
  //create dummy TLorentzVector
  TLorentzVector * dummy = new TLorentzVector(0,0,0,0); // 4-vector
  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kParton, // type 
                               "hadronic b quark",           // name 
                               0,                            // index of corresponding particle 
                               KLFitter::Particles::kB);     // b jet (truth)

  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kParton, 
                               "leptonic b quark",
                               1,                            // index of corresponding particle 
                               KLFitter::Particles::kB);     // b jet (truth)

  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kParton,
                               "light quark 1",
                               2,                            // index of corresponding particle 
                               KLFitter::Particles::kLight); // light jet (truth)

  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kParton,
                               "light quark 2",
                               3,                            // index of corresponding particle
                               KLFitter::Particles::kLight); // light jet (truth)
 
  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kParton, // type 
                               "Higgs b quark 1",           // name 
                               4,                            // index of corresponding particle 
                               KLFitter::Particles::kB);     // b jet (truth)

  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kParton, 
                               "Higgs b quark 2",
                               5,                            // index of corresponding particle 
                               KLFitter::Particles::kB);     // b jet (truth)
        
  if (fTypeLepton == kElectron) {
    fParticlesModel->AddParticle(dummy,
                                 KLFitter::Particles::kElectron,
                                 "electron"); 
  }
  else if (fTypeLepton == kMuon) {
    fParticlesModel->AddParticle(dummy,
                                 KLFitter::Particles::kMuon,
                                 "muon"); 
  }

  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kNeutrino, 
                               "neutrino"); 
  
  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kBoson, 
                               "hadronic W"); 
  
  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kBoson,
                               "leptonic W"); 

  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kParton,
                               "hadronic top");

  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kParton,
                               "leptonic top");
 
  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kBoson,
                               "Higgs"); 


  //free memory
  delete dummy; 

  // no error 
  return 1;
}

// --------------------------------------------------------- 
void KLFitter::LikelihoodTTHLeptonJets::DefineParameters()
{
  // add parameters of model
  AddParameter("energy hadronic b",       fPhysicsConstants->MassBottom(), 1000.0); // parBhadE
  AddParameter("energy leptonic b",       fPhysicsConstants->MassBottom(), 1000.0); // parBlepE
  AddParameter("energy light quark 1",    0.0, 1000.0);                             // parLQ1E
  AddParameter("energy light quark 2",    0.0, 1000.0);                             // parLQ2E
  AddParameter("energy lepton",           0.0, 1000.0);                             // parLepE
  AddParameter("p_x neutrino",        -1000.0, 1000.0);                             // parNuPx
  AddParameter("p_y neutrino",        -1000.0, 1000.0);                             // parNuPy
  AddParameter("p_z neutrino",        -1000.0, 1000.0);                             // parNuPz
  AddParameter("top mass",              100.0, 1000.0);                             // parTopM
  AddParameter("energy Higgs b quark 1",  fPhysicsConstants->MassBottom(), 1000.0); // parBHiggs1E
  AddParameter("energy Higgs b quark 2",  fPhysicsConstants->MassBottom(), 1000.0); // parBHiggs2E
  if (fFlagHiggsMassFixed)  AddParameter("Higgs mass",              100.0, 1000.0);   // parHiggsM
  //++++++++++++++++//
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTHLeptonJets::CalculateLorentzVectors(std::vector <double> const& parameters)
{
  static double scale;
  static double whad_fit_e;
  static double whad_fit_px;
  static double whad_fit_py;
  static double whad_fit_pz;
  static double wlep_fit_e;
  static double wlep_fit_px;
  static double wlep_fit_py;
  static double wlep_fit_pz;
  static double thad_fit_e;
  static double thad_fit_px;
  static double thad_fit_py;
  static double thad_fit_pz;
  static double tlep_fit_e;
  static double tlep_fit_px;
  static double tlep_fit_py;
  static double tlep_fit_pz;

  //++++++++++++++++//
  static double Higgs_fit_e;
  static double Higgs_fit_px;
  static double Higgs_fit_py;
  static double Higgs_fit_pz;
  //++++++++++++++++//


  // hadronic b quark 
  bhad_fit_e = parameters[parBhadE]; 
  scale = sqrt(bhad_fit_e*bhad_fit_e - bhad_meas_m*bhad_meas_m) / bhad_meas_p;
  bhad_fit_px = scale * bhad_meas_px;
  bhad_fit_py = scale * bhad_meas_py;
  bhad_fit_pz = scale * bhad_meas_pz;

  // leptonic b quark 
  blep_fit_e = parameters[parBlepE]; 
  scale = sqrt(blep_fit_e*blep_fit_e - blep_meas_m*blep_meas_m) / blep_meas_p;
  blep_fit_px = scale * blep_meas_px;
  blep_fit_py = scale * blep_meas_py;
  blep_fit_pz = scale * blep_meas_pz;

  // light quark 1 
  lq1_fit_e = parameters[parLQ1E]; 
  scale = sqrt(lq1_fit_e*lq1_fit_e - lq1_meas_m*lq1_meas_m) / lq1_meas_p;
  lq1_fit_px = scale * lq1_meas_px;
  lq1_fit_py = scale * lq1_meas_py;
  lq1_fit_pz = scale * lq1_meas_pz;

  // light quark 2 
  lq2_fit_e = parameters[parLQ2E]; 
  scale = sqrt(lq2_fit_e*lq2_fit_e - lq2_meas_m*lq2_meas_m) / lq2_meas_p;
  lq2_fit_px  = scale * lq2_meas_px;
  lq2_fit_py  = scale * lq2_meas_py;
  lq2_fit_pz  = scale * lq2_meas_pz;

 
  //++++++++++++++++//
  // Higgs b quark 1 
  BHiggs1_fit_e = parameters[parBHiggs1E]; 
  scale = sqrt(BHiggs1_fit_e*BHiggs1_fit_e - BHiggs1_meas_m*BHiggs1_meas_m) / BHiggs1_meas_p;
  BHiggs1_fit_px = scale * BHiggs1_meas_px;
  BHiggs1_fit_py = scale * BHiggs1_meas_py;
  BHiggs1_fit_pz = scale * BHiggs1_meas_pz;


  // Higgs b quark 2
  BHiggs2_fit_e = parameters[parBHiggs2E]; 
  scale = sqrt(BHiggs2_fit_e*BHiggs2_fit_e - BHiggs2_meas_m*BHiggs2_meas_m) / BHiggs2_meas_p;
  BHiggs2_fit_px = scale * BHiggs2_meas_px;
  BHiggs2_fit_py = scale * BHiggs2_meas_py;
  BHiggs2_fit_pz = scale * BHiggs2_meas_pz;
  //++++++++++++++++//
 

  // lepton
  lep_fit_e = parameters[parLepE];
  scale = lep_fit_e / lep_meas_e;
  lep_fit_px = scale * lep_meas_px;
  lep_fit_py = scale * lep_meas_py;
  lep_fit_pz = scale * lep_meas_pz;

  // neutrino 
  nu_fit_px = parameters[parNuPx];
  nu_fit_py = parameters[parNuPy];
  nu_fit_pz = parameters[parNuPz];
  nu_fit_e  = sqrt(nu_fit_px*nu_fit_px + nu_fit_py*nu_fit_py + nu_fit_pz*nu_fit_pz);

  // hadronic W 
  whad_fit_e  = lq1_fit_e +lq2_fit_e;
  whad_fit_px = lq1_fit_px+lq2_fit_px;
  whad_fit_py = lq1_fit_py+lq2_fit_py;
  whad_fit_pz = lq1_fit_pz+lq2_fit_pz;
  whad_fit_m = sqrt(whad_fit_e*whad_fit_e - (whad_fit_px*whad_fit_px + whad_fit_py*whad_fit_py + whad_fit_pz*whad_fit_pz));
  
  // leptonic W 
  wlep_fit_e  = lep_fit_e +nu_fit_e;
  wlep_fit_px = lep_fit_px+nu_fit_px;
  wlep_fit_py = lep_fit_py+nu_fit_py;
  wlep_fit_pz = lep_fit_pz+nu_fit_pz;
  wlep_fit_m = sqrt(wlep_fit_e*wlep_fit_e - (wlep_fit_px*wlep_fit_px + wlep_fit_py*wlep_fit_py + wlep_fit_pz*wlep_fit_pz));

  // hadronic top 
  thad_fit_e = whad_fit_e+bhad_fit_e;
  thad_fit_px = whad_fit_px+bhad_fit_px;
  thad_fit_py = whad_fit_py+bhad_fit_py;
  thad_fit_pz = whad_fit_pz+bhad_fit_pz;
  thad_fit_m = sqrt(thad_fit_e*thad_fit_e - (thad_fit_px*thad_fit_px + thad_fit_py*thad_fit_py + thad_fit_pz*thad_fit_pz));

  // leptonic top 
  tlep_fit_e = wlep_fit_e+blep_fit_e;
  tlep_fit_px = wlep_fit_px+blep_fit_px;
  tlep_fit_py = wlep_fit_py+blep_fit_py;
  tlep_fit_pz = wlep_fit_pz+blep_fit_pz;
  tlep_fit_m = sqrt(tlep_fit_e*tlep_fit_e - (tlep_fit_px*tlep_fit_px + tlep_fit_py*tlep_fit_py + tlep_fit_pz*tlep_fit_pz));


  //++++++++++++++++//
  // Higgs 
  Higgs_fit_e  = BHiggs1_fit_e  + BHiggs2_fit_e;
  Higgs_fit_px = BHiggs1_fit_px + BHiggs2_fit_px;
  Higgs_fit_py = BHiggs1_fit_py + BHiggs2_fit_py;
  Higgs_fit_pz = BHiggs1_fit_pz + BHiggs2_fit_pz;
  Higgs_fit_m  = sqrt(Higgs_fit_e*Higgs_fit_e - (Higgs_fit_px*Higgs_fit_px + Higgs_fit_py*Higgs_fit_py + Higgs_fit_pz*Higgs_fit_pz));
  //++++++++++++++++//


  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTHLeptonJets::Initialize()
{
  // error code 
  int err = 1; 

  // save the current permuted particles
  err *= SavePermutedParticles();

  // save the corresponding resolution functions
  err *= SaveResolutionFunctions();

  this->DefineParameters(); 

  // adjust parameter ranges 
  err *= AdjustParameterRanges(); 

  // set initial values
  // (only for Markov chains - initial parameters for other minimisation methods are set in Fitter.cxx)
  SetInitialParameters( GetInitialParameters() );       

  // return error code 
  return err; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTHLeptonJets::RemoveInvariantParticlePermutations()
{
  // error code 
  int err = 1; 

  // remove the permutation from the second and the third jet
  KLFitter::Particles::ParticleType ptype = KLFitter::Particles::kParton;
  std::vector<int> indexVector_Jets;
  indexVector_Jets.push_back(2);
  indexVector_Jets.push_back(3);
  err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Jets); 
        
  //++++++++++++++++//
  // remove the permutation from the two Higgs b jets
  KLFitter::Particles::ParticleType ptypeBB = KLFitter::Particles::kParton;
  std::vector<int> indexVector_JetsBB;
  indexVector_JetsBB.push_back(4);
  indexVector_JetsBB.push_back(5);
  err *= (*fPermutations)->InvariantParticlePermutations(ptypeBB, indexVector_JetsBB); 
   
  //remove invariant jet permutations of notevent jets
  KLFitter::Particles* particles = (*fPermutations)->Particles();
  indexVector_Jets.clear();
  for (int iPartons = 6; iPartons < particles->NPartons(); iPartons++)
    indexVector_Jets.push_back(iPartons);
  err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Jets);
  //++++++++++++++++//


  // remove the permutation from the other lepton
  if (fTypeLepton == kElectron)
    {
      ptype = KLFitter::Particles::kMuon;
      std::vector<int> indexVector_Muons;
      for (int iMuon = 0; iMuon < particles->NMuons(); iMuon++)
        indexVector_Muons.push_back(iMuon);
      err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Muons); 
    }
  else if (fTypeLepton == kMuon)
    {
      ptype = KLFitter::Particles::kElectron;
      std::vector<int> indexVector_Electrons;
      for (int iElectron = 0; iElectron < particles->NElectrons(); iElectron++)
        indexVector_Electrons.push_back(iElectron);
      err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Electrons); 
    }

  // return error code 
  return err; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTHLeptonJets::RemoveForbiddenParticlePermutations()
{
  // error code 
  int err = 1; 

  // only in b-tagging type kVetoNoFit
  if (!((fBTagMethod == kVetoNoFit)||(fBTagMethod == kVetoNoFitLight)||(fBTagMethod == kVetoNoFitBoth)))
    return err;

  // remove all permutations where a b-tagged jet is in the position of a model light quark
  KLFitter::Particles * particles = (*fPermutations)->Particles();
  int nPartons = particles->NPartons();

  KLFitter::Particles * particlesModel = fParticlesModel;
  int nPartonsModel = particlesModel->NPartons();
  for (int iParton(0); iParton < nPartons; ++iParton) {
    bool isBtagged = particles->IsBTagged(iParton);

    for (int iPartonModel(0); iPartonModel < nPartonsModel; ++iPartonModel) {

      KLFitter::Particles::TrueFlavorType trueFlavor = particlesModel->TrueFlavor(iPartonModel);
      if ((fBTagMethod == kVetoNoFit)&&((!isBtagged)||(trueFlavor != KLFitter::Particles::kLight)))
        continue;
      if ((fBTagMethod == kVetoNoFitLight)&&((isBtagged)||(trueFlavor != KLFitter::Particles::kB)))
        continue;
      if ((fBTagMethod == kVetoNoFitBoth)&&(((isBtagged)&&(trueFlavor != KLFitter::Particles::kLight))||((!isBtagged)&&(trueFlavor != KLFitter::Particles::kB))))
        continue;

      err *= (*fPermutations)->RemoveParticlePermutations(KLFitter::Particles::kParton, iParton, iPartonModel);
    }
  }
        
  // return error code 
  return err; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTHLeptonJets::AdjustParameterRanges()
{
  // adjust limits 
  double nsigmas_jet = 7.0; 
  double nsigmas_lepton = 2.0; 

  double E = (*fParticlesPermuted)->Parton(0)->E(); 
  double m = fPhysicsConstants->MassBottom(); 
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(0)->M()); 
  double Emin = std::max(m, E - nsigmas_jet* sqrt(E)); 
  double Emax  = E + nsigmas_jet* sqrt(E); 
  SetParameterRange(parBhadE, Emin, Emax); 

  E = (*fParticlesPermuted)->Parton(1)->E(); 
  m = fPhysicsConstants->MassBottom(); 
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(1)->M()); 
  Emin = std::max(m, E - nsigmas_jet* sqrt(E)); 
  Emax  = E + nsigmas_jet* sqrt(E); 
  SetParameterRange(parBlepE, Emin, Emax); 

  E = (*fParticlesPermuted)->Parton(2)->E(); 
  m = 0.001;
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(2)->M()); 
  Emin = std::max(m, E - nsigmas_jet* sqrt(E)); 
  Emax  = E + nsigmas_jet* sqrt(E); 
  SetParameterRange(parLQ1E, Emin, Emax); 

  E = (*fParticlesPermuted)->Parton(3)->E(); 
  m = 0.001;
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(3)->M()); 
  Emin = std::max(m, E - nsigmas_jet* sqrt(E)); 
  Emax  = E + nsigmas_jet* sqrt(E); 
  SetParameterRange(parLQ2E, Emin, Emax); 

  //++++++++++++++++//
  E = (*fParticlesPermuted)->Parton(4)->E(); 
  m = fPhysicsConstants->MassBottom(); 
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(4)->M()); 
  Emin = std::max(m, E - nsigmas_jet* sqrt(E)); 
  Emax  = E + nsigmas_jet* sqrt(E); 
  SetParameterRange(parBHiggs1E, Emin, Emax); 

  E = (*fParticlesPermuted)->Parton(5)->E(); 
  m = fPhysicsConstants->MassBottom(); 
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(5)->M()); 
  Emin = std::max(m, E - nsigmas_jet* sqrt(E)); 
  Emax  = E + nsigmas_jet* sqrt(E); 
  SetParameterRange(parBHiggs2E, Emin, Emax); 
  //++++++++++++++++//

  if (fTypeLepton == kElectron)
    {
      E = (*fParticlesPermuted)->Electron(0)->E();
      Emin = std::max(0.001, E - nsigmas_lepton* sqrt(E)); 
      Emax  = E + nsigmas_lepton* sqrt(E);  
    }
  else if (fTypeLepton == kMuon)
    {
      E = (*fParticlesPermuted)->Muon(0)->E(); 
      double sintheta= sin((*fParticlesPermuted)->Muon(0)->Theta());
      double sigrange=nsigmas_lepton* (E*E*sintheta);
      Emin=std::max(0.001,E -sigrange);
      Emax=E +sigrange;
    }
  SetParameterRange(parLepE, Emin, Emax); 

  // note: this is hard-coded in the momement 

  SetParameterRange(parNuPx, ETmiss_x-100.0, ETmiss_x+100);
  SetParameterRange(parNuPy, ETmiss_y-100.0, ETmiss_y+100);

  if (fFlagTopMassFixed)
    SetParameterRange(parTopM, fPhysicsConstants->MassTop(), fPhysicsConstants->MassTop()); 

  //++++++++++++++++//
  if (fFlagHiggsMassFixed)
    SetParameterRange(parHiggsM, fPhysicsConstants->MassHiggs(), fPhysicsConstants->MassHiggs());
  //++++++++++++++++//

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodTTHLeptonJets::LogLikelihood(const std::vector<double> & parameters)
{
  // calculate 4-vectors 
  CalculateLorentzVectors(parameters); 

  // define log of likelihood 
  double logprob(0.); 

  // temporary flag for a safe use of the transfer functions
  bool TFgoodTmp(true);

  // jet energy resolution terms 
  logprob += log( fResEnergyBhad->p(bhad_fit_e, bhad_meas_e, TFgoodTmp) );
  if (!TFgoodTmp) fTFgood = false;

  logprob += log( fResEnergyBlep->p(blep_fit_e, blep_meas_e, TFgoodTmp) );
  if (!TFgoodTmp) fTFgood = false;

  logprob += log( fResEnergyLQ1->p(lq1_fit_e, lq1_meas_e, TFgoodTmp) );
  if (!TFgoodTmp) fTFgood = false;

  logprob += log( fResEnergyLQ2->p(lq2_fit_e, lq2_meas_e, TFgoodTmp) );
  if (!TFgoodTmp) fTFgood = false;

  //++++++++++++++++//
  logprob += log( fResEnergyBHiggs1->p(BHiggs1_fit_e, BHiggs1_meas_e, TFgoodTmp) );
  if (!TFgoodTmp) fTFgood = false;

  logprob += log( fResEnergyBHiggs2->p(BHiggs2_fit_e, BHiggs2_meas_e, TFgoodTmp) );
  if (!TFgoodTmp) fTFgood = false;
  //++++++++++++++++//


  // lepton energy resolution terms 
  if (fTypeLepton == kElectron)
    logprob += log( fResLepton->p(lep_fit_e, lep_meas_e, TFgoodTmp) );
  else if (fTypeLepton == kMuon)
    logprob += log( fResLepton->p(lep_fit_e* lep_meas_sintheta, lep_meas_pt, TFgoodTmp) );
  if (!TFgoodTmp) fTFgood = false;

  // neutrino px and py
  logprob += log( fResMET->p(nu_fit_px, ETmiss_x, TFgoodTmp, SumET) );
  if (!TFgoodTmp) fTFgood = false;

  logprob += log( fResMET->p(nu_fit_py, ETmiss_y, TFgoodTmp, SumET) );
  if (!TFgoodTmp) fTFgood = false;

  // physics constants
  double massW = fPhysicsConstants->MassW();
  double gammaW = fPhysicsConstants->GammaW();
  double gammaTop = fPhysicsConstants->GammaTop();
 
  //++++++++++++++++//
  double gammaHiggs = fPhysicsConstants->GammaHiggs();
  //++++++++++++++++//

  // Breit-Wigner of hadronically decaying W-boson
  logprob += BCMath::LogBreitWignerRel(whad_fit_m, massW, gammaW);

  // Breit-Wigner of leptonically decaying W-boson
  logprob += BCMath::LogBreitWignerRel(wlep_fit_m, massW, gammaW);

  // Breit-Wigner of hadronically decaying top quark
  logprob += BCMath::LogBreitWignerRel(thad_fit_m, parameters[parTopM], gammaTop); 
        
  // Breit-Wigner of leptonically decaying top quark
  logprob += BCMath::LogBreitWignerRel(tlep_fit_m, parameters[parTopM], gammaTop); 

  //++++++++++++++++//
  // Breit-Wigner of Higgs decaying into 2 b-quark
  if (fFlagHiggsMassFixed) logprob += BCMath::LogBreitWignerRel(Higgs_fit_m, parameters[parHiggsM], gammaHiggs); 

  //++++++++++++++++//

  // return log of likelihood 
  return logprob; 
}

// --------------------------------------------------------- 
std::vector<double> KLFitter::LikelihoodTTHLeptonJets::GetInitialParameters()
{
  std::vector<double> values(GetInitialParametersWoNeutrinoPz());

  // check second neutrino solution
  std::vector<double> neutrino_pz_solutions = GetNeutrinoPzSolutions();
  if (int(neutrino_pz_solutions.size()) == 1)
    values[parNuPz] = neutrino_pz_solutions[0]; 
  else if(int(neutrino_pz_solutions.size()) == 2)
    {
      double sol1, sol2; 
      values[parNuPz] = neutrino_pz_solutions[0]; 
      sol1 = LogLikelihood(values); 
      values[parNuPz] = neutrino_pz_solutions[1]; 
      sol2 = LogLikelihood(values); 

      if (sol1 > sol2)
        values[parNuPz] = neutrino_pz_solutions[0]; 
    }

  return values;
}

// --------------------------------------------------------- 
std::vector<double> KLFitter::LikelihoodTTHLeptonJets::GetInitialParametersWoNeutrinoPz()
{
  std::vector<double> values(GetNParameters());

  // energies of the quarks
  values[parBhadE] = bhad_meas_e;
  values[parBlepE] = blep_meas_e;
  values[parLQ1E]  = lq1_meas_e;
  values[parLQ2E]  = lq2_meas_e;
  //++++++++++++++++//
  values[parBHiggs1E]  = BHiggs1_meas_e;
  values[parBHiggs2E]  = BHiggs2_meas_e;
  //++++++++++++++++//



  // energy of the lepton
  if (fTypeLepton == kElectron)
    values[parLepE] = (*fParticlesPermuted)->Electron(0)->E();
  else if (fTypeLepton == kMuon)
    values[parLepE] = (*fParticlesPermuted)->Muon(0)->E();

  // missing px and py
  values[parNuPx] = ETmiss_x;
  values[parNuPy] = ETmiss_y;

  // pz of the neutrino
  values[parNuPz] = 0.;

  // top mass
  double mtop = ( *(*fParticlesPermuted)->Parton(0) + *(*fParticlesPermuted)->Parton(2) + *(*fParticlesPermuted)->Parton(3) ).M(); 
  if (mtop < GetParameter(parTopM)->GetLowerLimit())
    mtop = GetParameter(parTopM)->GetLowerLimit(); 
  else if (mtop > GetParameter(parTopM)->GetUpperLimit())
    mtop = GetParameter(parTopM)->GetUpperLimit(); 
  values[parTopM] = mtop;

  //++++++++++++++++//
  // Higgs mass
  if (fFlagHiggsMassFixed){
	double mhiggs = ( *(*fParticlesPermuted)->Parton(4) + *(*fParticlesPermuted)->Parton(5)).M(); 
	if (mhiggs < GetParameter(parHiggsM)->GetLowerLimit())
	mhiggs = GetParameter(parHiggsM)->GetLowerLimit(); 
	else if (mhiggs > GetParameter(parHiggsM)->GetUpperLimit())
	mhiggs = GetParameter(parHiggsM)->GetUpperLimit(); 
	values[parHiggsM] = mhiggs;
  }
  //++++++++++++++++//


  // return the vector
  return values;
}

// --------------------------------------------------------- 
std::vector<double> KLFitter::LikelihoodTTHLeptonJets::GetNeutrinoPzSolutions() {
  return CalculateNeutrinoPzSolutions();
}

// --------------------------------------------------------- 
std::vector<double> KLFitter::LikelihoodTTHLeptonJets::CalculateNeutrinoPzSolutions(TLorentzVector* additionalParticle)
{
  std::vector<double> pz;

  KLFitter::PhysicsConstants constants;
  // electron mass
  double mE = 0.;

  double px_c = 0.0;
  double py_c = 0.0;
  double pz_c = 0.0;
  double Ec = 0.0; 

  if (fTypeLepton == kElectron)
    {
      px_c = (*fParticlesPermuted)->Electron(0)->Px();
      py_c = (*fParticlesPermuted)->Electron(0)->Py();
      pz_c = (*fParticlesPermuted)->Electron(0)->Pz();
      Ec = (*fParticlesPermuted)->Electron(0)->E();
    }
  else if (fTypeLepton == kMuon)
    {
      px_c = (*fParticlesPermuted)->Muon(0)->Px();
      py_c = (*fParticlesPermuted)->Muon(0)->Py();
      pz_c = (*fParticlesPermuted)->Muon(0)->Pz();
      Ec = (*fParticlesPermuted)->Muon(0)->E();
    }

  // add additional particle to "charged lepton" 4-vector
  if (additionalParticle) {
    px_c += additionalParticle->Px();
    py_c += additionalParticle->Py();
    pz_c += additionalParticle->Pz();
    Ec += additionalParticle->E();
  }

  double px_nu = ETmiss_x;
  double py_nu = ETmiss_y;
  double alpha = constants.MassW()*constants.MassW() - mE*mE + 2*(px_c*px_nu + py_c*py_nu);

  double a = pz_c*pz_c - Ec*Ec;
  double b = alpha* pz_c;
  double c = - Ec*Ec* (px_nu*px_nu + py_nu*py_nu) + alpha*alpha/4.;

  double discriminant = b*b - 4*a*c;
  if (discriminant < 0.)
    return pz;

  double pz_offset = - b / (2*a);

  double squareRoot = sqrt(discriminant);
  if (squareRoot < 1.e-6)
    pz.push_back(pz_offset);
  else {
    pz.push_back(pz_offset + squareRoot / (2*a));
    pz.push_back(pz_offset - squareRoot / (2*a));
  }

  return pz;
}

// --------------------------------------------------------- 
bool KLFitter::LikelihoodTTHLeptonJets::NoTFProblem(std::vector<double> parameters) {
  fTFgood = true;
  this->LogLikelihood(parameters);
  return fTFgood;
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTHLeptonJets::SavePermutedParticles() {
  bhad_meas_e      = (*fParticlesPermuted)->Parton(0)->E();
  bhad_meas_deteta = (*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kParton);
  bhad_meas_px     = (*fParticlesPermuted)->Parton(0)->Px();
  bhad_meas_py     = (*fParticlesPermuted)->Parton(0)->Py();
  bhad_meas_pz     = (*fParticlesPermuted)->Parton(0)->Pz();
  bhad_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(0)->M(), fPhysicsConstants->MassBottom(), bhad_meas_px, bhad_meas_py, bhad_meas_pz, bhad_meas_e);
  bhad_meas_p      = sqrt(bhad_meas_e*bhad_meas_e - bhad_meas_m*bhad_meas_m);

  blep_meas_e      = (*fParticlesPermuted)->Parton(1)->E();
  blep_meas_deteta = (*fParticlesPermuted)->DetEta(1, KLFitter::Particles::kParton);
  blep_meas_px     = (*fParticlesPermuted)->Parton(1)->Px();
  blep_meas_py     = (*fParticlesPermuted)->Parton(1)->Py();
  blep_meas_pz     = (*fParticlesPermuted)->Parton(1)->Pz();
  blep_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(1)->M(), fPhysicsConstants->MassBottom(), blep_meas_px, blep_meas_py, blep_meas_pz, blep_meas_e);
  blep_meas_p      = sqrt(blep_meas_e*blep_meas_e - blep_meas_m*blep_meas_m);

  lq1_meas_e      = (*fParticlesPermuted)->Parton(2)->E();
  lq1_meas_deteta = (*fParticlesPermuted)->DetEta(2, KLFitter::Particles::kParton);
  lq1_meas_px     = (*fParticlesPermuted)->Parton(2)->Px();
  lq1_meas_py     = (*fParticlesPermuted)->Parton(2)->Py();
  lq1_meas_pz     = (*fParticlesPermuted)->Parton(2)->Pz();
  lq1_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(2)->M(), 0., lq1_meas_px, lq1_meas_py, lq1_meas_pz, lq1_meas_e);
  lq1_meas_p      = sqrt(lq1_meas_e*lq1_meas_e - lq1_meas_m*lq1_meas_m);

  lq2_meas_e      = (*fParticlesPermuted)->Parton(3)->E();
  lq2_meas_deteta = (*fParticlesPermuted)->DetEta(3, KLFitter::Particles::kParton);
  lq2_meas_px     = (*fParticlesPermuted)->Parton(3)->Px();
  lq2_meas_py     = (*fParticlesPermuted)->Parton(3)->Py();
  lq2_meas_pz     = (*fParticlesPermuted)->Parton(3)->Pz();
  lq2_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(3)->M(), 0., lq2_meas_px, lq2_meas_py, lq2_meas_pz, lq2_meas_e);
  lq2_meas_p      = sqrt(lq2_meas_e*lq2_meas_e - lq2_meas_m*lq2_meas_m);

  //++++++++++++++++//
  BHiggs1_meas_e      = (*fParticlesPermuted)->Parton(4)->E();
  BHiggs1_meas_deteta = (*fParticlesPermuted)->DetEta(4, KLFitter::Particles::kParton);
  BHiggs1_meas_px     = (*fParticlesPermuted)->Parton(4)->Px();
  BHiggs1_meas_py     = (*fParticlesPermuted)->Parton(4)->Py();
  BHiggs1_meas_pz     = (*fParticlesPermuted)->Parton(4)->Pz();
  BHiggs1_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(4)->M(), fPhysicsConstants->MassBottom(), BHiggs1_meas_px, BHiggs1_meas_py, BHiggs1_meas_pz, BHiggs1_meas_e);
  BHiggs1_meas_p      = sqrt(BHiggs1_meas_e*BHiggs1_meas_e - BHiggs1_meas_m*BHiggs1_meas_m);
  
  BHiggs2_meas_e      = (*fParticlesPermuted)->Parton(5)->E();
  BHiggs2_meas_deteta = (*fParticlesPermuted)->DetEta(5, KLFitter::Particles::kParton);
  BHiggs2_meas_px     = (*fParticlesPermuted)->Parton(5)->Px();
  BHiggs2_meas_py     = (*fParticlesPermuted)->Parton(5)->Py();
  BHiggs2_meas_pz     = (*fParticlesPermuted)->Parton(5)->Pz();
  BHiggs2_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(5)->M(), fPhysicsConstants->MassBottom(), BHiggs2_meas_px, BHiggs2_meas_py, BHiggs2_meas_pz, BHiggs2_meas_e);
  BHiggs2_meas_p      = sqrt(BHiggs2_meas_e*BHiggs2_meas_e - BHiggs2_meas_m*BHiggs2_meas_m);
  //++++++++++++++++//
  
  TLorentzVector * lepton(0);
  if (fTypeLepton == kElectron) {
    lepton = (*fParticlesPermuted)->Electron(0);
    lep_meas_deteta = (*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kElectron);
  }
  else {
    lepton = (*fParticlesPermuted)->Muon(0);
    lep_meas_deteta = (*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kMuon);
  }
  lep_meas_e        = lepton->E();
  lep_meas_sintheta = sin(lepton->Theta());
  lep_meas_pt       = lepton->Pt();
  lep_meas_px       = lepton->Px();
  lep_meas_py       = lepton->Py();
  lep_meas_pz       = lepton->Pz();

  // no error
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTHLeptonJets::SaveResolutionFunctions() {

  fResEnergyBhad = (*fDetector)->ResEnergyBJet(bhad_meas_deteta);
  fResEnergyBlep = (*fDetector)->ResEnergyBJet(blep_meas_deteta);
  fResEnergyLQ1  = (*fDetector)->ResEnergyLightJet(lq1_meas_deteta);
  fResEnergyLQ2  = (*fDetector)->ResEnergyLightJet(lq2_meas_deteta);
  
  //++++++++++++++++//
  fResEnergyBHiggs1 = (*fDetector)->ResEnergyBJet(BHiggs1_meas_deteta);
  fResEnergyBHiggs2 = (*fDetector)->ResEnergyBJet(BHiggs2_meas_deteta);
  //++++++++++++++++//

  if (fTypeLepton == kElectron)
    fResLepton = (*fDetector)->ResEnergyElectron(lep_meas_deteta);
  else if (fTypeLepton == kMuon)
    fResLepton = (*fDetector)->ResEnergyMuon(lep_meas_deteta);
  
  fResMET = (*fDetector)->ResMissingET();

  // no error
  return 1;
}
// --------------------------------------------------------- 

int KLFitter::LikelihoodTTHLeptonJets::BuildModelParticles() {	
if (GetBestFitParameters().size() > 0) CalculateLorentzVectors(GetBestFitParameters());

  TLorentzVector * bhad = fParticlesModel->Parton(0);
  TLorentzVector * blep = fParticlesModel->Parton(1);
  TLorentzVector * lq1  = fParticlesModel->Parton(2);
  TLorentzVector * lq2  = fParticlesModel->Parton(3);
  
  //++++++++++++++++//
  TLorentzVector * BHiggs1  = fParticlesModel->Parton(4);
  TLorentzVector * BHiggs2  = fParticlesModel->Parton(5);
  //++++++++++++++++//

  TLorentzVector * lep(0);
  if (fTypeLepton == kElectron)
    lep  = fParticlesModel->Electron(0);
  else if (fTypeLepton == kMuon)
    lep  = fParticlesModel->Muon(0);
  TLorentzVector * nu   = fParticlesModel->Neutrino(0);
  TLorentzVector * whad  = fParticlesModel->Boson(0);
  TLorentzVector * wlep  = fParticlesModel->Boson(1);

  //++++++++++++++++//
  TLorentzVector * thad  = fParticlesModel->Parton(6);
  TLorentzVector * tlep  = fParticlesModel->Parton(7);
  TLorentzVector * Higgs = fParticlesModel->Boson(2);
  //++++++++++++++++//

  bhad->SetPxPyPzE(bhad_fit_px, bhad_fit_py, bhad_fit_pz, bhad_fit_e);
  blep->SetPxPyPzE(blep_fit_px, blep_fit_py, blep_fit_pz, blep_fit_e);
  lq1 ->SetPxPyPzE(lq1_fit_px,  lq1_fit_py,  lq1_fit_pz,  lq1_fit_e);
  lq2 ->SetPxPyPzE(lq2_fit_px,  lq2_fit_py,  lq2_fit_pz,  lq2_fit_e);

  //++++++++++++++++//
  BHiggs1->SetPxPyPzE(BHiggs1_fit_px, BHiggs1_fit_py, BHiggs1_fit_pz, BHiggs1_fit_e);
  BHiggs2->SetPxPyPzE(BHiggs2_fit_px, BHiggs2_fit_py, BHiggs2_fit_pz, BHiggs2_fit_e);
  //++++++++++++++++//

  lep ->SetPxPyPzE(lep_fit_px,  lep_fit_py,  lep_fit_pz,  lep_fit_e);
  nu  ->SetPxPyPzE(nu_fit_px,   nu_fit_py,   nu_fit_pz,   nu_fit_e);

  (*whad) = (*lq1)  + (*lq2);
  (*wlep) = (*lep)  + (*nu);
  (*thad) = (*whad) + (*bhad);
  (*tlep) = (*wlep) + (*blep);
  
  //++++++++++++++++//
  (*Higgs) = (*BHiggs1) + (*BHiggs2); 
  //++++++++++++++++//

  // no error
  return 1;
}
// --------------------------------------------------------- 

// --------------------------------------------------------- 
std::vector<double> KLFitter::LikelihoodTTHLeptonJets::LogLikelihoodComponents(std::vector<double> parameters)
{
std::vector<double> vecci;

  // calculate 4-vectors 
  CalculateLorentzVectors(parameters); 


  // temporary flag for a safe use of the transfer functions
  bool TFgoodTmp(true);

  // jet energy resolution terms 
  vecci.push_back(log( fResEnergyBhad->p(bhad_fit_e, bhad_meas_e, TFgoodTmp) )); //comp0
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log( fResEnergyBlep->p(blep_fit_e, blep_meas_e, TFgoodTmp) )); //comp1
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log( fResEnergyLQ1->p(lq1_fit_e, lq1_meas_e, TFgoodTmp) ));  //comp2
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log( fResEnergyLQ2->p(lq2_fit_e, lq2_meas_e, TFgoodTmp) ));  //comp3
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log( fResEnergyBHiggs1->p(BHiggs1_fit_e, BHiggs1_meas_e, TFgoodTmp) ));  //comp4
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log( fResEnergyBHiggs2->p(BHiggs2_fit_e, BHiggs2_meas_e, TFgoodTmp) ));  //comp5
  if (!TFgoodTmp) fTFgood = false;

  // lepton energy resolution terms 
  if (fTypeLepton == kElectron){
  vecci.push_back(log( fResLepton->p(lep_fit_e, lep_meas_e, TFgoodTmp) )); //comp6
  }
  else if (fTypeLepton == kMuon)
  vecci.push_back(log( fResLepton->p(lep_fit_e* lep_meas_sintheta, lep_meas_pt, TFgoodTmp) )); //comp6
  if (!TFgoodTmp) fTFgood = false;

  // neutrino px and py
  vecci.push_back(log( fResMET->p(nu_fit_px, ETmiss_x, TFgoodTmp, SumET) )); //comp7
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log( fResMET->p(nu_fit_py, ETmiss_y, TFgoodTmp, SumET) )); //comp8
  if (!TFgoodTmp) fTFgood = false;

  // physics constants
  double massW = fPhysicsConstants->MassW();
  double gammaW = fPhysicsConstants->GammaW();
  double gammaTop = fPhysicsConstants->GammaTop();
  double gammaHiggs = fPhysicsConstants->GammaHiggs();
  
  // Breit-Wigner of hadronically decaying W-boson
  vecci.push_back(BCMath::LogBreitWignerRel(whad_fit_m, massW, gammaW)); //comp9

  // Breit-Wigner of leptonically decaying W-boson
  vecci.push_back(BCMath::LogBreitWignerRel(wlep_fit_m, massW, gammaW)); //comp10

  // Breit-Wigner of hadronically decaying top quark
  vecci.push_back(BCMath::LogBreitWignerRel(thad_fit_m, parameters[parTopM], gammaTop)); //comp11
        
  // Breit-Wigner of leptonically decaying top quark
  vecci.push_back(BCMath::LogBreitWignerRel(tlep_fit_m, parameters[parTopM], gammaTop)); //comp12

  // Breit-Wigner of Higgs decaying into 2 b-quark
  if (fFlagHiggsMassFixed)  vecci.push_back(BCMath::LogBreitWignerRel(Higgs_fit_m, parameters[parHiggsM], gammaHiggs)); //comp13

  // return log of likelihood 
  return vecci; 
}
