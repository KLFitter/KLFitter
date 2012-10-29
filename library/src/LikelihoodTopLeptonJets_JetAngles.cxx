#include "LikelihoodTopLeptonJets_JetAngles.h" 
#include "ResolutionBase.h"
#include "Particles.h"
#include "Permutations.h"
#include "PhysicsConstants.h"
#include "DetectorBase.h"

#include <iostream> 
#include <algorithm> 

#include <BAT/BCMath.h> 

// --------------------------------------------------------- 
KLFitter::LikelihoodTopLeptonJets_JetAngles::LikelihoodTopLeptonJets_JetAngles() : KLFitter::LikelihoodBase::LikelihoodBase()
                                                             , fFlagTopMassFixed(false)
                                                             , fFlagUseJetMass(false)
                                                             , ETmiss_x(0.)
                                                             , ETmiss_y(0.)
                                                             , fTypeLepton(kElectron)
                                                             , fPi(3.14159265358979312)
                                                             , fTwoPi(6.28318530717958623)
                                                             , fTFgood(true)
{
  // define model particles 
  DefineModelParticles(); 

  // define parameters 
  DefineParameters(); 
}

// --------------------------------------------------------- 
KLFitter::LikelihoodTopLeptonJets_JetAngles::~LikelihoodTopLeptonJets_JetAngles()
{
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopLeptonJets_JetAngles::SetET_miss_XY(double etx, double ety)
{
  // set missing ET x and y component
  ETmiss_x = etx;
  ETmiss_y = ety;

  // no error
  return 1;
}

// --------------------------------------------------------- 
void KLFitter::LikelihoodTopLeptonJets_JetAngles::SetLeptonType(LeptonType leptontype)
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
void KLFitter::LikelihoodTopLeptonJets_JetAngles::SetLeptonType(int leptontype)
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
int KLFitter::LikelihoodTopLeptonJets_JetAngles::DefineModelParticles()
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
  //free memory
  delete dummy;
 
  // no error 
  return 1;
}

// --------------------------------------------------------- 
void KLFitter::LikelihoodTopLeptonJets_JetAngles::DefineParameters()
{
  // add parameters of model
  AddParameter("energy hadronic b",       fPhysicsConstants->MassBottom(), 1000.0); // par 0
  AddParameter("energy leptonic b",       fPhysicsConstants->MassBottom(), 1000.0); // par 1
  AddParameter("energy light quark 1",    0.0, 1000.0);                               // par 2
  AddParameter("energy light quark 2",    0.0, 1000.0);                               // par 3
  AddParameter("energy lepton",           0.0, 1000.0);                               // par 4
  AddParameter("p_z neutrino",        -1000.0, 1000.0);                               // par 5
  AddParameter("missPx",              -1000.0, 1000.0);                               // par 6
  AddParameter("missPy",              -1000.0, 1000.0);                               // par 7
  AddParameter("eta hadronic b",       -2.5, 2.5);                                    // par 8
  AddParameter("eta leptonic b",       -2.5, 2.5);                                    // par 9
  AddParameter("eta light quark 1",    -2.5, 2.5);                                    // par 10
  AddParameter("eta light quark 2",    -2.5, 2.5);                                    // par 11
  AddParameter("phi hadronic b",       -fPi, fPi);                                    // par 12
  AddParameter("phi leptonic b",       -fPi, fPi);                                    // par 13
  AddParameter("phi light quark 1",    -fPi, fPi);                                    // par 14
  AddParameter("phi light quark 2",    -fPi, fPi);                                    // par 15
  AddParameter("top mass",             100.0, 1000.0);                                // par 16
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopLeptonJets_JetAngles::CalculateLorentzVectors(std::vector <double> const& parameters)
{
  // variables
  double E, m, px, py, pz, pabs, theta, scale; 
  TLorentzVector* vect; 

  // short cut to model particles for CPU-time reasons
  TLorentzVector* parton0 = fParticlesModel->Parton(0); 
  TLorentzVector* parton1 = fParticlesModel->Parton(1); 
  TLorentzVector* parton2 = fParticlesModel->Parton(2); 
  TLorentzVector* parton3 = fParticlesModel->Parton(3); 
  TLorentzVector* neutrino = fParticlesModel->Neutrino(0); 

  // hadronic b quark 
  E = parameters[0]; 
  if (fFlagUseJetMass)
    m = std::max(0., (*fParticlesPermuted)->Parton(0)->M()); 
  else
    m = fPhysicsConstants->MassBottom(); 
  pabs = sqrt(E*E - m*m);
  theta = atan(exp(-parameters[8]))*2;
  px = pabs* cos(parameters[12])* sin(theta);
  py = pabs* sin(parameters[12])* sin(theta);
  pz = pabs* cos(theta);
  parton0->SetPxPyPzE(px, py, pz, E);

  // leptonic b quark 
  E = parameters[1]; 
  if (fFlagUseJetMass)
    m = std::max(0., (*fParticlesPermuted)->Parton(1)->M());
  else
    m = fPhysicsConstants->MassBottom(); 
  pabs = sqrt(E*E - m*m);
  theta = atan(exp(-parameters[9]))*2;
  px = pabs* cos(parameters[13])* sin(theta);
  py = pabs* sin(parameters[13])* sin(theta);
  pz = pabs* cos(theta);
  parton1->SetPxPyPzE(px, py, pz, E);

  // light quark 1 
  E = parameters[2]; 
  if (fFlagUseJetMass)
    m = std::max(0., (*fParticlesPermuted)->Parton(2)->M());
  else
    m = 0; 
  pabs = sqrt(E*E - m*m);
  theta = atan(exp(-parameters[10]))*2;
  px = pabs* cos(parameters[14])* sin(theta);
  py = pabs* sin(parameters[14])* sin(theta);
  pz = pabs* cos(theta);
  parton2->SetPxPyPzE(px, py, pz, E);

  // light quark 2 
  E = parameters[3]; 
  if (fFlagUseJetMass)
    m = std::max(0., (*fParticlesPermuted)->Parton(3)->M());
  else
    m = 0; 
  pabs = sqrt(E*E - m*m);
  theta = atan(exp(-parameters[11]))*2;
  px = pabs* cos(parameters[15])* sin(theta);
  py = pabs* sin(parameters[15])* sin(theta);
  pz = pabs* cos(theta);
  parton3->SetPxPyPzE(px, py, pz, E);

  // lepton
  if (fTypeLepton == kElectron)
    {
      vect = (*fParticlesPermuted)->Electron(0); 
      E = parameters[4]; 
      px = vect->Px(); 
      py = vect->Py(); 
      pz = vect->Pz(); 
      pabs = sqrt(px*px + py*py + pz*pz); 
      scale = E / vect->E(); 
      fParticlesModel->Electron(0)->SetPxPyPzE(scale* px, scale* py, scale* pz, E); 
    }
  else if (fTypeLepton == kMuon)
    {
      vect = (*fParticlesPermuted)->Muon(0); 
      E = parameters[4]; 
      px = vect->Px(); 
      py = vect->Py(); 
      pz = vect->Pz(); 
      pabs = sqrt(px*px + py*py + pz*pz); 
      scale = E / vect->E(); 
      fParticlesModel->Muon(0)->SetPxPyPzE(scale* px, scale* py, scale* pz, E); 
    }

  // neutrino 
  px = parameters[6]; // yes: (px, py, pz) = (par6, par7, par5) - this order (!)
  py = parameters[7];
  pz = parameters[5];
  E = sqrt(px*px + py*py + pz*pz);
  neutrino->SetPxPyPzE(px, py, pz, E); 

  // composite particles 

  // hadronic W 
  *(fParticlesModel->Boson(0)) = *(parton2) + *(parton3); 
        
  // leptonic W 
  if (fTypeLepton == kElectron)
    *(fParticlesModel->Boson(1)) = *(fParticlesModel->Electron(0)) + *(neutrino); 
  else if (fTypeLepton == kMuon)
    *(fParticlesModel->Boson(1)) = *(fParticlesModel->Muon(0)) + *(neutrino); 

  // hadronic top 
  *(fParticlesModel->Parton(4)) = *(fParticlesModel->Boson(0)) + *(parton0); 

  // leptonic top 
  *(fParticlesModel->Parton(5)) = *(fParticlesModel->Boson(1)) + *(parton1); 

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopLeptonJets_JetAngles::Initialize()
{
  // error code 
  int err = 1; 

  // adjust parameter ranges 
  err *= AdjustParameterRanges(); 

  // set initial values
  // (only for Markov chains - initial parameters for other minimisation methods are set in Fitter.cxx)
  SetInitialParameters( GetInitialParameters() );       

  // return error code 
  return err; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopLeptonJets_JetAngles::RemoveInvariantParticlePermutations()
{
  // error code 
  int err = 1; 

  // remove the permutation from the second and the third jet
  KLFitter::Particles::ParticleType ptype = KLFitter::Particles::kParton;
  std::vector<int> indexVector_Jets;
  indexVector_Jets.push_back(2);
  indexVector_Jets.push_back(3);
  err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Jets); 
        
  //remove invariant jet permutations of notevent jets
  KLFitter::Particles* particles = (*fPermutations)->Particles();
  indexVector_Jets.clear();
  for (int iPartons = 4; iPartons < particles->NPartons(); iPartons++)
    indexVector_Jets.push_back(iPartons);
  err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Jets);

  // remove the permutation from the other lepton
  if (fTypeLepton == kElectron)
    {
      ptype = KLFitter::Particles::kMuon;
      std::vector<int> indexVector_Muons;
      for (int iMuon = 0; iMuon < particles->NMuons(); iMuon++)
        indexVector_Muons.push_back(iMuon);
      err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Muons); 
    }
  if (fTypeLepton == kMuon)
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
int KLFitter::LikelihoodTopLeptonJets_JetAngles::RemoveForbiddenParticlePermutations()
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
int KLFitter::LikelihoodTopLeptonJets_JetAngles::AdjustParameterRanges()
{
  // adjust limits 
  double nsigmas_jet = 7.0; 
  double nsigmas_lepton = 2.0; 

  // parameter 0: energy of hadronic b quark 
  double E = (*fParticlesPermuted)->Parton(0)->E(); 
  double m = fPhysicsConstants->MassBottom(); 
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(0)->M()); 
  double Emin = std::max(m, E - nsigmas_jet* sqrt(E)); 
  double Emax  = E + nsigmas_jet* sqrt(E); 
  SetParameterRange(0, Emin, Emax); 

  // parameter 1: energy of leptonic b quark 
  E = (*fParticlesPermuted)->Parton(1)->E(); 
  m = fPhysicsConstants->MassBottom(); 
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(1)->M()); 
  Emin = std::max(m, E - nsigmas_jet* sqrt(E)); 
  Emax  = E + nsigmas_jet* sqrt(E); 
  SetParameterRange(1, Emin, Emax); 

  // parameter 2: energy of light quark 1 
  E = (*fParticlesPermuted)->Parton(2)->E(); 
  m = 0.001;
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(2)->M()); 
  Emin = std::max(m, E - nsigmas_jet* sqrt(E)); 
  Emax  = E + nsigmas_jet* sqrt(E); 
  SetParameterRange(2, Emin, Emax); 

  // parameter 3: energy of light quark2 
  E = (*fParticlesPermuted)->Parton(3)->E(); 
  m = 0.001;
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(3)->M()); 
  Emin = std::max(m, E - nsigmas_jet* sqrt(E)); 
  Emax  = E + nsigmas_jet* sqrt(E); 
  SetParameterRange(3, Emin, Emax); 

  // parameter 2: energy of lepton
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
  SetParameterRange(4, Emin, Emax); 

  // note: this is hard-coded in the momement 

  // missing ET 
  SetParameterRange(6, ETmiss_x-100.0, ETmiss_x+100);
  SetParameterRange(7, ETmiss_y-100.0, ETmiss_y+100);

  // eta
  double eta = (*fParticlesPermuted)->Parton(0)->Eta(); 
  double etamin = std::max(-2.5, eta - 0.2); 
  double etamax = std::min(2.5, eta + 0.2); 
  //SetParameterRange(8, etamin, etamax);
  SetParameterRange(8, eta, eta); // Fix temporarily eta & phi until proper TFs are available

  eta = (*fParticlesPermuted)->Parton(1)->Eta(); 
  etamin = std::max(-2.5, eta - 0.2); 
  etamax = std::min(2.5, eta + 0.2); 
  //SetParameterRange(9, etamin, etamax);
  SetParameterRange(9, eta, eta); // Fix temporarily eta & phi until proper TFs are available 

  eta = (*fParticlesPermuted)->Parton(2)->Eta(); 
  etamin = std::max(-2.5, eta - 0.2); 
  etamax = std::min(2.5, eta + 0.2); 
  //SetParameterRange(10, etamin, etamax);
  SetParameterRange(10, eta, eta); // Fix temporarily eta & phi until proper TFs are available 

  eta = (*fParticlesPermuted)->Parton(3)->Eta(); 
  etamin = std::max(-2.5, eta - 0.2); 
  etamax = std::min(2.5, eta + 0.2); 
  //SetParameterRange(11, etamin, etamax);
  SetParameterRange(11, eta, eta); // Fix temporarily eta & phi until proper TFs are available 

  // phi 
  double phi = (*fParticlesPermuted)->Parton(0)->Phi(); 
  double phimin = phi - 0.1;
  double phimax = phi + 0.1;
  //SetParameterRange(12, phimin, phimax);
  SetParameterRange(12, phi, phi); // Fix temporarily eta & phi until proper TFs are available 

  phi = (*fParticlesPermuted)->Parton(1)->Phi(); 
  phimin = phi - 0.1;
  phimax = phi + 0.1;
  //SetParameterRange(13, phimin, phimax); 
  SetParameterRange(13, phi, phi); // Fix temporarily eta & phi until proper TFs are available

  phi = (*fParticlesPermuted)->Parton(2)->Phi(); 
  phimin = phi - 0.1;
  phimax = phi + 0.1;
  //SetParameterRange(14, phimin, phimax);
  SetParameterRange(14, phi, phi); // Fix temporarily eta & phi until proper TFs are available 

  phi = (*fParticlesPermuted)->Parton(3)->Phi(); 
  phimin = phi - 0.1;
  phimax = phi + 0.1;
  //SetParameterRange(15, phimin, phimax);
  SetParameterRange(15, phi, phi); // Fix temporarily eta & phi until proper TFs are available
 

  // top mass 
  if (fFlagTopMassFixed)
    SetParameterRange(16, fPhysicsConstants->MassTop(), fPhysicsConstants->MassTop()); 

  // no error 
  return 1; 
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodTopLeptonJets_JetAngles::LogLikelihood(const std::vector<double> & parameters)
{
  //    // check if W mass is within range
  //    TLorentzVector Whad; 
  //    double px, py, pz, E; 
  //    px = (*fParticlesPermuted)->Parton(2)->Px() + (*fParticlesPermuted)->Parton(3)->Px(); 
  //    py = (*fParticlesPermuted)->Parton(2)->Py() + (*fParticlesPermuted)->Parton(3)->Py(); 
  //    pz = (*fParticlesPermuted)->Parton(2)->Pz() + (*fParticlesPermuted)->Parton(3)->Pz(); 
  //    E  = (*fParticlesPermuted)->Parton(2)->E()  + (*fParticlesPermuted)->Parton(3)->E(); 
  //    Whad.SetPxPyPzE(px, py, pz, E); 
  //
  //    if ( fabs(Whad.M() - fPhysicsConstants->MassW()) > 40.0 )
  //            return -500; 

  // short cut for particles for CPU-time reasons
  TLorentzVector* parton0 = (*fParticlesPermuted)->Parton(0);
  TLorentzVector* parton1 = (*fParticlesPermuted)->Parton(1);
  TLorentzVector* parton2 = (*fParticlesPermuted)->Parton(2);
  TLorentzVector* parton3 = (*fParticlesPermuted)->Parton(3);

  //Sometimes the Error ResEnergyBJet .Eta Range exceeded appears. As fast bug fix the next two if statements are implemented
  // - need only to check one variable (the first one to be used later one) - all parameters are 'nan'
  if(isnan(parameters[8]))
    {
      SetFlagIsNan(true);
      return -999999.999;
    }

  // calculate 4-vectors 
  CalculateLorentzVectors(parameters); 

  // define log of likelihood 
  double logprob = 0.; 

  fTFgood = true;
  bool TFgoodTmp(true);

  // jet energy resolution terms 
  logprob += log( (*fDetector)->ResEnergyBJet( (*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kParton) )->p( parameters[0], parton0->E(), TFgoodTmp) ); if (!TFgoodTmp) fTFgood = false;
  logprob += log( (*fDetector)->ResEnergyBJet( (*fParticlesPermuted)->DetEta(1, KLFitter::Particles::kParton) )->p( parameters[1], parton1->E(), TFgoodTmp) );  if (!TFgoodTmp) fTFgood = false;
  logprob += log( (*fDetector)->ResEnergyLightJet( (*fParticlesPermuted)->DetEta(2, KLFitter::Particles::kParton) )->p( parameters[2], parton2->E(), TFgoodTmp) );  if (!TFgoodTmp) fTFgood = false;
  logprob += log( (*fDetector)->ResEnergyLightJet( (*fParticlesPermuted)->DetEta(3, KLFitter::Particles::kParton) )->p( parameters[3], parton3->E(), TFgoodTmp) );  if (!TFgoodTmp) fTFgood = false;

  // lepton energy resolution terms 
  if (fTypeLepton == kElectron) {
    TLorentzVector* electron = (*fParticlesPermuted)->Electron(0);
    logprob += log( (*fDetector)->ResEnergyElectron( (*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kElectron) )->p( parameters[4], electron->E(), TFgoodTmp) );  if (!TFgoodTmp) fTFgood = false;
  }
  
  else if (fTypeLepton == kMuon) {
    //    logprob += log( (*fDetector)->ResEnergyMuon( (*fParticlesPermuted)->Muon(0)->Eta() )->p( parameters[4],               (*fParticlesPermuted)->Muon(0)->E(), TFgoodTmp) ); if (!TFgoodTmp) fTFgood = false;
    TLorentzVector* muon = (*fParticlesPermuted)->Muon(0);
    double pt_fit = parameters[4]* sin(muon->Theta());
    double pt_reco = muon->Pt();
    logprob += log( (*fDetector)->ResEnergyMuon( (*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kMuon) )->p( pt_fit, pt_reco, TFgoodTmp )); if (!TFgoodTmp) fTFgood = false;
  }

  //    ETmiss_x and ETmiss_y
  logprob += log( (*fDetector)->ResMissingET()->p( parameters[6], ETmiss_x, TFgoodTmp) );  if (!TFgoodTmp) fTFgood = false;
  logprob += log( (*fDetector)->ResMissingET()->p( parameters[7], ETmiss_y, TFgoodTmp) );  if (!TFgoodTmp) fTFgood = false;

  // eta resolution 
  logprob += log( (*fDetector)->ResEtaBJet((*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kParton))->p(parameters[8], parton0->Eta(), TFgoodTmp) );  if (!TFgoodTmp) fTFgood = false;
  logprob += log( (*fDetector)->ResEtaBJet((*fParticlesPermuted)->DetEta(1, KLFitter::Particles::kParton))->p(parameters[9], parton1->Eta(), TFgoodTmp) );  if (!TFgoodTmp) fTFgood = false;
  logprob += log( (*fDetector)->ResEtaLightJet((*fParticlesPermuted)->DetEta(2, KLFitter::Particles::kParton))->p(parameters[10], parton2->Eta(), TFgoodTmp) ); if (!TFgoodTmp) fTFgood = false;
  logprob += log( (*fDetector)->ResEtaLightJet((*fParticlesPermuted)->DetEta(3, KLFitter::Particles::kParton))->p(parameters[11], parton3->Eta(), TFgoodTmp) ); if (!TFgoodTmp) fTFgood = false;

  double temp_phi_angles[4];

  temp_phi_angles[0] = parameters[12];
  temp_phi_angles[1] = parameters[13];
  temp_phi_angles[2] = parameters[14];
  temp_phi_angles[3] = parameters[15];


  // check phi variable
  if (temp_phi_angles[0] < -fPi)
    temp_phi_angles[0] += fTwoPi; 
  if (temp_phi_angles[0] > fPi)
    temp_phi_angles[0] -= fTwoPi; 

  if (temp_phi_angles[1] < -fPi)
    temp_phi_angles[1] += fTwoPi; 
  if (temp_phi_angles[1] > fPi)
    temp_phi_angles[1] -= fTwoPi; 

  if (temp_phi_angles[2] < -fPi)
    temp_phi_angles[2] += fTwoPi; 
  if (temp_phi_angles[2] > fPi)
    temp_phi_angles[2] -= fTwoPi; 

  if (temp_phi_angles[3] < -fPi)
    temp_phi_angles[3] += fTwoPi; 
  if (temp_phi_angles[3] > fPi)
    temp_phi_angles[3] -= fTwoPi; 
        
  // phi resolution (this implementation assumes a symmetric phi resolution function
  double diff0 = fabs(temp_phi_angles[0]-parton0->Phi());
  double diff1 = fabs(temp_phi_angles[1]-parton1->Phi());
  double diff2 = fabs(temp_phi_angles[2]-parton2->Phi());
  double diff3 = fabs(temp_phi_angles[3]-parton3->Phi());

  if (diff0>fPi)
    diff0 -= fTwoPi;
  if (diff1>fPi)
    diff1 -= fTwoPi;
  if (diff2>fPi)
    diff2 -= fTwoPi;
  if (diff3>fPi)
    diff3 -= fTwoPi;

  logprob += log( (*fDetector)->ResPhiBJet((*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kParton))->p(0.0, diff0, TFgoodTmp) );  if (!TFgoodTmp) fTFgood = false;
  logprob += log( (*fDetector)->ResPhiBJet((*fParticlesPermuted)->DetEta(1, KLFitter::Particles::kParton))->p(0.0, diff1, TFgoodTmp) );  if (!TFgoodTmp) fTFgood = false;
  logprob += log( (*fDetector)->ResPhiLightJet((*fParticlesPermuted)->DetEta(2, KLFitter::Particles::kParton))->p(0.0, diff2, TFgoodTmp) );  if (!TFgoodTmp) fTFgood = false;
  logprob += log( (*fDetector)->ResPhiLightJet((*fParticlesPermuted)->DetEta(3, KLFitter::Particles::kParton))->p(0.0, diff3, TFgoodTmp) );  if (!TFgoodTmp) fTFgood = false;

  // Breit-Wigner of hadronically decaying W-boson
  logprob += BCMath::LogBreitWignerRel( fParticlesModel->Boson(0)->M(), fPhysicsConstants->MassW(), fPhysicsConstants->GammaW()); 

  // Breit-Wigner of leptonically decaying W-boson
  logprob += BCMath::LogBreitWignerRel( fParticlesModel->Boson(1)->M(), fPhysicsConstants->MassW(), fPhysicsConstants->GammaW()); 

  // note: top mass width should be made DEPENDENT on the top mass at a certain point
  //    fPhysicsConstants->SetMassTop(parameters[16]);

  // Breit-Wigner of hadronically decaying top quark
  logprob += BCMath::LogBreitWignerRel( fParticlesModel->Parton(4)->M(), parameters[16], fPhysicsConstants->GammaTop()); 
        
  // Breit-Wigner of leptonically decaying top quark
  logprob += BCMath::LogBreitWignerRel( fParticlesModel->Parton(5)->M(), parameters[16], fPhysicsConstants->GammaTop()); 

  // return log of likelihood 
  return logprob; 
}

// --------------------------------------------------------- 
std::vector<double> KLFitter::LikelihoodTopLeptonJets_JetAngles::GetInitialParameters()
{
  std::vector<double> values(GetInitialParametersWoNeutrinoPz());

  // check second neutrino solution
  std::vector<double> neutrino_pz_solutions = GetNeutrinoPzSolutions();
  if (int(neutrino_pz_solutions.size()) == 1)
    values[5] = neutrino_pz_solutions[0]; 
  else if(int(neutrino_pz_solutions.size()) == 2)
    {
      double sol1, sol2; 
      values[5] = neutrino_pz_solutions[0]; 
      sol1 = LogLikelihood(values); 
      values[5] = neutrino_pz_solutions[1]; 
      sol2 = LogLikelihood(values); 

      if (sol1 > sol2)
        values[5] = neutrino_pz_solutions[0]; 
    }

  return values;
}

// --------------------------------------------------------- 
std::vector<double> KLFitter::LikelihoodTopLeptonJets_JetAngles::GetInitialParametersWoNeutrinoPz()
{
  std::vector<double> values;

  // energies of the quarks
  values.push_back( (*fParticlesPermuted)->Parton(0)->E() );
  values.push_back( (*fParticlesPermuted)->Parton(1)->E() );
  values.push_back( (*fParticlesPermuted)->Parton(2)->E() );
  values.push_back( (*fParticlesPermuted)->Parton(3)->E() );

  // energy of the lepton
  if (fTypeLepton == kElectron)
    values.push_back( (*fParticlesPermuted)->Electron(0)->E() );
  else if (fTypeLepton == kMuon)
    values.push_back( (*fParticlesPermuted)->Muon(0)->E() );

  // pz of the neutrino
  values.push_back(0.0); 

  // missing px and py
  values.push_back(ETmiss_x);
  values.push_back(ETmiss_y);

  // eta and phi values of the jets
  values.push_back( (*fParticlesPermuted)->Parton(0)->Eta() );
  values.push_back( (*fParticlesPermuted)->Parton(1)->Eta() );
  values.push_back( (*fParticlesPermuted)->Parton(2)->Eta() );
  values.push_back( (*fParticlesPermuted)->Parton(3)->Eta() );
  values.push_back( (*fParticlesPermuted)->Parton(0)->Phi() );
  values.push_back( (*fParticlesPermuted)->Parton(1)->Phi() );
  values.push_back( (*fParticlesPermuted)->Parton(2)->Phi() );
  values.push_back( (*fParticlesPermuted)->Parton(3)->Phi() );

  // top mass
  double mtop = ( *(*fParticlesPermuted)->Parton(0) + *(*fParticlesPermuted)->Parton(2) + *(*fParticlesPermuted)->Parton(3) ).M(); 
  if (mtop < GetParameter(16)->GetLowerLimit())
    mtop = GetParameter(16)->GetLowerLimit(); 
  else if (mtop > GetParameter(16)->GetUpperLimit())
    mtop = GetParameter(16)->GetUpperLimit(); 
  values.push_back( mtop);

  // return the vector
  return values;
}

// --------------------------------------------------------- 
std::vector<double> KLFitter::LikelihoodTopLeptonJets_JetAngles::GetNeutrinoPzSolutions() {
  return CalculateNeutrinoPzSolutions();
}

// --------------------------------------------------------- 
std::vector<double> KLFitter::LikelihoodTopLeptonJets_JetAngles::CalculateNeutrinoPzSolutions(TLorentzVector* additionalParticle)
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
bool KLFitter::LikelihoodTopLeptonJets_JetAngles::NoTFProblem(std::vector<double> parameters) {
  fTFgood = true;
  this->LogLikelihood(parameters);
  return fTFgood;
}
// --------------------------------------------------------- 
