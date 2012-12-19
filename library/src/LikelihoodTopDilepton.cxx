#include "LikelihoodTopDilepton.h" 
#include "ResolutionBase.h"
#include "Particles.h"
#include "Permutations.h"
#include "PhysicsConstants.h"
#include "DetectorBase.h"
#include <assert.h>
#include <iostream> 
#include <algorithm> 
#include <cmath>

#include "BAT/BCModel.h"
#include "BAT/BCH1D.h"
#include <BAT/BCMath.h> 

// --------------------------------------------------------- 
KLFitter::LikelihoodTopDilepton::LikelihoodTopDilepton() : KLFitter::LikelihoodBase::LikelihoodBase()
                                                             , fFlagTopMassFixed(false)
                                                             , fFlagUseJetMass(false)
                                                             , ETmiss_x(0.)
                                                             , ETmiss_y(0.)
                                                             , SumET(0.)
                                                             , fTypeLepton_1(kElectron)
                                                             , fTypeLepton_2(kElectron)
							     , nueta_params(0.)
							     , doSumloglik(false)
                                                             , fTFgood(true)
							     , hist_mttbar(new TH1D())
							     , hist_costheta(new TH1D())
							     , fHistMttbar(new BCH1D())
                                                             , fHistCosTheta(new BCH1D())
{  


  // define model particles 
  this->DefineModelParticles();  

  // define parameters 
  this->DefineParameters(); 

  // define prior
  this->DefinePrior();

}

// --------------------------------------------------------- 
KLFitter::LikelihoodTopDilepton::~LikelihoodTopDilepton()
{
  delete fHistMttbar;
  delete fHistCosTheta;
 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopDilepton::SetET_miss_XY_SumET(double etx, double ety, double sumet)
{
  // set missing ET x and y component and the SumET
  ETmiss_x = etx;
  ETmiss_y = ety;
  SumET = sumet;

  // no error
  return 1;
}

// --------------------------------------------------------- 
void KLFitter::LikelihoodTopDilepton::SetLeptonType(LeptonType leptontype_1, LeptonType leptontype_2)
{
  if ((leptontype_1 != kElectron && leptontype_1 != kMuon) || (leptontype_2 != kElectron && leptontype_2 != kMuon))
    {
      std::cout << "KLFitter::SetLeptonTyp(). Warning: lepton type not defined. Set electron-electron as lepton type." << std::endl; 
      fTypeLepton_1 = kElectron;
      fTypeLepton_2 = kElectron;
    }
  else {
    fTypeLepton_1 = leptontype_1; 
    fTypeLepton_2 = leptontype_2; 
  }

  // define model particles 
  DefineModelParticles(); 

  // define histograms
  DefineHistograms();
}

// --------------------------------------------------------- 
void KLFitter::LikelihoodTopDilepton::SetLeptonType(int leptontype_1, int leptontype_2 )
{
  if ((leptontype_1 != 1 && leptontype_1 != 2) || (leptontype_2 != 1 && leptontype_2 != 2))
    {
      std::cout << "KLFitter::SetLeptonTyp(). Warning: lepton type not defined. Set electron as lepton type." << std::endl; 
      leptontype_1 = 1;
      leptontype_2 = 1;
    }

  if (leptontype_1 == 1 && leptontype_2 == 1)
    SetLeptonType(kElectron, kElectron);
  else if (leptontype_1 == 1 && leptontype_2 == 2)
    SetLeptonType(kElectron, kMuon);
  else if (leptontype_1 == 2 && leptontype_2 == 2)
    SetLeptonType(kMuon, kMuon);
  else
    std::cout << "LikelihoodTopDilepton::SetLeptonType Error: NOT POSSIBLE (2,1) conf! <---" << std::endl;
}



// --------------------------------------------------------- 
int KLFitter::LikelihoodTopDilepton::DefineModelParticles()
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
                               "b quark 1",                  // name 
                               0,                            // index of corresponding particle 
                               KLFitter::Particles::kB);     // b jet (truth)

  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kParton, 
                               "b quark 2",
                               1,                            // index of corresponding particle 
                               KLFitter::Particles::kB);     // b jet (truth)

  if (fTypeLepton_1 == kElectron && fTypeLepton_2 == kMuon) {
    fParticlesModel->AddParticle(dummy,
				 KLFitter::Particles::kElectron,
				 "electron"); 
    
    fParticlesModel->AddParticle(dummy,
				 KLFitter::Particles::kMuon,
				 "muon"); 
  }
  else if (fTypeLepton_1 == kElectron && fTypeLepton_2 == kElectron) {
    fParticlesModel->AddParticle(dummy,
				 KLFitter::Particles::kElectron,
				 "electron 1"); 
    fParticlesModel->AddParticle(dummy,
				 KLFitter::Particles::kElectron,
				 "electron 2"); 
  }
  else if (fTypeLepton_1 == kMuon && fTypeLepton_2 == kMuon) {
    fParticlesModel->AddParticle(dummy,
				 KLFitter::Particles::kMuon,
				 "muon 1"); 
    fParticlesModel->AddParticle(dummy,
				 KLFitter::Particles::kMuon,
				 "muon 2"); 
  }
  
  //free memory
  delete dummy; 

  // no error 
  return 1;
}


// --------------------------------------------------------- 
void KLFitter::LikelihoodTopDilepton::DefineParameters()
{
  // add parameters of model
  
  AddParameter("top mass",              100.0, 700.0);                             // parTopM  
  AddParameter("energy b1",       fPhysicsConstants->MassBottom(), 1000.0);        // parB1E
  AddParameter("energy b2",       fPhysicsConstants->MassBottom(), 1000.0);        // parB2E
  AddParameter("energy lepton1",           0.0, 1000.0);                           // parLep1E
  AddParameter("energy lepton2",           0.0, 1000.0);                           // parLep2E
  AddParameter("antinueta",               -5.0, 5.0);                           // parAntiNuEta
  AddParameter("nueta",                   -5.0, 5.0);                           // parNuEta
  

}

// --------------------------------------------------------- 

void KLFitter::LikelihoodTopDilepton::DefinePrior()
{
  // define sharp Gaussian prior for mtop 
  if (fFlagTopMassFixed)
    SetPriorGauss(0,fPhysicsConstants->MassTop(),fPhysicsConstants->MassTopUnc());
  
}

// --------------------------------------------------------- 

void KLFitter::LikelihoodTopDilepton::DefineHistograms()
{
  const char *channel = 0;

  if (fTypeLepton_1 == kElectron && fTypeLepton_2 == kMuon) {
    channel = "emu";
  }
  else if (fTypeLepton_1 == kElectron && fTypeLepton_2 == kElectron) {
    channel = "ee";
  }
  else if (fTypeLepton_1 == kMuon && fTypeLepton_2 == kMuon) {
    channel = "mumu";
  }

  // create a new ROOT and BAT histogram for mttbar
  hist_mttbar = new TH1D(TString::Format("hist_mttbar_%s", channel), "", 100, 100., 1500.);
  fHistMttbar->SetHistogram(hist_mttbar);
  
  // create a new ROOT and BAT histogram for costheta
  hist_costheta = new TH1D(TString::Format("hist_costheta_%s", channel), "", 100, 0., 1.);
  fHistCosTheta->SetHistogram(hist_costheta);
  
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopDilepton::CalculateLorentzVectors(std::vector <double> const& parameters)  
{
  static double scale;

  // b1 quark 
  b1_fit_e = parameters[parB1E]; 
  scale = sqrt(b1_fit_e*b1_fit_e - b1_meas_m*b1_meas_m) / b1_meas_p;
  b1_fit_px = scale * b1_meas_px;
  b1_fit_py = scale * b1_meas_py;
  b1_fit_pz = scale * b1_meas_pz;
  
  // b2 quark 
  b2_fit_e = parameters[parB2E]; 
  scale = sqrt(b2_fit_e*b2_fit_e - b2_meas_m*b2_meas_m) / b2_meas_p;
  b2_fit_px = scale * b2_meas_px;
  b2_fit_py = scale * b2_meas_py;
  b2_fit_pz = scale * b2_meas_pz;

  // lepton1
  lep1_fit_e = parameters[parLep1E];
  scale = lep1_fit_e / lep1_meas_e;
  lep1_fit_px = scale * lep1_meas_px;
  lep1_fit_py = scale * lep1_meas_py;
  lep1_fit_pz = scale * lep1_meas_pz;

  // lepton2
  lep2_fit_e = parameters[parLep2E];
  scale = lep2_fit_e / lep2_meas_e;
  lep2_fit_px = scale * lep2_meas_px;
  lep2_fit_py = scale * lep2_meas_py;
  lep2_fit_pz = scale * lep2_meas_pz;

  // no error 
  return 1; 
}


// --------------------------------------------------------- 
int KLFitter::LikelihoodTopDilepton::Initialize()
{
  // error code 
  int err = 1; 

  // save the current permuted particles
  err *= SavePermutedParticles(); 

  // save the corresponding resolution functions
  err *= SaveResolutionFunctions();

  // adjust parameter ranges 
  err *= AdjustParameterRanges(); 

  // set initial values
  // only for Markov chains: SPECIFY NChains!!
  SetInitialParametersNChains( GetInitialParameters(), 5 );     

  // return error code 
  return err; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopDilepton::RemoveInvariantParticlePermutations()
{
  // error code 
  int err = 1; 

  KLFitter::Particles::ParticleType ptype = KLFitter::Particles::kParton;
  std::vector<int> indexVector_Jets;

  //remove invariant jet permutations of notevent jets
  KLFitter::Particles* particles = (*fPermutations)->Particles();
  indexVector_Jets.clear();
  for (int iPartons = 2; iPartons < particles->NPartons(); iPartons++)
    indexVector_Jets.push_back(iPartons);
  err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Jets);
  
  //remove lepton permutations within the same kind
  if (fTypeLepton_1 == kElectron && fTypeLepton_2 == kElectron) {
    ptype = KLFitter::Particles::kElectron;
    std::vector<int> indexVector_Electrons;
    for (int iElectron = 0; iElectron < particles->NElectrons(); iElectron++)
      indexVector_Electrons.push_back(iElectron);
    err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Electrons); 
  }
  else if (fTypeLepton_1 == kMuon && fTypeLepton_2 == kMuon) {
    ptype = KLFitter::Particles::kMuon;
    std::vector<int> indexVector_Muons;
    for (int iMuon = 0; iMuon < particles->NMuons(); iMuon++)
      indexVector_Muons.push_back(iMuon);
    err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Muons); 
  }

  if(doSumloglik) {
    // remove the permutation from the 2 bjets
    ptype = KLFitter::Particles::kParton;
    indexVector_Jets.clear();
    indexVector_Jets.push_back(0);
    indexVector_Jets.push_back(1);
    err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Jets); 
  }

  // return error code 
  return err; 
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopDilepton::RemoveForbiddenParticlePermutations()
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
int KLFitter::LikelihoodTopDilepton::AdjustParameterRanges()
{
  // adjust limits 
  if (fFlagTopMassFixed)
    SetParameterRange(parTopM, fPhysicsConstants->MassTop()-3*fPhysicsConstants->MassTopUnc(), fPhysicsConstants->MassTop()+3*fPhysicsConstants->MassTopUnc()); 

  double nsigmas_jet = 7.0; 
  double nsigmas_lepton = 7.0; 
 
  // TF parameters for lepton1 and lepton2
  double par2_1=0.;
  double par3_1=0.;

  double par2_2=0.;
  double par3_2=0.;
 
  fResLepton1->Par(2, par2_1);
  fResLepton1->Par(3, par3_1);
  fResLepton2->Par(2, par2_2);
  fResLepton2->Par(3, par3_2);

  double E = (*fParticlesPermuted)->Parton(0)->E(); 
  double m = fPhysicsConstants->MassBottom(); 
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(0)->M()); 
  double Emin = std::max(m, E - nsigmas_jet* sqrt(E)); 
  double Emax  = E + nsigmas_jet* sqrt(E); 
  //std::cout << "Emin b1: " << Emin << " Emax b1: " << Emax << std::endl;
  SetParameterRange(parB1E, Emin, Emax); 

  E = (*fParticlesPermuted)->Parton(1)->E(); 
  m = fPhysicsConstants->MassBottom(); 
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(1)->M()); 
  Emin = std::max(m, E - nsigmas_jet* sqrt(E)); 
  Emax  = E + nsigmas_jet* sqrt(E); 
  //std::cout << "Emin b2: " << Emin << " Emax b2: " << Emax << std::endl;
  SetParameterRange(parB2E, Emin, Emax); 

  if (fTypeLepton_1 == kElectron && fTypeLepton_2 == kMuon) {
    E = (*fParticlesPermuted)->Electron(0)->E();
    double sigrange = nsigmas_lepton*std::max(par2_1*sqrt(E),par3_1*E);
    Emin = std::max(0.001, E - sigrange); 
    Emax  = E + sigrange;
    //std::cout << "Emin el: " << Emin << " Emax el: " << Emax << std::endl;
    SetParameterRange(parLep1E, Emin, Emax);
      
    E = (*fParticlesPermuted)->Muon(0)->E();
    double sintheta= sin((*fParticlesPermuted)->Muon(0)->Theta());
    sigrange = nsigmas_lepton*std::max(par2_2*E,par3_2*E*E*sintheta);
    Emin=std::max(0.001,E -sigrange);
    Emax=E +sigrange;
    //std::cout << "Emin mu: " << Emin << " Emax mu: " << Emax << std::endl;
    SetParameterRange(parLep2E, Emin, Emax);
  }
  else if (fTypeLepton_1 == kElectron && fTypeLepton_2 == kElectron) {
    E = (*fParticlesPermuted)->Electron(0)->E();
    double sigrange = nsigmas_lepton*std::max(par2_1*sqrt(E),par3_1*E);
    Emin = std::max(0.001, E - sigrange); 
    Emax  = E + sigrange;
    // std::cout << "Emin el 1: " << Emin << " Emax el 1: " << Emax << std::endl;
    SetParameterRange(parLep1E, Emin, Emax);
    
    E = (*fParticlesPermuted)->Electron(1)->E();
    sigrange = nsigmas_lepton*std::max(par2_2*sqrt(E),par3_2*E);
    Emin = std::max(0.001, E - sigrange); 
    Emax  = E + sigrange;
    //std::cout << "Emin el 2: " << Emin << " Emax el 2: " << Emax << std::endl;
    SetParameterRange(parLep2E, Emin, Emax);
  }
  else if (fTypeLepton_1 == kMuon && fTypeLepton_2 == kMuon) {
    E = (*fParticlesPermuted)->Muon(0)->E();
    double sintheta= sin((*fParticlesPermuted)->Muon(0)->Theta());
    double sigrange= nsigmas_lepton*std::max(par2_1*E,par3_1*E*E*sintheta);
    Emin=std::max(0.001,E -sigrange);
    Emax=E +sigrange;
    //std::cout << "Emin mu 1: " << Emin << " Emax mu 1: " << Emax << std::endl;
    SetParameterRange(parLep1E, Emin, Emax);
    
    E = (*fParticlesPermuted)->Muon(1)->E();
    sintheta= sin((*fParticlesPermuted)->Muon(1)->Theta());
    sigrange= nsigmas_lepton*std::max(par2_2*E,par3_2*E*E*sintheta);
    Emin=std::max(0.001,E -sigrange);
    Emax=E +sigrange;
    //std::cout << "Emin mu 2: " << Emin << " Emax mu 2: " << Emax << std::endl;
    SetParameterRange(parLep2E, Emin, Emax);
  }

  SetParameterRange(parAntiNuEta, -5.0, 5.0); 
  SetParameterRange(parNuEta, -5.0, 5.0); 
 
  
  // no error 
  return 1; 
}

// --------------------------------------------------------- 
double KLFitter::LikelihoodTopDilepton::LogLikelihood(const std::vector<double> & parameters)
{
  // define log of likelihood 
  double logweight(0.); 

  // calculate 4-vectors 
  CalculateLorentzVectors(parameters); 

  // temporary flag for a safe use of the transfer functions
  bool TFgoodTmp(true);


  ///////// NuWT likelihood term
  double nuwt_weight(0.);
  nuwt_weight = CalculateWeight(parameters);
  // ensure that likelihood components are not zero to avoid loglik=inf
  if(nuwt_weight == 0.){
    logweight = log(1e-99);
    return logweight;
  }
  else
    logweight += log(nuwt_weight);
  //std::cout << "NUWT  : " << logweight << std::endl;
  
  if(logweight + 10 == logweight) std::cout << "NUWT inf! : " << logweight << std::endl;
  

  //////// jet energy resolution terms 
  if( fResEnergyB1->p(b1_fit_e, b1_meas_e, TFgoodTmp) == 0. ){
    logweight = log(1e-99);
    return logweight;
  }
  else
    logweight += log( fResEnergyB1->p(b1_fit_e, b1_meas_e, TFgoodTmp) );
  if (!TFgoodTmp) fTFgood = false;
    //std::cout << "TF b1 : " << log( fResEnergyB1->p(b1_fit_e, b1_meas_e, TFgoodTmp) ) << std::endl;

  if(logweight + 10 == logweight) std::cout << "TF b1 inf! : " << log( fResEnergyB1->p(b1_fit_e, b1_meas_e, TFgoodTmp) ) << std::endl;
  
  
  if( fResEnergyB2->p(b2_fit_e, b2_meas_e, TFgoodTmp) == 0.){
    logweight = log(1e-99);
    return logweight;
  }
  else
    logweight += log( fResEnergyB2->p(b2_fit_e, b2_meas_e, TFgoodTmp) );
  if (!TFgoodTmp) fTFgood = false;
    //std::cout << "TF b2 : " << log( fResEnergyB2->p(b2_fit_e, b2_meas_e, TFgoodTmp) ) << std::endl;

  if(logweight + 10 == logweight) std::cout << "TF b2 inf! : " << log( fResEnergyB2->p(b2_fit_e, b2_meas_e, TFgoodTmp) ) << std::endl;
  
  
  //////// lepton energy resolution terms EM
  if (fTypeLepton_1 == kElectron && fTypeLepton_2 == kMuon) {
    if( fResLepton1->p(lep1_fit_e, lep1_meas_e, TFgoodTmp) == 0.){
      logweight = log(1e-99);
      return logweight;
    }
    else
      logweight += log( fResLepton1->p(lep1_fit_e, lep1_meas_e, TFgoodTmp) );
    
    if( fResLepton2->p(lep2_fit_e*lep2_meas_sintheta, lep2_meas_pt, TFgoodTmp) == 0.){
      logweight = log(1e-99);
      return logweight;
    }
    else
      logweight += log( fResLepton2->p(lep2_fit_e*lep2_meas_sintheta, lep2_meas_pt, TFgoodTmp) );
    if (!TFgoodTmp) fTFgood = false;
    
    //std::cout << "TF lep emu : " << log( fResLepton1->p(lep1_fit_e, lep1_meas_e, TFgoodTmp) ) << " and " << log( fResLepton2->p(lep2_fit_e*lep2_meas_sintheta, lep2_meas_pt, TFgoodTmp) ) << std::endl;

    if(logweight + 10 == logweight) std::cout << "TF lep emu inf! : "<< log(fResLepton1->p(lep1_fit_e, lep1_meas_e, TFgoodTmp)) <<" and "<< log(fResLepton2->p(lep2_fit_e*lep2_meas_sintheta, lep2_meas_pt, TFgoodTmp)) <<std::endl;    
  }

  //////// lepton energy resolution terms EE
  else if (fTypeLepton_1 == kElectron && fTypeLepton_2 == kElectron) {
    if( fResLepton1->p(lep1_fit_e, lep1_meas_e, TFgoodTmp) == 0.){
      logweight = log(1e-99);
      return logweight;
    }
    else
      logweight += log( fResLepton1->p(lep1_fit_e, lep1_meas_e, TFgoodTmp) );
    
    if( fResLepton2->p(lep2_fit_e, lep2_meas_e, TFgoodTmp) == 0.){
      logweight = log(1e-99);
      return logweight;
    }
    else
      logweight += log( fResLepton2->p(lep2_fit_e, lep2_meas_e, TFgoodTmp) );
    if (!TFgoodTmp) fTFgood = false;
    
    if(logweight + 10 == logweight) std::cout << "TF lep ee inf! : " << log( fResLepton1->p(lep1_fit_e, lep1_meas_e, TFgoodTmp) ) << " and " << log( fResLepton2->p(lep2_fit_e, lep2_meas_e, TFgoodTmp) ) << std::endl;
  }

  //////// lepton energy resolution terms MM
  else if (fTypeLepton_1 == kMuon && fTypeLepton_2 == kMuon) {
    if( fResLepton1->p(lep1_fit_e*lep1_meas_sintheta, lep1_meas_pt, TFgoodTmp) == 0.){
      logweight = log(1e-99);
      return logweight;
    }
    else
      logweight += log( fResLepton1->p(lep1_fit_e*lep1_meas_sintheta, lep1_meas_pt, TFgoodTmp) );
    
    if( fResLepton2->p(lep2_fit_e*lep2_meas_sintheta, lep2_meas_pt, TFgoodTmp) == 0.){
      logweight = log(1e-99);
      return logweight;
    }
    else
      logweight += log( fResLepton2->p(lep2_fit_e*lep2_meas_sintheta, lep2_meas_pt, TFgoodTmp) );
    if (!TFgoodTmp) fTFgood = false;
    
    if(logweight + 10 == logweight) std::cout << "TF lep mumu inf! : " << log(fResLepton1->p(lep1_fit_e*lep1_meas_sintheta, lep1_meas_pt, TFgoodTmp)) << " and " << log(fResLepton2->p(lep2_fit_e*lep2_meas_sintheta, lep2_meas_pt, TFgoodTmp)) << std::endl;
  }

  //////// Antineutrino eta term
  if (GaussAntiNuEta(parameters) == 0.){
    logweight = log(1e-99);
    return logweight;
  }
  else
    logweight += log( GaussAntiNuEta(parameters) );
  //std::cout << "Gauss AntiNuEta  : " << log( GaussAntiNuEta(parameters) ) << std::endl;

  if(logweight + 10 == logweight) std::cout << "Gauss AntiNuEta inf! : " << log( GaussAntiNuEta(parameters) ) << std::endl;


  //////// Neutrino eta term
  if (GaussNuEta(parameters) == 0.){
    logweight = log(1e-99);
    return logweight;
  }
  else
    logweight += log( GaussNuEta(parameters) );
  //std::cout << "Gauss NuEta : " << log( GaussNuEta(parameters) ) << std::endl;

  if(logweight + 10 == logweight) std::cout << "Gauss NuEta inf! : " << log( GaussNuEta(parameters) ) << std::endl;


  //////// Sum of invariant masses (lep,jet) term
  if (CalculateMLepJet(parameters) == 0.){
    logweight = log(1e-99);
    return logweight;
  }
  else 
    logweight += log( CalculateMLepJet(parameters) );


  if(logweight + 10 == logweight) std::cout << "CalculateMLepJet inf! : "  << log( CalculateMLepJet(parameters) ) << std::endl;

  // check if total logweight is inf
  if(logweight + 10 == logweight) std::cout << "logweight TOT = " << logweight << std::endl;

  // return log of weight
  return logweight; 
}

double KLFitter::LikelihoodTopDilepton::CalculateMLepJet(const std::vector<double> & parameters){

 double sumMinv=0.;
 double norm=1.;

 // tuning factor alpha
 double alpha=-2.;
 
 //charged leptons
 TLorentzVector l1(0.);
 TLorentzVector l2(0.);
 
 // include parLep1E, parLep2E
 l1.SetPxPyPzE(lep1_fit_px,  lep1_fit_py,  lep1_fit_pz,  lep1_fit_e);
 l2.SetPxPyPzE(lep2_fit_px,  lep2_fit_py,  lep2_fit_pz,  lep2_fit_e);
 
 // jet1 and jet2:
 TLorentzVector j1(0.);
 TLorentzVector j2(0.);
 
 // include parB1E, parB2E
 j1.SetPxPyPzE(b1_fit_px, b1_fit_py, b1_fit_pz, b1_fit_e);
 j2.SetPxPyPzE(b2_fit_px, b2_fit_py, b2_fit_pz, b2_fit_e);

 // normalized to the sum of all combinations of (lep,jet)
 if ((pow((l1+j1).M() + (l2+j2).M(),alpha) + pow((l2+j1).M() + (l1+j2).M(),alpha))!=0.)
   norm = 1/(pow((l1+j1).M() + (l2+j2).M(),alpha) + pow((l2+j1).M() + (l1+j2).M(),alpha));
 else
   std::cout << "Error LikelihoodTopDilepton::CalculateMLepJet: normalization is inf!" << std::endl;

 // ensure correctly (lepton, nu) pair according to lepton charge
 if (lep1_meas_charge == 1 && lep2_meas_charge == -1){
   //std::cout << "opt1: lep1_meas_charge: " << lep1_meas_charge << " and lep2_meas_charge: " << lep2_meas_charge  << std::endl;
   sumMinv = (l1+j1).M() + (l2+j2).M();
 }
 else if (lep1_meas_charge == -1 && lep2_meas_charge == 1){
   //std::cout << "opt2: lep1_meas_charge: " << lep1_meas_charge << " and lep2_meas_charge: " << lep2_meas_charge  << std::endl;
   sumMinv = (l2+j1).M() + (l1+j2).M();
 }
 else
   std::cout << "ERROR KLFitter::LikelihoodTopDilepton::CalculateMLepJet -------> NO VALID LEPTON CHARGE!!!" << std::endl;
 
 //  std::cout << "sumMinv: " << sumMinv << " pow(sumMinv," << alpha << "): " << pow(sumMinv,alpha) << std::endl;
 
 return norm*pow(sumMinv,alpha);
}


double KLFitter::LikelihoodTopDilepton::CalculateWeight(const std::vector<double> & parameters){
  
  double Weight=0.;
  
  //charged leptons
  TLorentzVector * l1 = new TLorentzVector();
  TLorentzVector * l2 = new TLorentzVector();
 
  // include parLep1E, parLep2E
  l1->SetPxPyPzE(lep1_fit_px,  lep1_fit_py,  lep1_fit_pz,  lep1_fit_e);
  l2->SetPxPyPzE(lep2_fit_px,  lep2_fit_py,  lep2_fit_pz,  lep2_fit_e);

  // jet1 and jet2:
  TLorentzVector * j1 = new TLorentzVector();
  TLorentzVector * j2 = new TLorentzVector();

  // include parB1E, parB2E
  j1->SetPxPyPzE(b1_fit_px, b1_fit_py, b1_fit_pz, b1_fit_e);
  j2->SetPxPyPzE(b2_fit_px, b2_fit_py, b2_fit_pz, b2_fit_e);

  //std::cout << "mtop= " << parameters[parTopM] << std::endl;
      
  Weight += CalculateWeightPerm(l1,l2,j1,j2,parameters);
  
  // if sumloglik option, sum over jet permutations
  if(doSumloglik)
    Weight += CalculateWeightPerm(l1,l2,j2,j1,parameters);
  
  delete l1;
  delete l2;
  delete j1;
  delete j2;

  return Weight;
}

double KLFitter::LikelihoodTopDilepton::CalculateWeightPerm(TLorentzVector * l1, TLorentzVector * l2, TLorentzVector * j1, TLorentzVector * j2, const std::vector<double> & parameters) {
  
  double weight=0.;
  int NSolutions=0;
  
  //set of kin solutions for neutrino/antineutrino
  NuSolutions nus;
  NuSolutions nubars;

  //std::cout << "parameters[parNuEta]: " << parameters[parNuEta] << " parameters[parAntiNuEta]: " << parameters[parAntiNuEta] << std::endl;
  
  // ensure correctly (lepton, nu) pair according to lepton charge
  if (lep1_meas_charge == 1 && lep2_meas_charge == -1){
    //std::cout << "opt1: lep1_meas_charge: " << lep1_meas_charge << " and lep2_meas_charge: " << lep2_meas_charge  << std::endl;
    nus = SolveForNuMom(l1,j1,parameters[parTopM],parameters[parNuEta]);
    nubars = SolveForNuMom(l2,j2,parameters[parTopM],parameters[parAntiNuEta]);
  }
  else if (lep1_meas_charge == -1 && lep2_meas_charge == 1){
    //std::cout << "opt2: lep1_meas_charge: " << lep1_meas_charge << " and lep2_meas_charge: " << lep2_meas_charge  << std::endl;
    nus = SolveForNuMom(l2,j1,parameters[parTopM],parameters[parNuEta]);
    nubars = SolveForNuMom(l1,j2,parameters[parTopM],parameters[parAntiNuEta]);
  }
  else
    std::cout << "ERROR KLFitter::LikelihoodTopDilepton::CalculateWeightPerm -------> NO VALID LEPTON CHARGE!!!" << std::endl;

  
  if (nus.NSolutions>0 && nubars.NSolutions>0) {
    weight += neutrino_weight(nus.nu1,nubars.nu1);  // ***solution 1
    ++NSolutions;
    
    if (nus.NSolutions==1 && nubars.NSolutions==2) {
      weight += neutrino_weight(nus.nu1,nubars.nu2);// ***solution 2
      
    }
    else if (nus.NSolutions==2 && nubars.NSolutions==1) {
      weight += neutrino_weight(nus.nu2,nubars.nu1);// ***solution 2
      
    }
    else if (nus.NSolutions==2 && nubars.NSolutions==2) {
      weight += neutrino_weight(nus.nu1,nubars.nu2);// ***solution 2
      
      weight += neutrino_weight(nus.nu2,nubars.nu1);// ***solution 3
      
      weight += neutrino_weight(nus.nu2,nubars.nu2);// ***solution 4
      
      NSolutions+=3;
    }
  }
  
  //std::cout << "final weight: " << weight << std::endl;
  return weight;
}

double KLFitter::LikelihoodTopDilepton::GaussNuEta(std::vector<double> parameters) {

  double weight=0.;

  double nueta_sigma=0.;
  double mtop_power=1.;
  double mtop=0.;

  mtop= parameters[parTopM];

  //build sigma of neutrino pseudorapidity dependency on mtop
  for (uint i=0;i!=nueta_params.size();++i) {
    nueta_sigma+=nueta_params[i]*mtop_power;
    mtop_power*=mtop;
  }

  weight = 1/(nueta_sigma*sqrt(2*M_PI))*exp( - parameters[parNuEta]*parameters[parNuEta]/(2*nueta_sigma*nueta_sigma));

  //std::cout << "parameters[parNuEta]: " << parameters[parNuEta] << " weight: " << weight << std::endl;
  //std::cout << "sigma1: " << nueta_sigma << " mtop: " << mtop << std::endl;

  return weight;
}

double KLFitter::LikelihoodTopDilepton::GaussAntiNuEta(std::vector<double> parameters) {

  double weight=0.;

  double nueta_sigma=0.;
  double mtop_power=1.;
  double mtop=0.;

  mtop= parameters[parTopM];

  //build sigma of antineutrino pseudorapidity dependency on mtop
  for (uint i=0;i!=nueta_params.size();++i) {
    nueta_sigma+=nueta_params[i]*mtop_power;
    mtop_power*=mtop;
  }

  weight =  1/(nueta_sigma*sqrt(2*M_PI))*exp( - parameters[parAntiNuEta]*parameters[parAntiNuEta]/(2*nueta_sigma*nueta_sigma));

  //std::cout << "parameters[parAntiNuEta]: " << parameters[parAntiNuEta] << " weight: " << weight << std::endl;
  //std::cout << "sigma2: " << nueta_sigma << " mtop: " << mtop << std::endl;

  return weight;
}

KLFitter::NuSolutions KLFitter::LikelihoodTopDilepton::SolveForNuMom(TLorentzVector * l, TLorentzVector * b, double mtop, double nueta){
  
  NuSolutions ret;
  double Wmass = fPhysicsConstants->MassW();
  double Wmass2=Wmass*Wmass;
  
  double bmass = fPhysicsConstants->MassBottom();
  
  double coshnueta;
  double sinhnueta;
  double Elprime;
  double Ebprime;
  double A;
  double B;
  double par1;
  double C;
  double par2;
  double D;
  double F;
  double det;
  double py;
  double px;
  double pT2;
  double pz;
  double tmp;
  double py1;
  double py2;
  double px1;
  double px2;
  double pT2_1;
  double pT2_2;
  double pz1;
  double pz2;
  double E1;
  double E2;


  coshnueta = cosh(nueta);
  sinhnueta = sinh(nueta);
  Elprime = l->E()*coshnueta-l->Pz()*sinhnueta;
  Ebprime = b->E()*coshnueta-b->Pz()*sinhnueta;

  //std::cout << "Elprime: " << Elprime << std::endl;
  //std::cout << "Ebprime: " << Ebprime << std::endl;

  A = (l->Py()*Ebprime-b->Py()*Elprime)/(b->Px()*Elprime-l->Px()*Ebprime);
  B = (Elprime*(mtop*mtop-Wmass2-bmass*bmass-2.*(*l)*(*b))-Ebprime*Wmass2)/(2.*(l->Px()*Ebprime-b->Px()*Elprime));

  par1 = (l->Px()*A+l->Py())/Elprime;
  C = A*A+1.-par1*par1;
  par2 = (Wmass2/2.+l->Px()*B)/Elprime;
  D = 2.*(A*B-par2*par1);
  F = B*B-par2*par2;
  det = D*D-4.*C*F;

  if (det<0.) {
    ret.NSolutions=0;
  }
  else if (det==0.) {
    ret.NSolutions=1;
    py=-D/(2.*C);
    px=A*py+B;
    pT2=px*px+py*py;
    pz=sqrt(pT2)*sinhnueta;
    ret.nu1.SetPxPyPzE(px,py,pz,sqrt(pT2+pz*pz));
  }
  else {
    ret.NSolutions=2;
    tmp=sqrt(det)/(2.*C);
    py1=-D/(2.*C)+tmp;
    py2=-D/(2.*C)-tmp;
    px1=A*py1+B;
    px2=A*py2+B;
    pT2_1=px1*px1+py1*py1;
    pT2_2=px2*px2+py2*py2;
    pz1=sqrt(pT2_1)*sinhnueta;
    pz2=sqrt(pT2_2)*sinhnueta;
    E1 = sqrt(pT2_1+pz1*pz1);
    E2 = sqrt(pT2_2+pz2*pz2);
    //std::cout << "py1= " << py1 << " py2= " << py2 << " px1= " << px1 << " px2= " << px2 << " pz1= " << pz1 << " pz2= " << pz2 << " E1= " << E1 << " E2= " << E2 << std::endl;
    //std::cout << "pT2_1= " << pT2_1 << " pT2_2= " << pT2_2 << " sqrt(pT2_1+pz1*pz1)= " << sqrt(pT2_1+pz1*pz1) << " sqrt(pT2_2+pz2*pz2)= " << sqrt(pT2_2+pz2*pz2) << std::endl;

    ret.nu1.SetPxPyPzE(px1,py1,pz1,E1);
    ret.nu2.SetPxPyPzE(px2,py2,pz2,E2);
    
  }
  if (ret.NSolutions==1) std::cout<<"KLFitter::LikelihoodTopDilepton::SolveForNuMom  NSolutions==1 !!!!!!!!!!!!!!!!!"<<std::endl;

  return ret;
}

double KLFitter::LikelihoodTopDilepton::neutrino_weight(TLorentzVector nu,TLorentzVector nubar) {

  static double sigmaX;
  static double sigmaY;
  static double dx;
  static double dy;
  
  // MET resolution in terms of SumET
  sigmaX = fResMET->GetSigma(SumET);
  sigmaY = fResMET->GetSigma(SumET);

  //std::cout << "sigmaX: " << sigmaX << " sigmaY: " << sigmaY << std::endl;

  dx = ETmiss_x-nu.Px()-nubar.Px();  // check!!
  dy = ETmiss_y-nu.Py()-nubar.Py();  // check!!
  
  //std::cout << "nu.Px(): " << nu.Px() << " nubar.Px(): " << nubar.Px() << std::endl;
  //std::cout << "nu.Py(): " << nu.Py() << " nubar.Py(): " << nubar.Py() << std::endl;

  //std::cout << "dx: " << dx << "|| dy: " << dy << std::endl;

  return exp( -dx*dx/(2.*sigmaX*sigmaX)  - dy*dy/(2.*sigmaY*sigmaY) );
}

// --------------------------------------------------------- 
std::vector<double> KLFitter::LikelihoodTopDilepton::GetInitialParameters()
{
  std::vector<double> values(GetNParameters());

  //std::cout << "GetNParameters: " << GetNParameters() << std::endl;

  // mtop
  values[parTopM] = fPhysicsConstants->MassTop();  
  
  // energies of the quarks
  values[parB1E] = b1_meas_e;
  values[parB2E] = b2_meas_e;

  // energy of the lepton
  if (fTypeLepton_1 == kElectron && fTypeLepton_2 == kMuon) {
    values[parLep1E] = (*fParticlesPermuted)->Electron(0)->E();
    values[parLep2E] = (*fParticlesPermuted)->Muon(0)->E();
  }
  else if (fTypeLepton_1 == kElectron && fTypeLepton_2 == kElectron) {
    values[parLep1E] = (*fParticlesPermuted)->Electron(0)->E();
    values[parLep2E] = (*fParticlesPermuted)->Electron(1)->E();
  }
  else if (fTypeLepton_1 == kMuon && fTypeLepton_2 == kMuon) {
    values[parLep1E] = (*fParticlesPermuted)->Muon(0)->E();
    values[parLep2E] = (*fParticlesPermuted)->Muon(1)->E();
  }
  
  values[parAntiNuEta] = 0.;
  values[parNuEta] = 0.;
 
  
  return values;
}

// --------------------------------------------------------- 
bool KLFitter::LikelihoodTopDilepton::NoTFProblem(std::vector<double> parameters) {
  fTFgood = true;
  this->LogLikelihood(parameters);
  return fTFgood;
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopDilepton::SavePermutedParticles() {
  b1_meas_e      = (*fParticlesPermuted)->Parton(0)->E();
  b1_meas_deteta = (*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kParton);
  b1_meas_px     = (*fParticlesPermuted)->Parton(0)->Px();
  b1_meas_py     = (*fParticlesPermuted)->Parton(0)->Py();
  b1_meas_pz     = (*fParticlesPermuted)->Parton(0)->Pz();
  b1_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(0)->M(), fPhysicsConstants->MassBottom(), b1_meas_px, b1_meas_py, b1_meas_pz, b1_meas_e);
  b1_meas_p      = sqrt(b1_meas_e*b1_meas_e - b1_meas_m*b1_meas_m);

  b2_meas_e      = (*fParticlesPermuted)->Parton(1)->E();
  b2_meas_deteta = (*fParticlesPermuted)->DetEta(1, KLFitter::Particles::kParton);
  b2_meas_px     = (*fParticlesPermuted)->Parton(1)->Px();
  b2_meas_py     = (*fParticlesPermuted)->Parton(1)->Py();
  b2_meas_pz     = (*fParticlesPermuted)->Parton(1)->Pz();
  b2_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(1)->M(), fPhysicsConstants->MassBottom(), b2_meas_px, b2_meas_py, b2_meas_pz, b2_meas_e);
  b2_meas_p      = sqrt(b2_meas_e*b2_meas_e - b2_meas_m*b2_meas_m);

  TLorentzVector * lepton_1(0);
  TLorentzVector * lepton_2(0);
  
  
  if (fTypeLepton_1 == kElectron && fTypeLepton_2 == kMuon){
    lepton_1 = (*fParticlesPermuted)->Electron(0);
    lep1_meas_deteta = (*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kElectron);
    lep1_meas_charge = (*fParticlesPermuted)->LeptonCharge(0, KLFitter::Particles::kElectron);
    lepton_2 = (*fParticlesPermuted)->Muon(0);
    lep2_meas_deteta = (*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kMuon);
    lep2_meas_charge = (*fParticlesPermuted)->LeptonCharge(0, KLFitter::Particles::kMuon);
  }
  else if (fTypeLepton_1 == kElectron && fTypeLepton_2 == kElectron) {
    lepton_1 = (*fParticlesPermuted)->Electron(0);
    lep1_meas_deteta = (*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kElectron);
    lep1_meas_charge = (*fParticlesPermuted)->LeptonCharge(0, KLFitter::Particles::kElectron);
    lepton_2 = (*fParticlesPermuted)->Electron(1);
    lep2_meas_deteta = (*fParticlesPermuted)->DetEta(1, KLFitter::Particles::kElectron);
    lep2_meas_charge = (*fParticlesPermuted)->LeptonCharge(1, KLFitter::Particles::kElectron);
  }
  else if (fTypeLepton_1 == kMuon && fTypeLepton_2 == kMuon) {
    lepton_1 = (*fParticlesPermuted)->Muon(0);
    lep1_meas_deteta = (*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kMuon);
    lep1_meas_charge = (*fParticlesPermuted)->LeptonCharge(0, KLFitter::Particles::kMuon);
    lepton_2 = (*fParticlesPermuted)->Muon(1);
    lep2_meas_deteta = (*fParticlesPermuted)->DetEta(1, KLFitter::Particles::kMuon);
    lep2_meas_charge = (*fParticlesPermuted)->LeptonCharge(1, KLFitter::Particles::kMuon);
  }
  
  lep1_meas_e        = lepton_1->E();
  lep1_meas_sintheta = sin(lepton_1->Theta());
  lep1_meas_pt       = lepton_1->Pt();
  lep1_meas_px       = lepton_1->Px();
  lep1_meas_py       = lepton_1->Py();
  lep1_meas_pz       = lepton_1->Pz();

  lep2_meas_e        = lepton_2->E();
  lep2_meas_sintheta = sin(lepton_2->Theta());
  lep2_meas_pt       = lepton_2->Pt();
  lep2_meas_px       = lepton_2->Px();
  lep2_meas_py       = lepton_2->Py();
  lep2_meas_pz       = lepton_2->Pz();

  // no error
  return 1;
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTopDilepton::SaveResolutionFunctions() {

  fResEnergyB1 = (*fDetector)->ResEnergyBJet(b1_meas_deteta);
  fResEnergyB2 = (*fDetector)->ResEnergyBJet(b2_meas_deteta);
  
  if(fTypeLepton_1 == kElectron && fTypeLepton_2 == kMuon){
    fResLepton1 = (*fDetector)->ResEnergyElectron(lep1_meas_deteta);
    fResLepton2 = (*fDetector)->ResEnergyMuon(lep2_meas_deteta);
  }
  else if (fTypeLepton_1 == kElectron && fTypeLepton_2 == kElectron) {
    fResLepton1 = (*fDetector)->ResEnergyElectron(lep1_meas_deteta);
    fResLepton2 = (*fDetector)->ResEnergyElectron(lep2_meas_deteta);
  }
  else if (fTypeLepton_1 == kMuon && fTypeLepton_2 == kMuon) {
    fResLepton1 = (*fDetector)->ResEnergyMuon(lep1_meas_deteta);
    fResLepton2 = (*fDetector)->ResEnergyMuon(lep2_meas_deteta);
  }
  
  fResMET = (*fDetector)->ResMissingET();
  
  // no error
  return 1;
}
// --------------------------------------------------------- 

int KLFitter::LikelihoodTopDilepton::BuildModelParticles() {	
  if (GetBestFitParameters().size() > 0) CalculateLorentzVectors(GetBestFitParameters());

  TLorentzVector * b1 = fParticlesModel->Parton(0);
  TLorentzVector * b2 = fParticlesModel->Parton(1);
  TLorentzVector * lep1(0);
  TLorentzVector * lep2(0);

  if(fTypeLepton_1 == kElectron && fTypeLepton_2 == kMuon){
    lep1  = fParticlesModel->Electron(0);
    lep2  = fParticlesModel->Muon(0);
  }
  else if (fTypeLepton_1 == kElectron && fTypeLepton_2 == kElectron) {
    lep1  = fParticlesModel->Electron(0);
    lep2  = fParticlesModel->Electron(1);
  }
  else {
    lep1  = fParticlesModel->Muon(0);
    lep2  = fParticlesModel->Muon(1);
  }
  
  b1   ->SetPxPyPzE(b1_fit_px, b1_fit_py, b1_fit_pz, b1_fit_e);
  b2   ->SetPxPyPzE(b2_fit_px, b2_fit_py, b2_fit_pz, b2_fit_e);
  lep1 ->SetPxPyPzE(lep1_fit_px,  lep1_fit_py,  lep1_fit_pz,  lep1_fit_e);
  lep2 ->SetPxPyPzE(lep2_fit_px,  lep2_fit_py,  lep2_fit_pz,  lep2_fit_e);
  
  // no error
  return 1;
}

// --------------------------------------------------------- 

std::vector<double> KLFitter::LikelihoodTopDilepton::LogLikelihoodComponents(std::vector<double> parameters)
{
  std::vector<double> vecci(0);

  // calculate 4-vectors 
  CalculateLorentzVectors(parameters); 

  // NuWT likelihood term
  double nuwt_weight(0.);
  
  // temporary flag for a safe use of the transfer functions
  bool TFgoodTmp(true);
  
  // NuWT weight
  nuwt_weight = CalculateWeight(parameters);
  if(nuwt_weight == 0.)
    vecci.push_back(log(1e-99));
  else
    vecci.push_back(log( CalculateWeight(parameters) )); //comp0

  // jet energy resolution terms 
  if( fResEnergyB1->p(b1_fit_e, b1_meas_e, TFgoodTmp) == 0. )
    vecci.push_back(log(1e-99));
  else
    vecci.push_back(log( fResEnergyB1->p(b1_fit_e, b1_meas_e, TFgoodTmp) ));   //comp1
  if (!TFgoodTmp) fTFgood = false;


  if( fResEnergyB2->p(b2_fit_e, b2_meas_e, TFgoodTmp) == 0.)
    vecci.push_back(log(1e-99));
  else
    vecci.push_back(log( fResEnergyB2->p(b2_fit_e, b2_meas_e, TFgoodTmp) ));   //comp2
  if (!TFgoodTmp) fTFgood = false;
  
  // lepton energy resolution terms 
  if (fTypeLepton_1 == kElectron && fTypeLepton_2 == kMuon) {
    if( fResLepton1->p(lep1_fit_e, lep1_meas_e, TFgoodTmp) == 0.) 
      vecci.push_back(log(1e-99));
    else
      vecci.push_back(log( fResLepton1->p(lep1_fit_e, lep1_meas_e, TFgoodTmp) ));                         //comp3

    if( fResLepton2->p(lep2_fit_e*lep2_meas_sintheta, lep2_meas_pt, TFgoodTmp) == 0.)
      vecci.push_back(log(1e-99));
    else
      vecci.push_back(log( fResLepton2->p(lep2_fit_e* lep2_meas_sintheta, lep2_meas_pt, TFgoodTmp) ));    //comp4
    if (!TFgoodTmp) fTFgood = false;
  }
  else if (fTypeLepton_1 == kElectron && fTypeLepton_2 == kElectron) {
    if( fResLepton1->p(lep1_fit_e, lep1_meas_e, TFgoodTmp) == 0.)
      vecci.push_back(log(1e-99));
    else
      vecci.push_back(log( fResLepton1->p(lep1_fit_e, lep1_meas_e, TFgoodTmp) )); //comp3

    if( fResLepton2->p(lep2_fit_e, lep2_meas_e, TFgoodTmp) == 0.)
      vecci.push_back(log(1e-99));
    else
      vecci.push_back(log( fResLepton2->p(lep2_fit_e, lep2_meas_e, TFgoodTmp) )); //comp4
    if (!TFgoodTmp) fTFgood = false;
  }
  else if (fTypeLepton_1 == kMuon && fTypeLepton_2 == kMuon) {
    if( fResLepton1->p(lep1_fit_e*lep1_meas_sintheta, lep1_meas_pt, TFgoodTmp) == 0.)
      vecci.push_back(log(1e-99));
    else
      vecci.push_back(log( fResLepton1->p(lep1_fit_e* lep1_meas_sintheta, lep1_meas_pt, TFgoodTmp) )); //comp3

    if( fResLepton2->p(lep2_fit_e*lep2_meas_sintheta, lep2_meas_pt, TFgoodTmp) == 0.)
      vecci.push_back(log(1e-99));
    else
      vecci.push_back(log( fResLepton2->p(lep2_fit_e* lep2_meas_sintheta, lep2_meas_pt, TFgoodTmp) )); //comp4
    if (!TFgoodTmp) fTFgood = false;
  }

  // nueta and antinueta terms
  
  if (GaussAntiNuEta(parameters) == 0.)
    vecci.push_back(log(1e-99));
  else
    vecci.push_back(log (GaussAntiNuEta(parameters)) ); //comp5

  if (GaussNuEta(parameters) == 0.)
    vecci.push_back(log(1e-99));
  else
    vecci.push_back(log (GaussNuEta(parameters)) ); //comp6

  // sum of invariant masses (lep,jet) term

 if (CalculateMLepJet(parameters) == 0.){
   vecci.push_back(log(1e-99));
 }
 else 
   vecci.push_back(log( CalculateMLepJet(parameters)) ); //comp7


  // return log of likelihood 
  return vecci; 
}

// --------------------------------------------------------- 

void KLFitter::LikelihoodTopDilepton::MCMCIterationInterface()
{
  TLorentzVector  MCMC_b1(0.);
  TLorentzVector  MCMC_b2(0.);
  TLorentzVector  MCMC_lep1(0.);
  TLorentzVector  MCMC_lep2(0.);
  TLorentzVector  MCMC_nu1(0.);
  TLorentzVector  MCMC_nu2(0.);

  TLorentzVector  MCMC_lep(0.);
  TLorentzVector  MCMC_antilep(0.);

  double scale_b1(0.);
  double scale_b2(0.);
  double scale_l1(0.);
  double scale_l2(0.);

  NuSolutions nus;
  NuSolutions nubars;

  // for mttbar
  double mttbar(0.);

  // for costheta
  std::pair<float, float> costheta(0.,0.);
  std::vector <TLorentzVector> * help_ParticleVector  = new std::vector<TLorentzVector>(0);

  // get number of chains
  int nchains = MCMCGetNChains();

  // get number of parameters
  int npar = GetNParameters();

  // loop over all chains and fill histogram
  for (int i = 0; i < nchains; ++i) {
    // get the current values of the KLFitter parameters. These are
    // stored in fMCMCx.
    double mtop = fMCMCx.at(i * npar + 0);
    double Eb1 = fMCMCx.at(i * npar + 1);
    double Eb2 = fMCMCx.at(i * npar + 2);
    double El1 = fMCMCx.at(i * npar + 3);
    double El2 = fMCMCx.at(i * npar + 4);
    double antinueta = fMCMCx.at(i * npar + 5);
    double nueta = fMCMCx.at(i * npar + 6);

    //         std::cout << "mtop: " << mtop << std::endl;
    //         std::cout << "Eb1: " << Eb1 << std::endl;
    //         std::cout << "Eb2: " << Eb2 << std::endl;
    //         std::cout << "El1: " << El1 << std::endl;
    //         std::cout << "El2: " << El2 << std::endl;
    // 	       std::cout << "antinueta: " << antinueta<< std::endl;
    //         std::cout << "nueta: " << nueta << std::endl;
    
    
    scale_b1 = sqrt(Eb1*Eb1 - b1_meas_m*b1_meas_m) / b1_meas_p;
    scale_b2 = sqrt(Eb2*Eb2 - b2_meas_m*b2_meas_m) / b2_meas_p;
    scale_l1 = El1 / lep1_meas_e;
    scale_l2 = El2 / lep2_meas_e;
  
    MCMC_b1  .SetPxPyPzE(scale_b1 * b1_meas_px,scale_b1 * b1_meas_py,scale_b1 * b1_meas_pz,Eb1);
    MCMC_b2  .SetPxPyPzE(scale_b2 * b2_meas_px,scale_b2 * b2_meas_py,scale_b2 * b2_meas_pz,Eb2);
    MCMC_lep1.SetPxPyPzE(scale_l1 * lep1_meas_px,scale_l1 * lep1_meas_py,scale_l1 * lep1_meas_pz,El1);
    MCMC_lep2.SetPxPyPzE(scale_l2 * lep2_meas_px,scale_l2 * lep2_meas_py,scale_l2 * lep2_meas_pz,El2);
   

    // ensure correctly (lepton, nu) pair according to lepton charge
    if (lep1_meas_charge == 1 && lep2_meas_charge == -1){
      //std::cout << "opt1: lep1_meas_charge: " << lep1_meas_charge << " and lep2_meas_charge: " << lep2_meas_charge  << std::endl;
      nus = SolveForNuMom(&MCMC_lep1,&MCMC_b1,mtop,nueta);
      nubars = SolveForNuMom(&MCMC_lep2,&MCMC_b2,mtop,antinueta);
      //for costheta
      MCMC_lep = MCMC_lep2;
      MCMC_antilep = MCMC_lep1;
    }
    else if (lep1_meas_charge == -1 && lep2_meas_charge == 1){
      //std::cout << "opt2: lep1_meas_charge: " << lep1_meas_charge << " and lep2_meas_charge: " << lep2_meas_charge  << std::endl;
      nus = SolveForNuMom(&MCMC_lep2,&MCMC_b1,mtop,nueta);
      nubars = SolveForNuMom(&MCMC_lep1,&MCMC_b2,mtop,antinueta);
      //for costheta
      MCMC_lep = MCMC_lep1;
      MCMC_antilep = MCMC_lep2;
    }
    else
      std::cout << "ERROR KLFitter::LikelihoodTopDilepton::MCMCIterationInterface -------> NO VALID LEPTON CHARGE!!!" << std::endl;
    

    if (nus.NSolutions > 0 && nubars.NSolutions > 0) { // 1 nu 1 nubar
      
      if (nus.nu1.M() >= 0 && nubars.nu1.M() >= 0) {
	//mttbar
	mttbar = (MCMC_b1 + MCMC_b2 + MCMC_lep1 + MCMC_lep2 + nus.nu1 + nubars.nu1).M();
	//std::cout << "1 mttbar: " << mttbar << std::endl;
	fHistMttbar->GetHistogram()->Fill(mttbar);
	//costheta
	help_ParticleVector->clear();
	help_ParticleVector -> push_back(MCMC_lep);
	help_ParticleVector -> push_back(MCMC_antilep);
	help_ParticleVector -> push_back(nus.nu1);
	help_ParticleVector -> push_back(nubars.nu1);
	help_ParticleVector -> push_back(MCMC_b1);
	help_ParticleVector -> push_back(MCMC_b2);
	costheta = CalculateCosTheta( help_ParticleVector );
	fHistCosTheta->GetHistogram()->Fill(costheta.first);
	fHistCosTheta->GetHistogram()->Fill(costheta.second);
	//std::cout << "costheta: " << costheta.first << std::endl;
      }

      if (nus.NSolutions == 1 && nubars.NSolutions == 2) { // 1 nu 2 nubar
	if (nus.nu1.M() >= 0 && nubars.nu2.M() >= 0) {
	  //mttbar
	  mttbar = (MCMC_b1 + MCMC_b2 + MCMC_lep1 + MCMC_lep2 + nus.nu1 + nubars.nu2).M();
	  //std::cout << "2 mttbar: " << mttbar << std::endl;
	  fHistMttbar->GetHistogram()->Fill(mttbar);
	  //costheta
	  help_ParticleVector->clear();
	  help_ParticleVector -> push_back(MCMC_lep);
	  help_ParticleVector -> push_back(MCMC_antilep);
	  help_ParticleVector -> push_back(nus.nu1);
	  help_ParticleVector -> push_back(nubars.nu2);
	  help_ParticleVector -> push_back(MCMC_b1);
	  help_ParticleVector -> push_back(MCMC_b2);
	  costheta = CalculateCosTheta( help_ParticleVector );
	  fHistCosTheta->GetHistogram()->Fill(costheta.first);
	  fHistCosTheta->GetHistogram()->Fill(costheta.second);
	  //std::cout << "costheta: " << costheta.first << std::endl;
	}	
      }// 1 nu 2 nubar
      else if (nus.NSolutions == 2 && nubars.NSolutions == 1){  // 2 nu 1 nubar
	if (nus.nu2.M() >= 0 && nubars.nu1.M() >= 0) {
	  //mttbar
	  mttbar = (MCMC_b1 + MCMC_b2 + MCMC_lep1 + MCMC_lep2 + nus.nu2 + nubars.nu1).M();
	  //std::cout << "3 mttbar: " << mttbar << std::endl;
	  fHistMttbar->GetHistogram()->Fill(mttbar);
	  //costheta
	  help_ParticleVector->clear();
	  help_ParticleVector -> push_back(MCMC_lep);
	  help_ParticleVector -> push_back(MCMC_antilep);
	  help_ParticleVector -> push_back(nus.nu2);
	  help_ParticleVector -> push_back(nubars.nu1);
	  help_ParticleVector -> push_back(MCMC_b1);
	  help_ParticleVector -> push_back(MCMC_b2);
	  costheta = CalculateCosTheta( help_ParticleVector );
	  fHistCosTheta->GetHistogram()->Fill(costheta.first);
	  fHistCosTheta->GetHistogram()->Fill(costheta.second);
	  //std::cout << "costheta: " << costheta.first << std::endl;
	}
      }// 2 nu 1 nubar
      else if (nus.NSolutions == 2 && nubars.NSolutions == 2){  // 2 nu 2 nubar
	
	if (nus.nu1.M() >= 0 && nubars.nu2.M() >= 0) {
	  //mttbar
	  mttbar = (MCMC_b1 + MCMC_b2 + MCMC_lep1 + MCMC_lep2 + nus.nu1 + nubars.nu2).M();
	  //std::cout << "4a mttbar: " << mttbar << std::endl;
	  fHistMttbar->GetHistogram()->Fill(mttbar);
	  //costheta
	  help_ParticleVector->clear();
	  help_ParticleVector -> push_back(MCMC_lep);
	  help_ParticleVector -> push_back(MCMC_antilep);
	  help_ParticleVector -> push_back(nus.nu1);
	  help_ParticleVector -> push_back(nubars.nu2);
	  help_ParticleVector -> push_back(MCMC_b1);
	  help_ParticleVector -> push_back(MCMC_b2);
	  costheta = CalculateCosTheta( help_ParticleVector );
	  fHistCosTheta->GetHistogram()->Fill(costheta.first);
	  fHistCosTheta->GetHistogram()->Fill(costheta.second);
	  //std::cout << "costheta: " << costheta.first << std::endl;
	}

	if (nus.nu2.M() >= 0 && nubars.nu1.M() >= 0) {
	  //mttbar
	  mttbar = (MCMC_b1 + MCMC_b2 + MCMC_lep1 + MCMC_lep2 + nus.nu2 + nubars.nu1).M();
	  //	std::cout << "4b mttbar: " << mttbar << std::endl;
	  fHistMttbar->GetHistogram()->Fill(mttbar);
	  //costheta
	  help_ParticleVector->clear();
	  help_ParticleVector -> push_back(MCMC_lep);
	  help_ParticleVector -> push_back(MCMC_antilep);
	  help_ParticleVector -> push_back(nus.nu2);
	  help_ParticleVector -> push_back(nubars.nu1);
	  help_ParticleVector -> push_back(MCMC_b1);
	  help_ParticleVector -> push_back(MCMC_b2);
	  costheta = CalculateCosTheta( help_ParticleVector );
	  fHistCosTheta->GetHistogram()->Fill(costheta.first);
	  fHistCosTheta->GetHistogram()->Fill(costheta.second);
	  //std::cout << "costheta: " << costheta.first << std::endl;
	}

	if (nus.nu2.M() >= 0 && nubars.nu2.M() >= 0) {
	  //mttbar
	  mttbar = (MCMC_b1 + MCMC_b2 + MCMC_lep1 + MCMC_lep2 + nus.nu2 + nubars.nu2).M();
	  //std::cout << "4c mttbar: " << mttbar << std::endl;
	  fHistMttbar->GetHistogram()->Fill(mttbar);
	  //costheta
	  help_ParticleVector->clear();
	  help_ParticleVector -> push_back(MCMC_lep);
	  help_ParticleVector -> push_back(MCMC_antilep);
	  help_ParticleVector -> push_back(nus.nu2);
	  help_ParticleVector -> push_back(nubars.nu2);
	  help_ParticleVector -> push_back(MCMC_b1);
	  help_ParticleVector -> push_back(MCMC_b2);
	  costheta = CalculateCosTheta( help_ParticleVector );
	  fHistCosTheta->GetHistogram()->Fill(costheta.first);
	  fHistCosTheta->GetHistogram()->Fill(costheta.second);
	  //std::cout << "costheta: " << costheta.first << std::endl;
	}
      }// 2 nu 2 nubar
    }//Nsol

  }//Nchains
  delete help_ParticleVector;
}

// ---------------------------------------------------------
std::pair<float, float> KLFitter::LikelihoodTopDilepton::CalculateCosTheta(std::vector <TLorentzVector> *particles)
{

  TLorentzVector         lep(0.0, 0.0, 0.0, 0.0);
  TLorentzVector     antilep(0.0, 0.0, 0.0, 0.0);
  TLorentzVector          nu(0.0, 0.0, 0.0, 0.0);
  TLorentzVector       nubar(0.0, 0.0, 0.0, 0.0);
  TLorentzVector           b(0.0, 0.0, 0.0, 0.0);
  TLorentzVector        bbar(0.0, 0.0, 0.0, 0.0);
 

  TLorentzVector       top(0.0, 0.0, 0.0, 0.0);
  TLorentzVector   antitop(0.0, 0.0, 0.0, 0.0);
  TLorentzVector     Wplus(0.0, 0.0, 0.0, 0.0);
  TLorentzVector    Wminus(0.0, 0.0, 0.0, 0.0);

  TVector3       lep3(0.0, 0.0, 0.0);
  TVector3   antilep3(0.0, 0.0, 0.0);
  TVector3   Wplus_bo(0.0, 0.0, 0.0);
  TVector3  Wminus_bo(0.0, 0.0, 0.0);
  TVector3         b3(0.0, 0.0, 0.0);
  TVector3      bbar3(0.0, 0.0, 0.0);
 

  lep         = particles -> at(0);
  antilep     = particles -> at(1);
  nu          = particles -> at(2);
  nubar       = particles -> at(3);
  b           = particles -> at(4);
  bbar        = particles -> at(5);
 
  //std::cout << "nu Px: " << nu.Px() << " nu Py: " << nu.Py() << std::endl;
  //std::cout << "nubar Px: " << nubar.Px() << " nubar Py: " << nubar.Py() << std::endl;
   

  Wplus     = antilep   +      nu;
  Wminus    = lep       +      nubar;
  top       = Wplus     +      b;
  antitop   = Wminus    +      bbar;

  //std::cout << "Wplus->Pt(): " << Wplus.Pt() << " Wminus->Pt(): " << Wminus.Pt() << std::endl;

  Wplus_bo  = Wplus.BoostVector();
  Wminus_bo = Wminus.BoostVector();

  b.Boost(-Wplus_bo);
  antilep.Boost(-Wplus_bo);

  bbar.Boost(-Wminus_bo);
  lep.Boost(-Wminus_bo);

  lep3.SetXYZ(lep.Px(), lep.Py(),     lep.Pz());
  antilep3.SetXYZ(antilep.Px(), antilep.Py(),     antilep.Pz());
  b3.SetXYZ(b.Px(), b.Py(), b.Pz());
  bbar3.SetXYZ(bbar.Px(), bbar.Py(), bbar.Pz());

 //  std::cout << "lep.Px(): " << lep.Px() << " lep.Py(): "<< lep.Py() << std::endl;
//   std::cout << "antilep.Px(): " << antilep.Px() << " antilep.Py(): "<< antilep.Py() << std::endl;
//   std::cout << "b.Px(): " << b.Px() << " b.Py(): "<< b.Py() << std::endl;
//   std::cout << "bbar.Px(): " << bbar.Px() << " bbar.Py(): "<< bbar.Py() << std::endl;

 
  float cos_theta_top          = cos(antilep3.Angle(-b3));
  float cos_theta_antitop      = cos(lep3.Angle(-bbar3));  

  std::pair<float, float> cos;

  cos.first  = cos_theta_top;
  cos.second = cos_theta_antitop;

  //std::cout << "cos_theta_top: "  << cos_theta_top << " and cos_theta_antitop: " << cos_theta_antitop << std::endl;

  return cos;

}
// ---------------------------------------------------------
