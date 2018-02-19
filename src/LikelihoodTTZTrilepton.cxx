/*
 * Copyright (c) 2009--2018, the KLFitter developer team
 *
 * This file is part of KLFitter.
 *
 * KLFitter is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * KLFitter is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with KLFitter. If not, see <http://www.gnu.org/licenses/>.
 */

#include "KLFitter/LikelihoodTTZTrilepton.h"
#include "KLFitter/ResolutionBase.h"
#include "KLFitter/Particles.h"
#include "KLFitter/Permutations.h"
#include "KLFitter/PhysicsConstants.h"
#include "KLFitter/DetectorBase.h"

#include <iostream>
#include <algorithm>

#include <BAT/BCMath.h>
#include "BAT/BCParameter.h"

// ---------------------------------------------------------
KLFitter::LikelihoodTTZTrilepton::LikelihoodTTZTrilepton() : KLFitter::LikelihoodBase::LikelihoodBase()
  , fFlagTopMassFixed(false)
  , fFlagUseJetMass(false)
  , fFlagGetParSigmasFromTFs(false)
  , ETmiss_x(0.)
  , ETmiss_y(0.)
  , SumET(0.)
  , fTypeLepton(kElectron)
  , fTFgood(true)
  , fInvMassCutoff(5.)
  , fOnShellFraction(0.869) {
  // define model particles
  this->DefineModelParticles();

  // define parameters
  this->DefineParameters();
}

// ---------------------------------------------------------
KLFitter::LikelihoodTTZTrilepton::~LikelihoodTTZTrilepton() {
}

// ---------------------------------------------------------
int KLFitter::LikelihoodTTZTrilepton::SetET_miss_XY_SumET(double etx, double ety, double sumet) {
  // set missing ET x and y component and the SumET
  ETmiss_x = etx;
  ETmiss_y = ety;
  SumET = sumet;

  // no error
  return 1;
}

// ---------------------------------------------------------
void KLFitter::LikelihoodTTZTrilepton::SetLeptonType(LeptonType leptontype) {
  if (leptontype != kElectron && leptontype != kMuon) {
    std::cout << "KLFitter::SetLeptonTyp(). Warning: lepton type not defined. Set electron as lepton type." << std::endl;
    fTypeLepton = kElectron;
  } else {
    fTypeLepton = leptontype;
  }

  // define model particles
  DefineModelParticles();
}

// ---------------------------------------------------------
void KLFitter::LikelihoodTTZTrilepton::SetLeptonType(int leptontype) {
  if (leptontype != 1 && leptontype != 2) {
    std::cout << "KLFitter::SetLeptonTyp(). Warning: lepton type not defined. Set electron as lepton type." << std::endl;
    leptontype = 1;
  }

  if (leptontype == 1) {
    SetLeptonType(kElectron);
  } else if (leptontype == 2) {
    SetLeptonType(kMuon);
  }
}

// ---------------------------------------------------------
int KLFitter::LikelihoodTTZTrilepton::DefineModelParticles() {
  // check if model particles and lorentz vector container exist and delete
  if (fParticlesModel) {
    delete fParticlesModel;
    fParticlesModel = 0;
  }

  // create the particles of the model
  fParticlesModel = new KLFitter::Particles();

  // add model particles
  // create dummy TLorentzVector
  TLorentzVector * dummy = new TLorentzVector(0,0,0,0);  // 4-vector
  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kParton,  // type
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
                               KLFitter::Particles::kLight);  // light jet (truth)

  fParticlesModel->AddParticle(dummy,
                               KLFitter::Particles::kParton,
                               "light quark 2",
                               3,                            // index of corresponding particle
                               KLFitter::Particles::kLight);  // light jet (truth)

  if (fTypeLepton == kElectron) {
    fParticlesModel->AddParticle(dummy,
                                 KLFitter::Particles::kElectron,
                                 "electron");
    fParticlesModel->AddParticle(dummy,
                                 KLFitter::Particles::kElectron,
                                 "electron Z1");
    fParticlesModel->AddParticle(dummy,
                                 KLFitter::Particles::kElectron,
                                 "electron Z2");
  } else if (fTypeLepton == kMuon) {
    fParticlesModel->AddParticle(dummy,
                                 KLFitter::Particles::kMuon,
                                 "muon");
    fParticlesModel->AddParticle(dummy,
                                 KLFitter::Particles::kMuon,
                                 "muon Z1");
    fParticlesModel->AddParticle(dummy,
                                 KLFitter::Particles::kMuon,
                                 "muon Z2");
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
                               "Z boson");

  // free memory
  delete dummy;

  // no error
  return 1;
}

// ---------------------------------------------------------
void KLFitter::LikelihoodTTZTrilepton::DefineParameters() {
  // add parameters of model
  AddParameter("energy hadronic b",       fPhysicsConstants->MassBottom(), 1000.0);  // parBhadE
  AddParameter("energy leptonic b",       fPhysicsConstants->MassBottom(), 1000.0);  // parBlepE
  AddParameter("energy light quark 1",    0.0, 1000.0);                             // parLQ1E
  AddParameter("energy light quark 2",    0.0, 1000.0);                             // parLQ2E
  AddParameter("energy lepton",           0.0, 1000.0);                             // parLepE
  AddParameter("energy Z lepton 1",       0.0, 1000.0);                             // parLepZ1E
  AddParameter("energy Z lepton 2",       0.0, 1000.0);                             // parLepZ2E
  AddParameter("p_x neutrino",        -1000.0, 1000.0);                             // parNuPx
  AddParameter("p_y neutrino",        -1000.0, 1000.0);                             // parNuPy
  AddParameter("p_z neutrino",        -1000.0, 1000.0);                             // parNuPz
  AddParameter("top mass",              100.0, 1000.0);                             // parTopM
  AddParameter("Z mass",                  0.0, 1000.0);                             // parZM
}

// ---------------------------------------------------------
double KLFitter::LikelihoodTTZTrilepton::LogBreitWignerRelNorm(const double& x, const double& mean, const double& gamma) {
  double g = std::sqrt(pow(mean, 2) * (pow(mean, 2) + pow(gamma, 2)));
  double k = (2 * std::sqrt(2) * mean * gamma * g) / (M_PI * std::sqrt(pow(mean, 2) + g));
  double f = k / (pow(pow(x, 2) - pow(mean, 2), 2) + pow(mean * gamma, 2));
  return log(f);
}

// ---------------------------------------------------------
double KLFitter::LikelihoodTTZTrilepton::LogZCombinedDistribution(const double& x, const double& mean, const double& gamma) {
  // note: This catches exceptions when the variables are set to non-sensible
  // values. If there is any common way to handle exceptions, it should be
  // implemented here.
  const auto fraction = fOnShellFraction;
  if (fraction < 0 || fraction > 1) throw;
  if (fInvMassCutoff < 0) throw;

  double g = std::sqrt(pow(mean, 2) * (pow(mean, 2) + pow(gamma, 2)));
  double k = (2 * std::sqrt(2) * mean * gamma * g) / (M_PI * std::sqrt(pow(mean, 2) + g));
  double on_shell = k / (pow(pow(x, 2) - pow(mean, 2), 2) + pow(mean * gamma, 2));
  double off_shell = fInvMassCutoff / x / x;
  return log(on_shell * fraction + off_shell * (1 - fraction));
}

// ---------------------------------------------------------
int KLFitter::LikelihoodTTZTrilepton::CalculateLorentzVectors(std::vector <double> const& parameters) {
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

  static double Z_fit_e;
  static double Z_fit_px;
  static double Z_fit_py;
  static double Z_fit_pz;

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

  // Z lepton 1
  lepZ1_fit_e = parameters[parLepZ1E];
  scale = lepZ1_fit_e / lepZ1_meas_e;
  lepZ1_fit_px = scale * lepZ1_meas_px;
  lepZ1_fit_py = scale * lepZ1_meas_py;
  lepZ1_fit_pz = scale * lepZ1_meas_pz;

  // Z lepton 2
  lepZ2_fit_e = parameters[parLepZ2E];
  scale = lepZ2_fit_e / lepZ2_meas_e;
  lepZ2_fit_px = scale * lepZ2_meas_px;
  lepZ2_fit_py = scale * lepZ2_meas_py;
  lepZ2_fit_pz = scale * lepZ2_meas_pz;

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

  // Z boson
  Z_fit_e  = lepZ1_fit_e  + lepZ2_fit_e;
  Z_fit_px = lepZ1_fit_px + lepZ2_fit_px;
  Z_fit_py = lepZ1_fit_py + lepZ2_fit_py;
  Z_fit_pz = lepZ1_fit_pz + lepZ2_fit_pz;
  Z_fit_m  = sqrt(Z_fit_e*Z_fit_e - (Z_fit_px*Z_fit_px + Z_fit_py*Z_fit_py + Z_fit_pz*Z_fit_pz));

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::LikelihoodTTZTrilepton::Initialize() {
  // error code
  int err = 1;

  // save the current permuted particles
  err *= SavePermutedParticles();

  // save the corresponding resolution functions
  err *= SaveResolutionFunctions();

  // adjust parameter ranges
  err *= AdjustParameterRanges();

  // set initial values
  // (only for Markov chains - initial parameters for other minimisation methods are set in Fitter.cxx)
  SetInitialParameters(GetInitialParameters());

  // return error code
  return err;
}

// ---------------------------------------------------------
int KLFitter::LikelihoodTTZTrilepton::RemoveInvariantParticlePermutations() {
  // error code
  int err = 1;

  // remove the permutation from the second and the third jet
  KLFitter::Particles::ParticleType ptype = KLFitter::Particles::kParton;
  std::vector<int> indexVector_Jets;
  indexVector_Jets.push_back(2);
  indexVector_Jets.push_back(3);
  err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Jets);

  // remove the permutation from the two Z leptons
  KLFitter::Particles::ParticleType ptypeLepZ;
  if (fTypeLepton == kElectron) {
    ptypeLepZ = KLFitter::Particles::kElectron;
  } else {
    ptypeLepZ = KLFitter::Particles::kMuon;
  }

  std::vector<int> indexVector_LepZ;
  indexVector_LepZ.push_back(1);
  indexVector_LepZ.push_back(2);
  err *= (*fPermutations)->InvariantParticlePermutations(ptypeLepZ, indexVector_LepZ);

  // remove invariant jet permutations of notevent jets
  KLFitter::Particles* particles = (*fPermutations)->Particles();
  indexVector_Jets.clear();
  for (int iPartons = 4; iPartons < particles->NPartons(); iPartons++) {
    indexVector_Jets.push_back(iPartons);
  }
  err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Jets);

  // remove the permutation from the other lepton
  if (fTypeLepton == kElectron) {
    ptype = KLFitter::Particles::kMuon;
    std::vector<int> indexVector_Muons;
    for (int iMuon = 0; iMuon < particles->NMuons(); iMuon++) {
      indexVector_Muons.push_back(iMuon);
    }
    err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Muons);
  } else if (fTypeLepton == kMuon) {
    ptype = KLFitter::Particles::kElectron;
    std::vector<int> indexVector_Electrons;
    for (int iElectron = 0; iElectron < particles->NElectrons(); iElectron++) {
      indexVector_Electrons.push_back(iElectron);
    }
    err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Electrons);
  }

  // return error code
  return err;
}

// ---------------------------------------------------------
int KLFitter::LikelihoodTTZTrilepton::RemoveForbiddenParticlePermutations() {
  // error code
  int err = 1;

  // only in b-tagging type kVetoNoFit
  if (!((fBTagMethod == kVetoNoFit) || (fBTagMethod == kVetoNoFitLight) || (fBTagMethod == kVetoNoFitBoth)))
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
      if ((fBTagMethod == kVetoNoFit)&&((!isBtagged) || (trueFlavor != KLFitter::Particles::kLight)))
        continue;
      if ((fBTagMethod == kVetoNoFitLight)&&((isBtagged) || (trueFlavor != KLFitter::Particles::kB)))
        continue;
      if ((fBTagMethod == kVetoNoFitBoth)&&(((isBtagged)&&(trueFlavor != KLFitter::Particles::kLight)) || ((!isBtagged)&&(trueFlavor != KLFitter::Particles::kB))))
        continue;

      err *= (*fPermutations)->RemoveParticlePermutations(KLFitter::Particles::kParton, iParton, iPartonModel);
    }
  }

  // return error code
  return err;
}

// ---------------------------------------------------------
int KLFitter::LikelihoodTTZTrilepton::AdjustParameterRanges() {
  // adjust limits
  double nsigmas_jet    = fFlagGetParSigmasFromTFs ? 10 : 7;
  double nsigmas_lepton = fFlagGetParSigmasFromTFs ? 10 : 2;
  double nsigmas_met    = fFlagGetParSigmasFromTFs ? 10 : 1;

  double E = (*fParticlesPermuted)->Parton(0)->E();
  double m = fPhysicsConstants->MassBottom();
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(0)->M());
  double sigma = fFlagGetParSigmasFromTFs ? fResEnergyBhad->GetSigma(E) : sqrt(E);
  double Emin = std::max(m, E - nsigmas_jet* sigma);
  double Emax  = E + nsigmas_jet* sigma;
  SetParameterRange(parBhadE, Emin, Emax);

  E = (*fParticlesPermuted)->Parton(1)->E();
  m = fPhysicsConstants->MassBottom();
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(1)->M());
  sigma = fFlagGetParSigmasFromTFs ? fResEnergyBlep->GetSigma(E) : sqrt(E);
  Emin = std::max(m, E - nsigmas_jet* sigma);
  Emax  = E + nsigmas_jet* sigma;
  SetParameterRange(parBlepE, Emin, Emax);

  E = (*fParticlesPermuted)->Parton(2)->E();
  m = 0.001;
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(2)->M());
  sigma = fFlagGetParSigmasFromTFs ? fResEnergyLQ1->GetSigma(E) : sqrt(E);
  Emin = std::max(m, E - nsigmas_jet* sigma);
  Emax  = E + nsigmas_jet* sigma;
  SetParameterRange(parLQ1E, Emin, Emax);

  E = (*fParticlesPermuted)->Parton(3)->E();
  m = 0.001;
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(3)->M());
  sigma = fFlagGetParSigmasFromTFs ? fResEnergyLQ2->GetSigma(E) : sqrt(E);
  Emin = std::max(m, E - nsigmas_jet* sigma);
  Emax  = E + nsigmas_jet* sigma;
  SetParameterRange(parLQ2E, Emin, Emax);

  if (fTypeLepton == kElectron) {
    E = (*fParticlesPermuted)->Electron(0)->E();
    sigma = fFlagGetParSigmasFromTFs ? fResLepton->GetSigma(E) : sqrt(E);
    Emin = std::max(0.001, E - nsigmas_lepton* sigma);
    Emax  = E + nsigmas_lepton* sigma;
  } else if (fTypeLepton == kMuon) {
    E = (*fParticlesPermuted)->Muon(0)->E();
    double sintheta = sin((*fParticlesPermuted)->Muon(0)->Theta());
    sigma = fFlagGetParSigmasFromTFs ? fResLepton->GetSigma(E*sintheta)/sintheta : E*E*sintheta;
    double sigrange = nsigmas_lepton* sigma;
    Emin = std::max(0.001,E -sigrange);
    Emax = E +sigrange;
  }
  SetParameterRange(parLepE, Emin, Emax);

  if (fTypeLepton == kElectron) {
    E = (*fParticlesPermuted)->Electron(1)->E();
    sigma = fFlagGetParSigmasFromTFs ? fResLeptonZ1->GetSigma(E) : sqrt(E);
    Emin = std::max(0.001, E - nsigmas_lepton* sigma);
    Emax  = E + nsigmas_lepton* sigma;
  } else if (fTypeLepton == kMuon) {
    E = (*fParticlesPermuted)->Muon(1)->E();
    double sintheta = sin((*fParticlesPermuted)->Muon(1)->Theta());
    sigma = fFlagGetParSigmasFromTFs ? fResLeptonZ1->GetSigma(E*sintheta)/sintheta : E*E*sintheta;
    double sigrange = nsigmas_lepton* sigma;
    Emin = std::max(0.001,E -sigrange);
    Emax = E +sigrange;
  }
  SetParameterRange(parLepZ1E, Emin, Emax);

  if (fTypeLepton == kElectron) {
    E = (*fParticlesPermuted)->Electron(2)->E();
    sigma = fFlagGetParSigmasFromTFs ? fResLeptonZ2->GetSigma(E) : sqrt(E);
    Emin = std::max(0.001, E - nsigmas_lepton* sigma);
    Emax  = E + nsigmas_lepton* sigma;
  } else if (fTypeLepton == kMuon) {
    E = (*fParticlesPermuted)->Muon(2)->E();
    double sintheta = sin((*fParticlesPermuted)->Muon(2)->Theta());
    sigma = fFlagGetParSigmasFromTFs ? fResLeptonZ2->GetSigma(E*sintheta)/sintheta : E*E*sintheta;
    double sigrange = nsigmas_lepton* sigma;
    Emin = std::max(0.001,E -sigrange);
    Emax = E +sigrange;
  }
  SetParameterRange(parLepZ2E, Emin, Emax);

  // note: this is hard-coded at the moment

  sigma = fFlagGetParSigmasFromTFs ? fResMET->GetSigma(SumET) : 100;
  double sigrange = nsigmas_met*sigma;
  SetParameterRange(parNuPx, ETmiss_x-sigrange, ETmiss_x+sigrange);
  SetParameterRange(parNuPy, ETmiss_y-sigrange, ETmiss_y+sigrange);

  if (fFlagTopMassFixed)
    SetParameterRange(parTopM, fPhysicsConstants->MassTop(), fPhysicsConstants->MassTop());

  SetParameterRange(parZM, fPhysicsConstants->MassZ(), fPhysicsConstants->MassZ());

  // no error
  return 1;
}

// ---------------------------------------------------------
double KLFitter::LikelihoodTTZTrilepton::LogLikelihood(const std::vector<double> & parameters) {
  // calculate 4-vectors
  CalculateLorentzVectors(parameters);

  // define log of likelihood
  double logprob(0.);

  // temporary flag for a safe use of the transfer functions
  bool TFgoodTmp(true);

  // jet energy resolution terms
  logprob += log(fResEnergyBhad->p(bhad_fit_e, bhad_meas_e, TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;

  logprob += log(fResEnergyBlep->p(blep_fit_e, blep_meas_e, TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;

  logprob += log(fResEnergyLQ1->p(lq1_fit_e, lq1_meas_e, TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;

  logprob += log(fResEnergyLQ2->p(lq2_fit_e, lq2_meas_e, TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;

  // lepton energy resolution terms
  if (fTypeLepton == kElectron) {
    logprob += log(fResLepton->p(lep_fit_e, lep_meas_e, TFgoodTmp));
  } else if (fTypeLepton == kMuon) {
    logprob += log(fResLepton->p(lep_fit_e* lep_meas_sintheta, lep_meas_pt, TFgoodTmp));
  }
  if (!TFgoodTmp) fTFgood = false;

  if (fTypeLepton == kElectron) {
    logprob += log(fResLeptonZ1->p(lepZ1_fit_e, lepZ1_meas_e, TFgoodTmp));
  } else if (fTypeLepton == kMuon) {
    logprob += log(fResLeptonZ1->p(lepZ1_fit_e* lepZ1_meas_sintheta, lepZ1_meas_pt, TFgoodTmp));
  }
  if (!TFgoodTmp) fTFgood = false;

  if (fTypeLepton == kElectron) {
    logprob += log(fResLeptonZ2->p(lepZ2_fit_e, lepZ2_meas_e, TFgoodTmp));
  } else if (fTypeLepton == kMuon) {
    logprob += log(fResLeptonZ2->p(lepZ2_fit_e* lepZ2_meas_sintheta, lepZ2_meas_pt, TFgoodTmp));
  }
  if (!TFgoodTmp) fTFgood = false;

  // neutrino px and py
  logprob += log(fResMET->p(nu_fit_px, ETmiss_x, TFgoodTmp, SumET));
  if (!TFgoodTmp) fTFgood = false;

  logprob += log(fResMET->p(nu_fit_py, ETmiss_y, TFgoodTmp, SumET));
  if (!TFgoodTmp) fTFgood = false;

  // physics constants
  double massW = fPhysicsConstants->MassW();
  double gammaW = fPhysicsConstants->GammaW();
  // note: top mass width should be made DEPENDENT on the top mass at a certain point
  //    fPhysicsConstants->SetMassTop(parameters[parTopM]);
  // (this will also set the correct width for the top)
  double gammaTop = fPhysicsConstants->GammaTop();

  double gammaZ = fPhysicsConstants->GammaZ();

  // note: as opposed to the LikelihoodTopLeptonJets class, we use a normalised
  // version of the Breit-Wigner here to make sure that weightings between
  // functions are handled correctly.

  // Breit-Wigner of hadronically decaying W-boson
  logprob += LogBreitWignerRelNorm(whad_fit_m, massW, gammaW);

  // Breit-Wigner of leptonically decaying W-boson
  logprob += LogBreitWignerRelNorm(wlep_fit_m, massW, gammaW);

  // Breit-Wigner of hadronically decaying top quark
  logprob += LogBreitWignerRelNorm(thad_fit_m, parameters[parTopM], gammaTop);

  // Breit-Wigner of leptonically decaying top quark
  logprob += LogBreitWignerRelNorm(tlep_fit_m, parameters[parTopM], gammaTop);

  // Breit-Wigner of Z boson decaying into two leptons
  logprob += LogZCombinedDistribution(Z_fit_m, parameters[parZM], gammaZ);

  // return log of likelihood
  return logprob;
}

// ---------------------------------------------------------
std::vector<double> KLFitter::LikelihoodTTZTrilepton::GetInitialParameters() {
  std::vector<double> values(GetInitialParametersWoNeutrinoPz());

  // check second neutrino solution
  std::vector<double> neutrino_pz_solutions = GetNeutrinoPzSolutions();
  if (int(neutrino_pz_solutions.size()) == 1) {
    values[parNuPz] = neutrino_pz_solutions[0];
  } else if (int(neutrino_pz_solutions.size()) == 2) {
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
std::vector<double> KLFitter::LikelihoodTTZTrilepton::GetInitialParametersWoNeutrinoPz() {
  std::vector<double> values(GetNParameters());

  // energies of the quarks
  values[parBhadE] = bhad_meas_e;
  values[parBlepE] = blep_meas_e;
  values[parLQ1E]  = lq1_meas_e;
  values[parLQ2E]  = lq2_meas_e;
  values[parLepZ1E]  = lepZ1_meas_e;
  values[parLepZ2E]  = lepZ2_meas_e;

  // energy of the lepton
  if (fTypeLepton == kElectron) {
    values[parLepE] = (*fParticlesPermuted)->Electron(0)->E();
  } else if (fTypeLepton == kMuon) {
    values[parLepE] = (*fParticlesPermuted)->Muon(0)->E();
  }

  if (fTypeLepton == kElectron) {
    values[parLepZ1E] = (*fParticlesPermuted)->Electron(1)->E();
  } else if (fTypeLepton == kMuon) {
    values[parLepZ1E] = (*fParticlesPermuted)->Muon(1)->E();
  }

  if (fTypeLepton == kElectron) {
    values[parLepZ2E] = (*fParticlesPermuted)->Electron(2)->E();
  } else if (fTypeLepton == kMuon) {
    values[parLepZ2E] = (*fParticlesPermuted)->Muon(2)->E();
  }

  // missing px and py
  values[parNuPx] = ETmiss_x;
  values[parNuPy] = ETmiss_y;

  // pz of the neutrino
  values[parNuPz] = 0.;

  // top mass
  double mtop = (*(*fParticlesPermuted)->Parton(0) + *(*fParticlesPermuted)->Parton(2) + *(*fParticlesPermuted)->Parton(3)).M();
  if (mtop < GetParameter(parTopM)->GetLowerLimit()) {
    mtop = GetParameter(parTopM)->GetLowerLimit();
  } else if (mtop > GetParameter(parTopM)->GetUpperLimit()) {
    mtop = GetParameter(parTopM)->GetUpperLimit();
  }
  values[parTopM] = mtop;

  // Z mass
  double mz;
  if (fTypeLepton == kElectron) {
    mz = (*(*fParticlesPermuted)->Electron(1) + *(*fParticlesPermuted)->Electron(2)).M();
  } else {
    mz = (*(*fParticlesPermuted)->Muon(1) + *(*fParticlesPermuted)->Muon(2)).M();
  }
  if (mz < GetParameter(parZM)->GetLowerLimit()) {
    mz = GetParameter(parZM)->GetLowerLimit();
  } else if (mz > GetParameter(parZM)->GetUpperLimit()) {
    mz = GetParameter(parZM)->GetUpperLimit();
  }
  values[parZM] = mz;

  // return the vector
  return values;
}

// ---------------------------------------------------------
std::vector<double> KLFitter::LikelihoodTTZTrilepton::GetNeutrinoPzSolutions() {
  return CalculateNeutrinoPzSolutions();
}

// ---------------------------------------------------------
std::vector<double> KLFitter::LikelihoodTTZTrilepton::CalculateNeutrinoPzSolutions(TLorentzVector* additionalParticle) {
  std::vector<double> pz;

  KLFitter::PhysicsConstants constants;
  // electron mass
  double mE = 0.;

  double px_c = 0.0;
  double py_c = 0.0;
  double pz_c = 0.0;
  double Ec = 0.0;

  if (fTypeLepton == kElectron) {
    px_c = (*fParticlesPermuted)->Electron(0)->Px();
    py_c = (*fParticlesPermuted)->Electron(0)->Py();
    pz_c = (*fParticlesPermuted)->Electron(0)->Pz();
    Ec = (*fParticlesPermuted)->Electron(0)->E();
  } else if (fTypeLepton == kMuon) {
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
  if (squareRoot < 1.e-6) {
    pz.push_back(pz_offset);
  } else {
    pz.push_back(pz_offset + squareRoot / (2*a));
    pz.push_back(pz_offset - squareRoot / (2*a));
  }

  return pz;
}

// ---------------------------------------------------------
bool KLFitter::LikelihoodTTZTrilepton::NoTFProblem(std::vector<double> parameters) {
  fTFgood = true;
  this->LogLikelihood(parameters);
  return fTFgood;
}

// ---------------------------------------------------------
int KLFitter::LikelihoodTTZTrilepton::SavePermutedParticles() {
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

  TLorentzVector * leptonZ1(0);
  TLorentzVector * leptonZ2(0);
  if (fTypeLepton == kElectron) {
    leptonZ1 = (*fParticlesPermuted)->Electron(1);
    lepZ1_meas_deteta = (*fParticlesPermuted)->DetEta(1, KLFitter::Particles::kElectron);
    leptonZ2 = (*fParticlesPermuted)->Electron(2);
    lepZ2_meas_deteta = (*fParticlesPermuted)->DetEta(2, KLFitter::Particles::kElectron);
  } else {
    leptonZ1 = (*fParticlesPermuted)->Muon(1);
    lepZ1_meas_deteta = (*fParticlesPermuted)->DetEta(1, KLFitter::Particles::kMuon);
    leptonZ2 = (*fParticlesPermuted)->Muon(2);
    lepZ2_meas_deteta = (*fParticlesPermuted)->DetEta(2, KLFitter::Particles::kMuon);
  }

  lepZ1_meas_e        = leptonZ1->E();
  lepZ1_meas_sintheta = sin(leptonZ1->Theta());
  lepZ1_meas_pt       = leptonZ1->Pt();
  lepZ1_meas_px       = leptonZ1->Px();
  lepZ1_meas_py       = leptonZ1->Py();
  lepZ1_meas_pz       = leptonZ1->Pz();

  lepZ2_meas_e        = leptonZ2->E();
  lepZ2_meas_sintheta = sin(leptonZ2->Theta());
  lepZ2_meas_pt       = leptonZ2->Pt();
  lepZ2_meas_px       = leptonZ2->Px();
  lepZ2_meas_py       = leptonZ2->Py();
  lepZ2_meas_pz       = leptonZ2->Pz();

  TLorentzVector * lepton(0);
  if (fTypeLepton == kElectron) {
    lepton = (*fParticlesPermuted)->Electron(0);
    lep_meas_deteta = (*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kElectron);
  } else {
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
int KLFitter::LikelihoodTTZTrilepton::SaveResolutionFunctions() {
  fResEnergyBhad = (*fDetector)->ResEnergyBJet(bhad_meas_deteta);
  fResEnergyBlep = (*fDetector)->ResEnergyBJet(blep_meas_deteta);
  fResEnergyLQ1  = (*fDetector)->ResEnergyLightJet(lq1_meas_deteta);
  fResEnergyLQ2  = (*fDetector)->ResEnergyLightJet(lq2_meas_deteta);
  if (fTypeLepton == kElectron) {
    fResLepton = (*fDetector)->ResEnergyElectron(lep_meas_deteta);
  } else if (fTypeLepton == kMuon) {
    fResLepton = (*fDetector)->ResEnergyMuon(lep_meas_deteta);
  }
  fResMET = (*fDetector)->ResMissingET();

  if (fTypeLepton == kElectron) {
    fResLeptonZ1 = (*fDetector)->ResEnergyElectron(lepZ1_meas_deteta);
  } else if (fTypeLepton == kMuon) {
    fResLeptonZ1 = (*fDetector)->ResEnergyMuon(lepZ1_meas_deteta);
  }

  if (fTypeLepton == kElectron) {
    fResLeptonZ2 = (*fDetector)->ResEnergyElectron(lepZ2_meas_deteta);
  } else if (fTypeLepton == kMuon) {
    fResLeptonZ2 = (*fDetector)->ResEnergyMuon(lepZ2_meas_deteta);
  }
  fResMET = (*fDetector)->ResMissingET();

  // no error
  return 1;
}
// ---------------------------------------------------------

int KLFitter::LikelihoodTTZTrilepton::BuildModelParticles() {
  if (GetBestFitParameters().size() > 0) CalculateLorentzVectors(GetBestFitParameters());

  TLorentzVector * bhad = fParticlesModel->Parton(0);
  TLorentzVector * blep = fParticlesModel->Parton(1);
  TLorentzVector * lq1  = fParticlesModel->Parton(2);
  TLorentzVector * lq2  = fParticlesModel->Parton(3);
  TLorentzVector * lep(0);
  if (fTypeLepton == kElectron) {
    lep  = fParticlesModel->Electron(0);
  } else if (fTypeLepton == kMuon) {
    lep  = fParticlesModel->Muon(0);
  }
  TLorentzVector * lepZ1(0);
  if (fTypeLepton == kElectron) {
    lepZ1  = fParticlesModel->Electron(1);
  } else if (fTypeLepton == kMuon) {
    lepZ1  = fParticlesModel->Muon(1);
  }
  TLorentzVector * lepZ2(0);
  if (fTypeLepton == kElectron) {
    lepZ2  = fParticlesModel->Electron(2);
  } else if (fTypeLepton == kMuon) {
    lepZ2  = fParticlesModel->Muon(2);
  }
  TLorentzVector * nu   = fParticlesModel->Neutrino(0);
  TLorentzVector * whad  = fParticlesModel->Boson(0);
  TLorentzVector * wlep  = fParticlesModel->Boson(1);
  TLorentzVector * thad  = fParticlesModel->Parton(4);
  TLorentzVector * tlep  = fParticlesModel->Parton(5);

  TLorentzVector * Z = fParticlesModel->Boson(2);

  bhad->SetPxPyPzE(bhad_fit_px, bhad_fit_py, bhad_fit_pz, bhad_fit_e);
  blep->SetPxPyPzE(blep_fit_px, blep_fit_py, blep_fit_pz, blep_fit_e);
  lq1 ->SetPxPyPzE(lq1_fit_px,  lq1_fit_py,  lq1_fit_pz,  lq1_fit_e);
  lq2 ->SetPxPyPzE(lq2_fit_px,  lq2_fit_py,  lq2_fit_pz,  lq2_fit_e);
  lepZ1->SetPxPyPzE(lepZ1_fit_px, lepZ1_fit_py, lepZ1_fit_pz, lepZ1_fit_e);
  lepZ2->SetPxPyPzE(lepZ2_fit_px, lepZ2_fit_py, lepZ2_fit_pz, lepZ2_fit_e);
  lep ->SetPxPyPzE(lep_fit_px,  lep_fit_py,  lep_fit_pz,  lep_fit_e);
  nu  ->SetPxPyPzE(nu_fit_px,   nu_fit_py,   nu_fit_pz,   nu_fit_e);

  (*whad) = (*lq1)  + (*lq2);
  (*wlep) = (*lep)  + (*nu);
  (*thad) = (*whad) + (*bhad);
  (*tlep) = (*wlep) + (*blep);

  (*Z) = (*lepZ1) + (*lepZ2);

  // no error
  return 1;
}
// ---------------------------------------------------------

// ---------------------------------------------------------
std::vector<double> KLFitter::LikelihoodTTZTrilepton::LogLikelihoodComponents(std::vector<double> parameters) {
  std::vector<double> vecci;

  // calculate 4-vectors
  CalculateLorentzVectors(parameters);


  // temporary flag for a safe use of the transfer functions
  bool TFgoodTmp(true);

  // jet energy resolution terms
  vecci.push_back(log(fResEnergyBhad->p(bhad_fit_e, bhad_meas_e, TFgoodTmp)));  // comp0
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log(fResEnergyBlep->p(blep_fit_e, blep_meas_e, TFgoodTmp)));  // comp1
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log(fResEnergyLQ1->p(lq1_fit_e, lq1_meas_e, TFgoodTmp)));  // comp2
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log(fResEnergyLQ2->p(lq2_fit_e, lq2_meas_e, TFgoodTmp)));  // comp3
  if (!TFgoodTmp) fTFgood = false;

  // lepton energy resolution terms
  if (fTypeLepton == kElectron) {
    vecci.push_back(log(fResLepton->p(lep_fit_e, lep_meas_e, TFgoodTmp)));  // comp4
  } else if (fTypeLepton == kMuon) {
    vecci.push_back(log(fResLepton->p(lep_fit_e* lep_meas_sintheta, lep_meas_pt, TFgoodTmp)));  // comp4
  }
  if (!TFgoodTmp) fTFgood = false;

  if (fTypeLepton == kElectron) {
    vecci.push_back(log(fResLeptonZ1->p(lepZ1_fit_e, lepZ1_meas_e, TFgoodTmp)));  // comp4
  } else if (fTypeLepton == kMuon) {
    vecci.push_back(log(fResLeptonZ1->p(lepZ1_fit_e* lepZ1_meas_sintheta, lepZ1_meas_pt, TFgoodTmp)));  // comp4
  }
  if (!TFgoodTmp) fTFgood = false;

  if (fTypeLepton == kElectron) {
    vecci.push_back(log(fResLeptonZ2->p(lepZ2_fit_e, lepZ2_meas_e, TFgoodTmp)));  // comp4
  } else if (fTypeLepton == kMuon) {
    vecci.push_back(log(fResLeptonZ2->p(lepZ2_fit_e* lepZ2_meas_sintheta, lepZ2_meas_pt, TFgoodTmp)));  // comp4
  }
  if (!TFgoodTmp) fTFgood = false;

  // neutrino px and py
  vecci.push_back(log(fResMET->p(nu_fit_px, ETmiss_x, TFgoodTmp, SumET)));  // comp5
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log(fResMET->p(nu_fit_py, ETmiss_y, TFgoodTmp, SumET)));  // comp6
  if (!TFgoodTmp) fTFgood = false;

  // physics constants
  double massW = fPhysicsConstants->MassW();
  double gammaW = fPhysicsConstants->GammaW();
  // note: top mass width should be made DEPENDENT on the top mass at a certain point
  //    fPhysicsConstants->SetMassTop(parameters[parTopM]);
  // (this will also set the correct width for the top)
  double gammaTop = fPhysicsConstants->GammaTop();

  double gammaZ = fPhysicsConstants->GammaZ();

  // note: as opposed to the LikelihoodTopLeptonJets class, we use a normalised
  // version of the Breit-Wigner here to make sure that weightings between
  // functions are handled correctly.

  // Breit-Wigner of hadronically decaying W-boson
  vecci.push_back(LogBreitWignerRelNorm(whad_fit_m, massW, gammaW));  // comp7

  // Breit-Wigner of leptonically decaying W-boson
  vecci.push_back(LogBreitWignerRelNorm(wlep_fit_m, massW, gammaW));  // comp8

  // Breit-Wigner of hadronically decaying top quark
  vecci.push_back(LogBreitWignerRelNorm(thad_fit_m, parameters[parTopM], gammaTop));  // comp9

  // Breit-Wigner of leptonically decaying top quark
  vecci.push_back(LogBreitWignerRelNorm(tlep_fit_m, parameters[parTopM], gammaTop));  // comp10

  // Breit-Wigner of Z decaying into 2 leptons
  vecci.push_back(LogZCombinedDistribution(Z_fit_m, parameters[parZM], gammaZ));  // comp11

  // return log of likelihood
  return vecci;
}
