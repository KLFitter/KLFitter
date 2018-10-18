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

#include "KLFitter/BoostedLikelihoodTopLeptonJets.h"

#include <algorithm>
#include <iostream>

#include "BAT/BCMath.h"
#include "BAT/BCParameter.h"
#include "KLFitter/DetectorBase.h"
#include "KLFitter/ParticleCollection.h"
#include "KLFitter/Permutations.h"
#include "KLFitter/PhysicsConstants.h"
#include "KLFitter/ResolutionBase.h"
#include "TLorentzVector.h"

// ---------------------------------------------------------
KLFitter::BoostedLikelihoodTopLeptonJets::BoostedLikelihoodTopLeptonJets()
  : KLFitter::LikelihoodBase::LikelihoodBase()
  , fFlagTopMassFixed(false)
  , fFlagGetParSigmasFromTFs(false)
  , ETmiss_x(0.)
  , ETmiss_y(0.)
  , SumET(0.)
  , fTypeLepton(kElectron) {
  // define model particles
  this->DefineModelParticles();

  // define parameters
  this->DefineParameters();
}

// ---------------------------------------------------------
KLFitter::BoostedLikelihoodTopLeptonJets::~BoostedLikelihoodTopLeptonJets() = default;

// ---------------------------------------------------------
int KLFitter::BoostedLikelihoodTopLeptonJets::SetET_miss_XY_SumET(double etx, double ety, double sumet) {
  // set missing ET x and y component and the SumET
  ETmiss_x = etx;
  ETmiss_y = ety;
  SumET = sumet;

  // no error
  return 1;
}

// ---------------------------------------------------------
void KLFitter::BoostedLikelihoodTopLeptonJets::SetLeptonType(LeptonType leptontype) {
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
void KLFitter::BoostedLikelihoodTopLeptonJets::SetLeptonType(int leptontype) {
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
int KLFitter::BoostedLikelihoodTopLeptonJets::DefineModelParticles() {
  // create the particles of the model
  fParticlesModel.reset(new KLFitter::ParticleCollection{});

  // add model particles
  // create dummy TLorentzVector
  TLorentzVector dummy{0, 0, 0, 0};  // 4-vector

  fParticlesModel->AddParticle(&dummy,
                               Particle::Type::kParton,  // type
                               "hadronic b quark",            // name
                               0,                             // index of corresponding particle
                               Particle::JetTrueFlavor::kB);      // b jet (truth)

  fParticlesModel->AddParticle(&dummy,
                               Particle::Type::kParton,
                               "leptonic b quark",
                               1,                             // index of corresponding particle
                               Particle::JetTrueFlavor::kB);      // b jet (truth)

  fParticlesModel->AddParticle(&dummy,
                               Particle::Type::kParton,
                               "light quarks",
                               2,                             // index of corresponding particle
                               Particle::JetTrueFlavor::kLight);  // merged light jet (truth)

  if (fTypeLepton == kElectron) {
    fParticlesModel->AddParticle(&dummy,
                                 Particle::Type::kElectron,
                                 "electron");
  } else if (fTypeLepton == kMuon) {
    fParticlesModel->AddParticle(&dummy,
                                 Particle::Type::kMuon,
                                 "muon");
  }

  fParticlesModel->AddParticle(&dummy,
                               Particle::Type::kNeutrino,
                               "neutrino");

  fParticlesModel->AddParticle(&dummy,
                               Particle::Type::kBoson,
                               "leptonic W");

  fParticlesModel->AddParticle(&dummy,
                               Particle::Type::kParton,
                               "hadronic top");

  fParticlesModel->AddParticle(&dummy,
                               Particle::Type::kParton,
                               "leptonic top");

  // no error
  return 1;
}

// ---------------------------------------------------------
void KLFitter::BoostedLikelihoodTopLeptonJets::DefineParameters() {
  // add parameters of model
  AddParameter("energy hadronic b",       fPhysicsConstants.MassBottom(), 1000.0);  // parBhadE
  AddParameter("energy leptonic b",       fPhysicsConstants.MassBottom(), 1000.0);  // parBlepE
  AddParameter("energy light quarks",     fPhysicsConstants.MassW(), 1000.0);       // parLQE
  AddParameter("energy lepton",           0.0, 1000.0);                              // parLepE
  AddParameter("p_x neutrino",        -1000.0, 1000.0);                              // parNuPx
  AddParameter("p_y neutrino",        -1000.0, 1000.0);                              // parNuPy
  AddParameter("p_z neutrino",        -1000.0, 1000.0);                              // parNuPz
  AddParameter("top mass",              100.0, 1000.0);                              // parTopM
}

// ---------------------------------------------------------
int KLFitter::BoostedLikelihoodTopLeptonJets::CalculateLorentzVectors(std::vector <double> const& parameters) {
  static double scale;
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

  // light quarks
  lq_fit_e = parameters[parLQE];
  scale = sqrt(lq_fit_e*lq_fit_e - lq_meas_m*lq_meas_m) / lq_meas_p;
  lq_fit_px = scale * lq_meas_px;
  lq_fit_py = scale * lq_meas_py;
  lq_fit_pz = scale * lq_meas_pz;

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

  // leptonic W
  wlep_fit_e  = lep_fit_e +nu_fit_e;
  wlep_fit_px = lep_fit_px+nu_fit_px;
  wlep_fit_py = lep_fit_py+nu_fit_py;
  wlep_fit_pz = lep_fit_pz+nu_fit_pz;
  wlep_fit_m = sqrt(wlep_fit_e*wlep_fit_e - (wlep_fit_px*wlep_fit_px + wlep_fit_py*wlep_fit_py + wlep_fit_pz*wlep_fit_pz));

  // hadronic top
  thad_fit_e = lq_fit_e+bhad_fit_e;
  thad_fit_px = lq_fit_px+bhad_fit_px;
  thad_fit_py = lq_fit_py+bhad_fit_py;
  thad_fit_pz = lq_fit_pz+bhad_fit_pz;
  thad_fit_m = sqrt(thad_fit_e*thad_fit_e - (thad_fit_px*thad_fit_px + thad_fit_py*thad_fit_py + thad_fit_pz*thad_fit_pz));

  // leptonic top
  tlep_fit_e = wlep_fit_e+blep_fit_e;
  tlep_fit_px = wlep_fit_px+blep_fit_px;
  tlep_fit_py = wlep_fit_py+blep_fit_py;
  tlep_fit_pz = wlep_fit_pz+blep_fit_pz;
  tlep_fit_m = sqrt(tlep_fit_e*tlep_fit_e - (tlep_fit_px*tlep_fit_px + tlep_fit_py*tlep_fit_py + tlep_fit_pz*tlep_fit_pz));

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::BoostedLikelihoodTopLeptonJets::RemoveInvariantParticlePermutations() {
  // error code
  int err = 1;

  // there are no invariant particle permutations with this likelihood
  Particle::Type ptype = Particle::Type::kParton;
  std::vector<int> indexVector_Jets;

  // remove invariant jet permutations of notevent jets
  KLFitter::ParticleCollection* particles = (*fPermutations)->Particles();
  indexVector_Jets.clear();
  for (int iPartons = 3; iPartons < particles->NPartons(); iPartons++) {
    indexVector_Jets.push_back(iPartons);
  }
  err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Jets);

  // remove the permutation from the other lepton
  if (fTypeLepton == kElectron) {
    ptype = Particle::Type::kMuon;
    std::vector<int> indexVector_Muons;
    for (int iMuon = 0; iMuon < particles->NMuons(); iMuon++) {
      indexVector_Muons.push_back(iMuon);
    }
    err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Muons);
  } else if (fTypeLepton == kMuon) {
    ptype = Particle::Type::kElectron;
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
int KLFitter::BoostedLikelihoodTopLeptonJets::AdjustParameterRanges() {
  // adjust limits
  double nsigmas_jet    = fFlagGetParSigmasFromTFs ? 10 : 7;
  double nsigmas_lepton = fFlagGetParSigmasFromTFs ? 10 : 2;
  double nsigmas_met    = fFlagGetParSigmasFromTFs ? 10 : 1;

  double E = (*fParticlesPermuted)->Parton(0)->E();
  double m = fPhysicsConstants.MassBottom();
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(0)->M());
  double sigma = fFlagGetParSigmasFromTFs ? fResEnergyBhad->GetSigma(E) : sqrt(E);
  double Emin = std::max(m, E - nsigmas_jet* sigma);
  double Emax  = E + nsigmas_jet* sigma;
  SetParameterRange(parBhadE, Emin, Emax);

  E = (*fParticlesPermuted)->Parton(1)->E();
  m = fPhysicsConstants.MassBottom();
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(1)->M());
  sigma = fFlagGetParSigmasFromTFs ? fResEnergyBlep->GetSigma(E) : sqrt(E);
  Emin = std::max(m, E - nsigmas_jet* sigma);
  Emax  = E + nsigmas_jet* sigma;
  SetParameterRange(parBlepE, Emin, Emax);

  E = (*fParticlesPermuted)->Parton(2)->E();
  m = fPhysicsConstants.MassW();
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(2)->M());
  sigma = fFlagGetParSigmasFromTFs ? fResEnergyLQ->GetSigma(E) : sqrt(E);
  Emin = std::max(m, E - nsigmas_jet* sigma);
  Emax  = E + nsigmas_jet* sigma;
  SetParameterRange(parLQE, Emin, Emax);

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
    Emin = std::max(0.001, E -sigrange);
    Emax = E +sigrange;
  }
  SetParameterRange(parLepE, Emin, Emax);

  // note: this is hard-coded in the momement
  sigma = fFlagGetParSigmasFromTFs ? fResMET->GetSigma(SumET) : 100;
  double sigrange = nsigmas_met*sigma;
  SetParameterRange(parNuPx, ETmiss_x-sigrange, ETmiss_x+sigrange);
  SetParameterRange(parNuPy, ETmiss_y-sigrange, ETmiss_y+sigrange);

  if (fFlagTopMassFixed)
    SetParameterRange(parTopM, fPhysicsConstants.MassTop(), fPhysicsConstants.MassTop());

  // no error
  return 1;
}

// ---------------------------------------------------------
void KLFitter::BoostedLikelihoodTopLeptonJets::RequestResolutionFunctions() {
  (*fDetector)->RequestResolutionType(ResolutionType::EnergyLightJet);
  (*fDetector)->RequestResolutionType(ResolutionType::EnergyBJet);
  (*fDetector)->RequestResolutionType(ResolutionType::EnergyElectron);
  (*fDetector)->RequestResolutionType(ResolutionType::EnergyMuon);
  (*fDetector)->RequestResolutionType(ResolutionType::MissingET);
}

// ---------------------------------------------------------
double KLFitter::BoostedLikelihoodTopLeptonJets::LogLikelihood(const std::vector<double> & parameters) {
  // calculate 4-vectors
  CalculateLorentzVectors(parameters);

  // define log of likelihood
  double logprob(0.);

  // temporary flag for a safe use of the transfer functions
  bool TFgoodTmp(true);

  // jet energy resolution terms
  logprob += log(fResEnergyBhad->p(bhad_fit_e, bhad_meas_e, &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;

  logprob += log(fResEnergyBlep->p(blep_fit_e, blep_meas_e, &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;

  logprob += log(fResEnergyLQ->p(lq_fit_e, lq_meas_e, &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;

  // lepton energy resolution terms
  if (fTypeLepton == kElectron) {
    logprob += log(fResLepton->p(lep_fit_e, lep_meas_e, &TFgoodTmp));
  } else if (fTypeLepton == kMuon) {
    logprob += log(fResLepton->p(lep_fit_e* lep_meas_sintheta, lep_meas_pt, &TFgoodTmp));
  }
  if (!TFgoodTmp) fTFgood = false;

  // neutrino px and py
  logprob += log(fResMET->p(nu_fit_px, ETmiss_x, &TFgoodTmp, SumET));
  if (!TFgoodTmp) fTFgood = false;

  logprob += log(fResMET->p(nu_fit_py, ETmiss_y, &TFgoodTmp, SumET));
  if (!TFgoodTmp) fTFgood = false;

  // physics constants
  double massW = fPhysicsConstants.MassW();
  double gammaW = fPhysicsConstants.GammaW();
  // note: top mass width should be made DEPENDENT on the top mass at a certain point
  //    fPhysicsConstants.SetMassTop(parameters[parTopM]);
  // (this will also set the correct width for the top)
  double gammaTop = fPhysicsConstants.GammaTop();

  // Breit-Wigner of leptonically decaying W-boson
  logprob += BCMath::LogBreitWignerRel(wlep_fit_m, massW, gammaW);

  // Breit-Wigner of hadronically decaying top quark
  logprob += BCMath::LogBreitWignerRel(thad_fit_m, parameters[parTopM], gammaTop);

  // Breit-Wigner of leptonically decaying top quark
  logprob += BCMath::LogBreitWignerRel(tlep_fit_m, parameters[parTopM], gammaTop);

  // return log of likelihood
  return logprob;
}

// ---------------------------------------------------------
std::vector<double> KLFitter::BoostedLikelihoodTopLeptonJets::GetInitialParameters() {
  std::vector<double> values(GetInitialParametersWoNeutrinoPz());

  // check second neutrino solution
  std::vector<double> neutrino_pz_solutions = GetNeutrinoPzSolutions();
  if (neutrino_pz_solutions.size() == 1) {
    values[parNuPz] = neutrino_pz_solutions[0];
  } else if (neutrino_pz_solutions.size() == 2) {
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
std::vector<double> KLFitter::BoostedLikelihoodTopLeptonJets::GetInitialParametersWoNeutrinoPz() {
  std::vector<double> values(GetNParameters());

  // energies of the quarks
  values[parBhadE] = bhad_meas_e;
  values[parBlepE] = blep_meas_e;
  values[parLQE]  = lq_meas_e;

  // energy of the lepton
  if (fTypeLepton == kElectron) {
    values[parLepE] = (*fParticlesPermuted)->Electron(0)->E();
  } else if (fTypeLepton == kMuon) {
    values[parLepE] = (*fParticlesPermuted)->Muon(0)->E();
  }

  // missing px and py
  values[parNuPx] = ETmiss_x;
  values[parNuPy] = ETmiss_y;

  // pz of the neutrino
  values[parNuPz] = 0.;

  // top mass
  double mtop = (*(*fParticlesPermuted)->Parton(0) + *(*fParticlesPermuted)->Parton(2)).M();
  if (mtop < GetParameter(parTopM)->GetLowerLimit()) {
    mtop = GetParameter(parTopM)->GetLowerLimit();
  } else if (mtop > GetParameter(parTopM)->GetUpperLimit()) {
    mtop = GetParameter(parTopM)->GetUpperLimit();
  }
  values[parTopM] = mtop;

  // return the vector
  return values;
}

// ---------------------------------------------------------
std::vector<double> KLFitter::BoostedLikelihoodTopLeptonJets::GetNeutrinoPzSolutions() {
  return CalculateNeutrinoPzSolutions();
}

// ---------------------------------------------------------
std::vector<double> KLFitter::BoostedLikelihoodTopLeptonJets::CalculateNeutrinoPzSolutions(TLorentzVector* additionalParticle) {
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
int KLFitter::BoostedLikelihoodTopLeptonJets::SavePermutedParticles() {
  bhad_meas_e      = (*fParticlesPermuted)->Parton(0)->E();
  bhad_meas_deteta = (*fParticlesPermuted)->DetEta(0, Particle::Type::kParton);
  bhad_meas_px     = (*fParticlesPermuted)->Parton(0)->Px();
  bhad_meas_py     = (*fParticlesPermuted)->Parton(0)->Py();
  bhad_meas_pz     = (*fParticlesPermuted)->Parton(0)->Pz();
  bhad_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(0)->M(), fPhysicsConstants.MassBottom(), &bhad_meas_px, &bhad_meas_py, &bhad_meas_pz, bhad_meas_e);
  bhad_meas_p      = sqrt(bhad_meas_e*bhad_meas_e - bhad_meas_m*bhad_meas_m);

  blep_meas_e      = (*fParticlesPermuted)->Parton(1)->E();
  blep_meas_deteta = (*fParticlesPermuted)->DetEta(1, Particle::Type::kParton);
  blep_meas_px     = (*fParticlesPermuted)->Parton(1)->Px();
  blep_meas_py     = (*fParticlesPermuted)->Parton(1)->Py();
  blep_meas_pz     = (*fParticlesPermuted)->Parton(1)->Pz();
  blep_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(1)->M(), fPhysicsConstants.MassBottom(), &blep_meas_px, &blep_meas_py, &blep_meas_pz, blep_meas_e);
  blep_meas_p      = sqrt(blep_meas_e*blep_meas_e - blep_meas_m*blep_meas_m);

  lq_meas_e      = (*fParticlesPermuted)->Parton(2)->E();
  lq_meas_deteta = (*fParticlesPermuted)->DetEta(2, Particle::Type::kParton);
  lq_meas_px     = (*fParticlesPermuted)->Parton(2)->Px();
  lq_meas_py     = (*fParticlesPermuted)->Parton(2)->Py();
  lq_meas_pz     = (*fParticlesPermuted)->Parton(2)->Pz();
  lq_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(2)->M(), 0., &lq_meas_px, &lq_meas_py, &lq_meas_pz, lq_meas_e);
  lq_meas_p      = sqrt(lq_meas_e*lq_meas_e - lq_meas_m*lq_meas_m);

  TLorentzVector * lepton(0);
  if (fTypeLepton == kElectron) {
    lepton = (*fParticlesPermuted)->Electron(0);
    lep_meas_deteta = (*fParticlesPermuted)->DetEta(0, Particle::Type::kElectron);
  } else {
    lepton = (*fParticlesPermuted)->Muon(0);
    lep_meas_deteta = (*fParticlesPermuted)->DetEta(0, Particle::Type::kMuon);
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
int KLFitter::BoostedLikelihoodTopLeptonJets::SaveResolutionFunctions() {
  fResEnergyBhad = (*fDetector)->ResEnergyBJet(bhad_meas_deteta);
  fResEnergyBlep = (*fDetector)->ResEnergyBJet(blep_meas_deteta);
  fResEnergyLQ  = (*fDetector)->ResEnergyLightJet(lq_meas_deteta);
  if (fTypeLepton == kElectron) {
    fResLepton = (*fDetector)->ResEnergyElectron(lep_meas_deteta);
  } else if (fTypeLepton == kMuon) {
    fResLepton = (*fDetector)->ResEnergyMuon(lep_meas_deteta);
  }
  fResMET = (*fDetector)->ResMissingET();

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::BoostedLikelihoodTopLeptonJets::BuildModelParticles() {
  if (GetBestFitParameters().size() > 0) CalculateLorentzVectors(GetBestFitParameters());

  TLorentzVector * bhad = fParticlesModel->Parton(0);
  TLorentzVector * blep = fParticlesModel->Parton(1);
  TLorentzVector * lq  = fParticlesModel->Parton(2);
  TLorentzVector * lep(0);
  if (fTypeLepton == kElectron) {
    lep  = fParticlesModel->Electron(0);
  } else if (fTypeLepton == kMuon) {
    lep  = fParticlesModel->Muon(0);
  }
  TLorentzVector * nu   = fParticlesModel->Neutrino(0);
  TLorentzVector * wlep  = fParticlesModel->Boson(0);
  TLorentzVector * thad  = fParticlesModel->Parton(3);
  TLorentzVector * tlep  = fParticlesModel->Parton(4);

  bhad->SetPxPyPzE(bhad_fit_px, bhad_fit_py, bhad_fit_pz, bhad_fit_e);
  blep->SetPxPyPzE(blep_fit_px, blep_fit_py, blep_fit_pz, blep_fit_e);
  lq ->SetPxPyPzE(lq_fit_px,  lq_fit_py,  lq_fit_pz,  lq_fit_e);
  lep ->SetPxPyPzE(lep_fit_px,  lep_fit_py,  lep_fit_pz,  lep_fit_e);
  nu  ->SetPxPyPzE(nu_fit_px,   nu_fit_py,   nu_fit_pz,   nu_fit_e);

  (*wlep) = (*lep)  + (*nu);
  (*thad) = (*lq) + (*bhad);
  (*tlep) = (*wlep) + (*blep);

  // no error
  return 1;
}

// ---------------------------------------------------------
std::vector<double> KLFitter::BoostedLikelihoodTopLeptonJets::LogLikelihoodComponents(std::vector<double> parameters) {
  std::vector<double> vecci;

  // calculate 4-vectors
  CalculateLorentzVectors(parameters);

  // temporary flag for a safe use of the transfer functions
  bool TFgoodTmp(true);

  // jet energy resolution terms
  vecci.push_back(log(fResEnergyBhad->p(bhad_fit_e, bhad_meas_e, &TFgoodTmp)));  // comp0
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log(fResEnergyBlep->p(blep_fit_e, blep_meas_e, &TFgoodTmp)));  // comp1
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log(fResEnergyLQ->p(lq_fit_e, lq_meas_e, &TFgoodTmp)));  // comp2
  if (!TFgoodTmp) fTFgood = false;

  // lepton energy resolution terms
  if (fTypeLepton == kElectron) {
    vecci.push_back(log(fResLepton->p(lep_fit_e, lep_meas_e, &TFgoodTmp)));  // comp3
  } else if (fTypeLepton == kMuon) {
    vecci.push_back(log(fResLepton->p(lep_fit_e* lep_meas_sintheta, lep_meas_pt, &TFgoodTmp)));  // comp3
  }
  if (!TFgoodTmp) fTFgood = false;

  // neutrino px and py
  vecci.push_back(log(fResMET->p(nu_fit_px, ETmiss_x, &TFgoodTmp, SumET)));  // comp4
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log(fResMET->p(nu_fit_py, ETmiss_y, &TFgoodTmp, SumET)));  // comp5
  if (!TFgoodTmp) fTFgood = false;

  // physics constants
  double massW = fPhysicsConstants.MassW();
  double gammaW = fPhysicsConstants.GammaW();
  // note: top mass width should be made DEPENDENT on the top mass at a certain point
  //    fPhysicsConstants.SetMassTop(parameters[parTopM]);
  // (this will also set the correct width for the top)
  double gammaTop = fPhysicsConstants.GammaTop();

  // Breit-Wigner of leptonically decaying W-boson
  vecci.push_back(BCMath::LogBreitWignerRel(wlep_fit_m, massW, gammaW));  // comp6

  // Breit-Wigner of hadronically decaying top quark
  vecci.push_back(BCMath::LogBreitWignerRel(thad_fit_m, parameters[parTopM], gammaTop));  // comp7

  // Breit-Wigner of leptonically decaying top quark
  vecci.push_back(BCMath::LogBreitWignerRel(tlep_fit_m, parameters[parTopM], gammaTop));  // comp8

  // return log of likelihood
  return vecci;
}
