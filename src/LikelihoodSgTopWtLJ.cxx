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

#include "KLFitter/LikelihoodSgTopWtLJ.h"

#include <cmath>
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
KLFitter::LikelihoodSgTopWtLJ::LikelihoodSgTopWtLJ()
  : KLFitter::LikelihoodBase::LikelihoodBase()
  , fHadronicTop(true)
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
KLFitter::LikelihoodSgTopWtLJ::~LikelihoodSgTopWtLJ() = default;

// ---------------------------------------------------------
TLorentzVector* KLFitter::LikelihoodSgTopWtLJ::GetLepton(KLFitter::ParticleCollection* particles) {
  if (!particles) {
    std::cout << "KLFitter::LikelihoodSgTopWtLJ::GetLepton():\tERROR\t Null pointer to particles object." << std::endl;
    return 0;
  }
  TLorentzVector* lepton = 0;
  if (fTypeLepton == kElectron) {
    lepton = particles->GetP4(Particles::Type::kElectron, 0);
  } else if (fTypeLepton == kMuon) {
    lepton = particles->GetP4(Particles::Type::kMuon, 0);
  } else {
    std::cout << "KLFitter::LikelihoodSgTopWtLJ::GetLepton():\tERROR\tInvalid letpon type: " << fTypeLepton << std::endl;
  }

  return lepton;
}

// ---------------------------------------------------------
int KLFitter::LikelihoodSgTopWtLJ::GetLeptonType() {
  return fTypeLepton;
}

// ---------------------------------------------------------
int KLFitter::LikelihoodSgTopWtLJ::SetET_miss_XY_SumET(double etx, double ety, double sumet) {
  // set missing ET x and y component and the SumET
  ETmiss_x = etx;
  ETmiss_y = ety;
  SumET = sumet;

  // no error
  return 1;
}

// ---------------------------------------------------------
void KLFitter::LikelihoodSgTopWtLJ::RequestResolutionFunctions() {
  (*fDetector)->RequestResolutionType(ResolutionType::EnergyLightJet);
  (*fDetector)->RequestResolutionType(ResolutionType::EnergyBJet);
  (*fDetector)->RequestResolutionType(ResolutionType::EnergyElectron);
  (*fDetector)->RequestResolutionType(ResolutionType::EnergyMuon);
  (*fDetector)->RequestResolutionType(ResolutionType::MissingET);
}

// ---------------------------------------------------------
void KLFitter::LikelihoodSgTopWtLJ::SetLeptonType(int leptontype) {
  if (leptontype != kElectron && leptontype != kMuon) {
    std::cout << "KLFitter::SetLeptonType()\tWARNING\t lepton type not defined: " << leptontype << ". Set electron as lepton type." << std::endl;
    fTypeLepton = kElectron;
  } else {
    fTypeLepton = leptontype;
  }

  // define model particles
  DefineModelParticles();
}

// ---------------------------------------------------------
int KLFitter::LikelihoodSgTopWtLJ::DefineModelParticles() {
  // create the particles of the model
  fParticlesModel.reset(new KLFitter::ParticleCollection{});

  Particles::Parton parton0{"b quark", TLorentzVector{}};
  parton0.SetIdentifier(0);
  parton0.SetTrueFlavor(Particles::PartonTrueFlavor::kB);
  fParticlesModel->AddParticle(parton0);

  Particles::Parton parton1{"light quark 1", TLorentzVector{}};
  parton1.SetIdentifier(1);
  parton1.SetTrueFlavor(Particles::PartonTrueFlavor::kLight);
  fParticlesModel->AddParticle(parton1);

  Particles::Parton parton2{"light quark 2", TLorentzVector{}};
  parton2.SetIdentifier(2);
  parton2.SetTrueFlavor(Particles::PartonTrueFlavor::kLight);
  fParticlesModel->AddParticle(parton2);

  if (fTypeLepton == kElectron) {
    fParticlesModel->AddParticle(Particles::Electron{"electron", TLorentzVector{}});
  } else if (fTypeLepton == kMuon) {
    fParticlesModel->AddParticle(Particles::Muon{"muon", TLorentzVector{}});
  }

  fParticlesModel->AddParticle(Particles::Neutrino{"neutrino", TLorentzVector{}});
  fParticlesModel->AddParticle(Particles::Boson{"hadronic W", TLorentzVector{}});
  fParticlesModel->AddParticle(Particles::Boson{"leptonic W", TLorentzVector{}});
  fParticlesModel->AddParticle(Particles::Parton{"hadronic top", TLorentzVector{}});
  fParticlesModel->AddParticle(Particles::Parton{"leptonic top", TLorentzVector{}});

  // no error
  return 1;
}

// ---------------------------------------------------------
void KLFitter::LikelihoodSgTopWtLJ::DefineParameters() {
  // add parameters of model
  this->AddParameter("energy b",                fPhysicsConstants.MassBottom(), 1000.0);   // parBE
  this->AddParameter("energy light quark 1",    0.0, 1000.0);                               // parLQ1E
  this->AddParameter("energy light quark 2",    0.0, 1000.0);                               // parLQ2E
  this->AddParameter("energy lepton",           0.0, 1000.0);                               // parLepE
  this->AddParameter("missPx",              -1000.0, 1000.0);                               // parNuPx
  this->AddParameter("missPy",              -1000.0, 1000.0);                               // parNuPy
  this->AddParameter("p_z neutrino",        -1000.0, 1000.0);                               // parNuPz
}

// ---------------------------------------------------------
int KLFitter::LikelihoodSgTopWtLJ::CalculateLorentzVectors(const std::vector <double>& parameters) {
  // variables
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

  // b quark
  b_fit_e = parameters[parBE];
  scale = sqrt(b_fit_e*b_fit_e - b_meas_m*b_meas_m) / b_meas_p;
  b_fit_px = scale * b_meas_px;
  b_fit_py = scale * b_meas_py;
  b_fit_pz = scale * b_meas_pz;

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

  // composite particles

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
  thad_fit_e  = whad_fit_e+b_fit_e;
  thad_fit_px = whad_fit_px+b_fit_px;
  thad_fit_py = whad_fit_py+b_fit_py;
  thad_fit_pz = whad_fit_pz+b_fit_pz;
  thad_fit_m = sqrt(thad_fit_e*thad_fit_e - (thad_fit_px*thad_fit_px + thad_fit_py*thad_fit_py + thad_fit_pz*thad_fit_pz));

  // leptonic top
  tlep_fit_e = wlep_fit_e+b_fit_e;
  tlep_fit_px = wlep_fit_px+b_fit_px;
  tlep_fit_py = wlep_fit_py+b_fit_py;
  tlep_fit_pz = wlep_fit_pz+b_fit_pz;
  tlep_fit_m = sqrt(tlep_fit_e*tlep_fit_e - (tlep_fit_px*tlep_fit_px + tlep_fit_py*tlep_fit_py + tlep_fit_pz*tlep_fit_pz));

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::LikelihoodSgTopWtLJ::Initialize() {
  // error code
  int err = 1;

  ResetResults();

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
int KLFitter::LikelihoodSgTopWtLJ::RemoveInvariantParticlePermutations() {
  // error code
  int err = 1;

  // remove the permutation from the first and the second jet
  Particles::Type ptype = Particles::Type::kParton;
  std::vector<int> indexVector_Jets;
  indexVector_Jets.push_back(1);
  indexVector_Jets.push_back(2);
  err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Jets);

  // remove the permutation from all jet not used for the likelihood
  const KLFitter::ParticleCollection* particles = (*fPermutations)->Particles();
  indexVector_Jets.clear();
  for (size_t i = 3; i < particles->partons.size(); i++) {
    indexVector_Jets.push_back(i);
  }
  err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Jets);

  // return error code
  return err;
}

// ---------------------------------------------------------
int KLFitter::LikelihoodSgTopWtLJ::AdjustParameterRanges() {
  // adjust limits
  double nsigmas_jet = 7.0;
  double nsigmas_lepton = 2.0;

  // energy of b quark
  double E = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 0)->E();
  double m = fPhysicsConstants.MassBottom();
  if (fFlagUseJetMass)
    m = TMath::Max(0.0, (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 0)->M());
  double Emin = TMath::Max(m, E - nsigmas_jet * sqrt(E));
  double Emax  = E + nsigmas_jet * sqrt(E);
  SetParameterRange(parBE, Emin, Emax);

  // energy of light quark 1
  E = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 1)->E();
  m = 0.001;
  if (fFlagUseJetMass)
    m = TMath::Max(0.0, (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 1)->M());
  Emin = TMath::Max(m, E - nsigmas_jet * sqrt(E));
  Emax  = E + nsigmas_jet * sqrt(E);
  SetParameterRange(parLQ1E, Emin, Emax);

  // energy of light quark2
  E = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 2)->E();
  m = 0.001;
  if (fFlagUseJetMass)
    m = TMath::Max(0.0, (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 2)->M());
  Emin = TMath::Max(m, E - nsigmas_jet * sqrt(E));
  Emax  = E + nsigmas_jet * sqrt(E);
  SetParameterRange(parLQ2E, Emin, Emax);

  // energy of lepton
  if (fTypeLepton == kElectron) {
    E = (*fParticlesPermuted)->GetP4(Particles::Type::kElectron, 0)->E();
    Emin = TMath::Max(0.001, E - nsigmas_lepton * sqrt(E));
    Emax = E + nsigmas_lepton * sqrt(E);
  } else if (fTypeLepton == kMuon) {
    E = (*fParticlesPermuted)->GetP4(Particles::Type::kMuon, 0)->E();
    double sintheta = sin((*fParticlesPermuted)->GetP4(Particles::Type::kMuon, 0)->Theta());
    double sigrange = nsigmas_lepton * (E * E * sintheta);
    Emin = std::max(0.001, E - sigrange);
    Emax = E + sigrange;
  }
  SetParameterRange(parLepE, Emin, Emax);

  // note: this is hard-coded in the momement

  // missing ET
  SetParameterRange(parNuPx, ETmiss_x-100.0, ETmiss_x+100);
  SetParameterRange(parNuPy, ETmiss_y-100.0, ETmiss_y+100);

  // no error
  return 1;
}

// ---------------------------------------------------------
double KLFitter::LikelihoodSgTopWtLJ::LogLikelihood(const std::vector<double> & parameters) {
  // calculate 4-vectors
  this->CalculateLorentzVectors(parameters);

  // define log of likelihood
  double logprob(0.);

  // temporary flag for a safe use of the transfer functions
  bool TFgoodTmp(true);

  // jet energy resolution terms
  logprob += fResEnergyB->logp(b_fit_e, b_meas_e, &TFgoodTmp);
  if (!TFgoodTmp) fTFgood = false;
  logprob += fResEnergyLQ1->logp(lq1_fit_e, lq1_meas_e, &TFgoodTmp);
  if (!TFgoodTmp) fTFgood = false;
  logprob += fResEnergyLQ2->logp(lq2_fit_e, lq2_meas_e, &TFgoodTmp);
  if (!TFgoodTmp) fTFgood = false;

  // lepton energy resolution terms
  if (fTypeLepton == kElectron) {
    logprob += fResLepton->logp(lep_fit_e, lep_meas_e, &TFgoodTmp);
  } else if (fTypeLepton == kMuon) {
    logprob += fResLepton->logp(lep_fit_e* lep_meas_sintheta, lep_meas_pt, &TFgoodTmp);
  }
  if (!TFgoodTmp) fTFgood = false;

  // neutrino px and py
  logprob += fResMET->logp(nu_fit_px, ETmiss_x, &TFgoodTmp, SumET);
  if (!TFgoodTmp) fTFgood = false;
  logprob += fResMET->logp(nu_fit_py, ETmiss_y, &TFgoodTmp, SumET);
  if (!TFgoodTmp) fTFgood = false;

  // physics constants
  double massW = fPhysicsConstants.MassW();
  double gammaW = fPhysicsConstants.GammaW();
  double massTop = fPhysicsConstants.MassTop();
  double gammaTop = fPhysicsConstants.GammaTop();

  // Breit-Wigner of hadronically decaying W-boson
  logprob += BCMath::LogBreitWignerRel(whad_fit_m, massW, gammaW);

  // Breit-Wigner of leptonically decaying W-boson
  logprob += BCMath::LogBreitWignerRel(wlep_fit_m, massW, gammaW);

  if (fHadronicTop) {
    logprob += BCMath::LogBreitWignerRel(thad_fit_m, massTop, gammaTop);
  } else {
    logprob += BCMath::LogBreitWignerRel(tlep_fit_m, massTop, gammaTop);
  }

  // return log of likelihood
  return logprob;
}

// ---------------------------------------------------------
std::vector<double> KLFitter::LikelihoodSgTopWtLJ::GetInitialParameters() {
  std::vector<double> values(GetNParameters());

  // energies of the quarks
  values[parBE] = b_meas_e;
  values[parLQ1E]  = lq1_meas_e;
  values[parLQ2E]  = lq2_meas_e;

  // energy of the lepton
  values[parLepE] = lep_meas_e;

  // missing px and py
  values[parNuPx] = ETmiss_x;
  values[parNuPy] = ETmiss_y;

  // pz of the neutrino
  values[parNuPz] = 0;

  // check neutrino solutions
  std::vector<double> neutrino_pz_solutions = GetNeutrinoPzSolutions();
  if (neutrino_pz_solutions.size() == 1) {
    values[parNuPz] = neutrino_pz_solutions.at(0);
  } else if (neutrino_pz_solutions.size() == 2) {
    if (fHadronicTop) {
      if (std::fabs(neutrino_pz_solutions[0]) < std::fabs(neutrino_pz_solutions[1])) {
        values[parNuPz] = neutrino_pz_solutions[0];
      } else {
        values[parNuPz] = neutrino_pz_solutions[1];
      }
      SetParameterRange(parNuPx, values[parNuPx], values[parNuPx]);
      SetParameterRange(parNuPy, values[parNuPy], values[parNuPy]);
      SetParameterRange(parNuPz, values[parNuPz], values[parNuPz]);
      SetParameterRange(parLepE, values[parLepE], values[parLepE]);
    } else {
      double sol1, sol2;
      values[parNuPz] = neutrino_pz_solutions.at(0);
      sol1 = LogLikelihood(values);
      values[parNuPz] = neutrino_pz_solutions.at(1);
      sol2 = LogLikelihood(values);
      if (sol1 > sol2) {
        values[parNuPz] = neutrino_pz_solutions.at(0);
      }
    }
  } else {
    // scaleMET method
    TLorentzVector* lepton = GetLepton(*fParticlesPermuted);
    KLFitter::PhysicsConstants constants;
    double mW = constants.MassW();
    double cosDeltaPhi = (lepton->Px()*ETmiss_x + lepton->Py()*ETmiss_y) / (lepton->Perp() * sqrt(ETmiss_x * ETmiss_x + ETmiss_y * ETmiss_y));
    values[parNuPz] = mW * mW * lepton->Pz() / (2 * lepton->Perp2()*(1 - cosDeltaPhi));
    if (fHadronicTop)
      SetParameterRange(parNuPz, values[parNuPz], values[parNuPz]);
  }

  // return the vector
  return values;
}

// ---------------------------------------------------------
std::vector<double> KLFitter::LikelihoodSgTopWtLJ::GetNeutrinoPzSolutions() {
  std::vector<double> pz;

  KLFitter::PhysicsConstants constants;
  // electron mass
  double mE = 0.;

  TLorentzVector* lepton = GetLepton(*fParticlesPermuted);

  double px_c = lepton->Px();
  double py_c = lepton->Py();
  double pz_c = lepton->Pz();
  double Ec = lepton->E();

  double px_nu = ETmiss_x;
  double py_nu = ETmiss_y;
  double alpha = constants.MassW()*constants.MassW() - mE*mE + 2*(px_c*px_nu + py_c*py_nu);

  double a = pz_c*pz_c - Ec*Ec;
  double b = alpha * pz_c;
  double c = - Ec*Ec * (px_nu*px_nu + py_nu*py_nu) + alpha*alpha/4.;

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
int KLFitter::LikelihoodSgTopWtLJ::SavePermutedParticles() {
  b_meas_e      = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 0)->E();
  b_meas_deteta = (*fParticlesPermuted)->partons.at(0).GetDetEta();
  b_meas_px     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 0)->Px();
  b_meas_py     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 0)->Py();
  b_meas_pz     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 0)->Pz();
  b_meas_m      = SetPartonMass((*fParticlesPermuted)->GetP4(Particles::Type::kParton, 0)->M(), fPhysicsConstants.MassBottom(), &b_meas_px, &b_meas_py, &b_meas_pz, b_meas_e);
  b_meas_p      = sqrt(b_meas_e*b_meas_e - b_meas_m*b_meas_m);

  lq1_meas_e      = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 1)->E();
  lq1_meas_deteta = (*fParticlesPermuted)->partons.at(1).GetDetEta();
  lq1_meas_px     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 1)->Px();
  lq1_meas_py     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 1)->Py();
  lq1_meas_pz     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 1)->Pz();
  lq1_meas_m      = SetPartonMass((*fParticlesPermuted)->GetP4(Particles::Type::kParton, 1)->M(), 0., &lq1_meas_px, &lq1_meas_py, &lq1_meas_pz, lq1_meas_e);
  lq1_meas_p      = sqrt(lq1_meas_e*lq1_meas_e - lq1_meas_m*lq1_meas_m);

  lq2_meas_e      = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 2)->E();
  lq2_meas_deteta = (*fParticlesPermuted)->partons.at(2).GetDetEta();
  lq2_meas_px     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 2)->Px();
  lq2_meas_py     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 2)->Py();
  lq2_meas_pz     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 2)->Pz();
  lq2_meas_m      = SetPartonMass((*fParticlesPermuted)->GetP4(Particles::Type::kParton, 2)->M(), 0., &lq2_meas_px, &lq2_meas_py, &lq2_meas_pz, lq2_meas_e);
  lq2_meas_p      = sqrt(lq2_meas_e*lq2_meas_e - lq2_meas_m*lq2_meas_m);

  TLorentzVector * lepton(0);
  if (fTypeLepton == kElectron) {
    lepton = (*fParticlesPermuted)->GetP4(Particles::Type::kElectron, 0);
    lep_meas_deteta = (*fParticlesPermuted)->electrons.at(0).GetDetEta();
  } else {
    lepton = (*fParticlesPermuted)->GetP4(Particles::Type::kMuon, 0);
    lep_meas_deteta = (*fParticlesPermuted)->muons.at(0).GetDetEta();
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
int KLFitter::LikelihoodSgTopWtLJ::SaveResolutionFunctions() {
  fResEnergyB = (*fDetector)->ResEnergyBJet(b_meas_deteta);
  fResEnergyLQ1  = (*fDetector)->ResEnergyLightJet(lq1_meas_deteta);
  fResEnergyLQ2  = (*fDetector)->ResEnergyLightJet(lq2_meas_deteta);
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
int KLFitter::LikelihoodSgTopWtLJ::BuildModelParticles() {
  if (GetBestFitParameters().size() > 0) CalculateLorentzVectors(GetBestFitParameters());

  TLorentzVector * b = fParticlesModel->GetP4(Particles::Type::kParton, 0);
  TLorentzVector * lq1  = fParticlesModel->GetP4(Particles::Type::kParton, 1);
  TLorentzVector * lq2  = fParticlesModel->GetP4(Particles::Type::kParton, 2);
  TLorentzVector * lep = GetLepton(fParticlesModel.get());
  TLorentzVector * nu   = fParticlesModel->GetP4(Particles::Type::kNeutrino, 0);
  TLorentzVector * whad  = fParticlesModel->GetP4(Particles::Type::kBoson, 0);
  TLorentzVector * wlep  = fParticlesModel->GetP4(Particles::Type::kBoson, 1);
  TLorentzVector * thad  = fParticlesModel->GetP4(Particles::Type::kParton, 3);
  TLorentzVector * tlep  = fParticlesModel->GetP4(Particles::Type::kParton, 4);

  b->SetPxPyPzE(b_fit_px, b_fit_py, b_fit_pz, b_fit_e);
  lq1 ->SetPxPyPzE(lq1_fit_px,  lq1_fit_py,  lq1_fit_pz,  lq1_fit_e);
  lq2 ->SetPxPyPzE(lq2_fit_px,  lq2_fit_py,  lq2_fit_pz,  lq2_fit_e);
  lep ->SetPxPyPzE(lep_fit_px,  lep_fit_py,  lep_fit_pz,  lep_fit_e);
  nu  ->SetPxPyPzE(nu_fit_px,   nu_fit_py,   nu_fit_pz,   nu_fit_e);

  (*whad) = (*lq1)  + (*lq2);
  (*wlep) = (*lep)  + (*nu);
  (*thad) = (*whad) + (*b);
  (*tlep) = (*wlep) + (*b);

  // no error
  return 1;
}
