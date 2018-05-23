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

#include "KLFitter/LikelihoodTopLeptonJets_JetAngles.h"

#include <algorithm>
#include <iostream>

#include "BAT/BCMath.h"
#include "BAT/BCParameter.h"
#include "KLFitter/DetectorBase.h"
#include "KLFitter/Particles.h"
#include "KLFitter/Permutations.h"
#include "KLFitter/PhysicsConstants.h"
#include "KLFitter/ResolutionBase.h"
#include "TLorentzVector.h"

// ---------------------------------------------------------
KLFitter::LikelihoodTopLeptonJets_JetAngles::LikelihoodTopLeptonJets_JetAngles() = default;

// ---------------------------------------------------------
KLFitter::LikelihoodTopLeptonJets_JetAngles::~LikelihoodTopLeptonJets_JetAngles() = default;

// ---------------------------------------------------------
void KLFitter::LikelihoodTopLeptonJets_JetAngles::DefineParameters() {
  // add parameters of model
  AddParameter("energy hadronic b",       fPhysicsConstants.MassBottom(), 1000.0);  // parBhadE
  AddParameter("energy leptonic b",       fPhysicsConstants.MassBottom(), 1000.0);  // parBlepE
  AddParameter("energy light quark 1",    0.0, 1000.0);                              // parLQ1E
  AddParameter("energy light quark 2",    0.0, 1000.0);                              // parLQ2E
  AddParameter("energy lepton",           0.0, 1000.0);                              // parLepE
  AddParameter("p_x neutrino",        -1000.0, 1000.0);                              // parNuPx
  AddParameter("p_y neutrino",        -1000.0, 1000.0);                              // parNuPy
  AddParameter("p_z neutrino",        -1000.0, 1000.0);                              // parNuPz
  AddParameter("eta hadronic b",       -2.5, 2.5);                                   // parBhadEta
  AddParameter("eta leptonic b",       -2.5, 2.5);                                   // parBlepEta
  AddParameter("eta light quark 1",    -2.5, 2.5);                                   // parLQ1Eta
  AddParameter("eta light quark 2",    -2.5, 2.5);                                   // parLQ2Eta
  AddParameter("phi hadronic b",       -TMath::Pi(), TMath::Pi());                   // parBhadPhi
  AddParameter("phi leptonic b",       -TMath::Pi(), TMath::Pi());                   // parBlepPhi
  AddParameter("phi light quark 1",    -TMath::Pi(), TMath::Pi());                   // parLQ1Phi
  AddParameter("phi light quark 2",    -TMath::Pi(), TMath::Pi());                   // parLQ2Phi
  AddParameter("top mass",              100.0, 1000.0);                              // parTopM
}

// ---------------------------------------------------------
int KLFitter::LikelihoodTopLeptonJets_JetAngles::CalculateLorentzVectors(std::vector <double> const& parameters) {
  static double scale;
  static double thad_fit_e;
  static double thad_fit_px;
  static double thad_fit_py;
  static double thad_fit_pz;
  static double tlep_fit_e;
  static double tlep_fit_px;
  static double tlep_fit_py;
  static double tlep_fit_pz;

  static TLorentzVector v;

  // hadronic b quark
  v.SetPtEtaPhiE(sqrt(parameters[parBhadE]*parameters[parBhadE]-bhad_meas_m*bhad_meas_m)/cosh(parameters[parBhadEta]), parameters[parBhadEta], parameters[parBhadPhi], parameters[parBhadE]);
  bhad_fit_e = v.E();
  bhad_fit_px = v.Px();
  bhad_fit_py = v.Py();
  bhad_fit_pz = v.Pz();

  // leptonic b quark
  v.SetPtEtaPhiE(sqrt(parameters[parBlepE]*parameters[parBlepE]-blep_meas_m*blep_meas_m)/cosh(parameters[parBlepEta]), parameters[parBlepEta], parameters[parBlepPhi], parameters[parBlepE]);
  blep_fit_e = v.E();
  blep_fit_px = v.Px();
  blep_fit_py = v.Py();
  blep_fit_pz = v.Pz();

  // light quark 1
  v.SetPtEtaPhiE(sqrt(parameters[parLQ1E]*parameters[parLQ1E]-lq1_meas_m*lq1_meas_m)/cosh(parameters[parLQ1Eta]), parameters[parLQ1Eta], parameters[parLQ1Phi], parameters[parLQ1E]);
  lq1_fit_e = v.E();
  lq1_fit_px = v.Px();
  lq1_fit_py = v.Py();
  lq1_fit_pz = v.Pz();

  // light quark 2
  v.SetPtEtaPhiE(sqrt(parameters[parLQ2E]*parameters[parLQ2E]-lq2_meas_m*lq2_meas_m)/cosh(parameters[parLQ2Eta]), parameters[parLQ2Eta], parameters[parLQ2Phi], parameters[parLQ2E]);
  lq2_fit_e = v.E();
  lq2_fit_px = v.Px();
  lq2_fit_py = v.Py();
  lq2_fit_pz = v.Pz();

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

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::LikelihoodTopLeptonJets_JetAngles::AdjustParameterRanges() {
  // adjust limits
  double nsigmas_jet    = fFlagGetParSigmasFromTFs ? 10 : 7;
  double nsigmas_lepton = fFlagGetParSigmasFromTFs ? 10 : 2;
  double nsigmas_met    = fFlagGetParSigmasFromTFs ? 10 : 1;

  double E = (*fParticlesPermuted)->Parton(0)->E();
  double m = fPhysicsConstants.MassBottom();
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(0)->M());
  double sigma = fFlagGetParSigmasFromTFs ? (*fDetector)->ResEnergyBJet((*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kParton))->GetSigma(E) : sqrt(E);
  double Emin = std::max(m, E - nsigmas_jet* sigma);
  double Emax  = E + nsigmas_jet* sigma;
  SetParameterRange(parBhadE, Emin, Emax);

  E = (*fParticlesPermuted)->Parton(1)->E();
  m = fPhysicsConstants.MassBottom();
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(1)->M());
  sigma = fFlagGetParSigmasFromTFs ? (*fDetector)->ResEnergyBJet((*fParticlesPermuted)->DetEta(1, KLFitter::Particles::kParton))->GetSigma(E) : sqrt(E);
  Emin = std::max(m, E - nsigmas_jet* sigma);
  Emax  = E + nsigmas_jet* sigma;
  SetParameterRange(parBlepE, Emin, Emax);

  E = (*fParticlesPermuted)->Parton(2)->E();
  m = 0.001;
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(2)->M());
  sigma = fFlagGetParSigmasFromTFs ? (*fDetector)->ResEnergyLightJet((*fParticlesPermuted)->DetEta(2, KLFitter::Particles::kParton))->GetSigma(E) : sqrt(E);
  Emin = std::max(m, E - nsigmas_jet* sigma);
  Emax  = E + nsigmas_jet* sigma;
  SetParameterRange(parLQ1E, Emin, Emax);

  E = (*fParticlesPermuted)->Parton(3)->E();
  m = 0.001;
  if (fFlagUseJetMass)
    m = std::max(0.0, (*fParticlesPermuted)->Parton(3)->M());
  sigma = fFlagGetParSigmasFromTFs ? (*fDetector)->ResEnergyLightJet((*fParticlesPermuted)->DetEta(3, KLFitter::Particles::kParton))->GetSigma(E) : sqrt(E);
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
    Emin = std::max(0.001, E -sigrange);
    Emax = E +sigrange;
  }
  SetParameterRange(parLepE, Emin, Emax);

  // note: this is hard-coded in the momement

  sigma = fFlagGetParSigmasFromTFs ? fResMET->GetSigma(SumET) : 100;
  double sigrange = nsigmas_met*sigma;
  SetParameterRange(parNuPx, ETmiss_x-sigrange, ETmiss_x+sigrange);
  SetParameterRange(parNuPy, ETmiss_y-sigrange, ETmiss_y+sigrange);

  // eta
  double eta = (*fParticlesPermuted)->Parton(0)->Eta();
  double etamin = std::max(-2.5, eta - 0.2);
  double etamax = std::min(2.5, eta + 0.2);
  SetParameterRange(parBhadEta, etamin, etamax);

  eta = (*fParticlesPermuted)->Parton(1)->Eta();
  etamin = std::max(-2.5, eta - 0.2);
  etamax = std::min(2.5, eta + 0.2);
  SetParameterRange(parBlepEta, etamin, etamax);

  eta = (*fParticlesPermuted)->Parton(2)->Eta();
  etamin = std::max(-2.5, eta - 0.2);
  etamax = std::min(2.5, eta + 0.2);
  SetParameterRange(parLQ1Eta, etamin, etamax);

  eta = (*fParticlesPermuted)->Parton(3)->Eta();
  etamin = std::max(-2.5, eta - 0.2);
  etamax = std::min(2.5, eta + 0.2);
  SetParameterRange(parLQ2Eta, etamin, etamax);

  // phi
  double phi = (*fParticlesPermuted)->Parton(0)->Phi();
  double phimin = phi - 0.1;
  double phimax = phi + 0.1;
  SetParameterRange(parBhadPhi, phimin, phimax);

  phi = (*fParticlesPermuted)->Parton(1)->Phi();
  phimin = phi - 0.1;
  phimax = phi + 0.1;
  SetParameterRange(parBlepPhi, phimin, phimax);

  phi = (*fParticlesPermuted)->Parton(2)->Phi();
  phimin = phi - 0.1;
  phimax = phi + 0.1;
  SetParameterRange(parLQ1Phi, phimin, phimax);

  phi = (*fParticlesPermuted)->Parton(3)->Phi();
  phimin = phi - 0.1;
  phimax = phi + 0.1;
  SetParameterRange(parLQ2Phi, phimin, phimax);

  if (fFlagTopMassFixed)
    SetParameterRange(parTopM, fPhysicsConstants.MassTop(), fPhysicsConstants.MassTop());

  // no error
  return 1;
}

// ---------------------------------------------------------
double KLFitter::LikelihoodTopLeptonJets_JetAngles::LogLikelihood(const std::vector<double> & parameters) {
  // calculate 4-vectors
  CalculateLorentzVectors(parameters);

  // define log of likelihood
  double logprob(0.);

  // temporary flag for a safe use of the transfer functions
  bool TFgoodTmp(true);

  // jet energy resolution terms
  logprob += log((*fDetector)->ResEnergyBJet((*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kParton))->p(bhad_fit_e, bhad_meas_e, &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;

  logprob += log((*fDetector)->ResEnergyBJet((*fParticlesPermuted)->DetEta(1, KLFitter::Particles::kParton))->p(blep_fit_e, blep_meas_e, &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;

  logprob += log((*fDetector)->ResEnergyLightJet((*fParticlesPermuted)->DetEta(2, KLFitter::Particles::kParton))->p(lq1_fit_e, lq1_meas_e, &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;

  logprob += log((*fDetector)->ResEnergyLightJet((*fParticlesPermuted)->DetEta(3, KLFitter::Particles::kParton))->p(lq2_fit_e, lq2_meas_e, &TFgoodTmp));
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

  // eta resolution
  logprob += log((*fDetector)->ResEtaBJet((*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kParton))->p(parameters[parBhadEta], (*fParticlesPermuted)->Parton(0)->Eta(), &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;
  logprob += log((*fDetector)->ResEtaBJet((*fParticlesPermuted)->DetEta(1, KLFitter::Particles::kParton))->p(parameters[parBlepEta], (*fParticlesPermuted)->Parton(1)->Eta(), &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;
  logprob += log((*fDetector)->ResEtaLightJet((*fParticlesPermuted)->DetEta(2, KLFitter::Particles::kParton))->p(parameters[parLQ1Eta], (*fParticlesPermuted)->Parton(2)->Eta(), &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;
  logprob += log((*fDetector)->ResEtaLightJet((*fParticlesPermuted)->DetEta(3, KLFitter::Particles::kParton))->p(parameters[parLQ2Eta], (*fParticlesPermuted)->Parton(3)->Eta(), &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;

  // transform all phi values, so that they are centered around zero, and not around the measured phi

  // phi resolution
  logprob += log((*fDetector)->ResPhiBJet((*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kParton))->p(diffPhi(parameters[parBhadPhi], (*fParticlesPermuted)->Parton(0)->Phi()), 0., &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;
  logprob += log((*fDetector)->ResPhiBJet((*fParticlesPermuted)->DetEta(1, KLFitter::Particles::kParton))->p(diffPhi(parameters[parBlepPhi], (*fParticlesPermuted)->Parton(1)->Phi()), 0., &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;
  logprob += log((*fDetector)->ResPhiLightJet((*fParticlesPermuted)->DetEta(2, KLFitter::Particles::kParton))->p(diffPhi(parameters[parLQ1Phi], (*fParticlesPermuted)->Parton(2)->Phi()), 0., &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;
  logprob += log((*fDetector)->ResPhiLightJet((*fParticlesPermuted)->DetEta(3, KLFitter::Particles::kParton))->p(diffPhi(parameters[parLQ2Phi], (*fParticlesPermuted)->Parton(3)->Phi()), 0., &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;

  // physics constants
  double massW = fPhysicsConstants.MassW();
  double gammaW = fPhysicsConstants.GammaW();
  // note: top mass width should be made DEPENDENT on the top mass at a certain point
  //    fPhysicsConstants.SetMassTop(parameters[parTopM]);
  // (this will also set the correct width for the top)
  double gammaTop = fPhysicsConstants.GammaTop();

  // Breit-Wigner of hadronically decaying W-boson
  logprob += BCMath::LogBreitWignerRel(whad_fit_m, massW, gammaW);

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
std::vector<double> KLFitter::LikelihoodTopLeptonJets_JetAngles::GetInitialParametersWoNeutrinoPz() {
  std::vector<double> values(GetNParameters());

  // energies of the quarks
  values[parBhadE] = bhad_meas_e;
  values[parBlepE] = blep_meas_e;
  values[parLQ1E]  = lq1_meas_e;
  values[parLQ2E]  = lq2_meas_e;

  // energy of the lepton
  if (fTypeLepton == kElectron) {
    values[parLepE] = (*fParticlesPermuted)->Electron(0)->E();
  } else if (fTypeLepton == kMuon) {
    values[parLepE] = (*fParticlesPermuted)->Muon(0)->E();
  }

  // missing px and py
  values[parNuPx] = ETmiss_x;
  values[parNuPy] = ETmiss_y;

  // eta of the quarks
  values[parBhadEta] = bhad_meas_eta;
  values[parBlepEta] = blep_meas_eta;
  values[parLQ1Eta]  = lq1_meas_eta;
  values[parLQ2Eta]  = lq2_meas_eta;

  // phi of the quarks
  values[parBhadPhi] = bhad_meas_phi;
  values[parBlepPhi] = blep_meas_phi;
  values[parLQ1Phi]  = lq1_meas_phi;
  values[parLQ2Phi]  = lq2_meas_phi;

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

  // return the vector
  return values;
}

// ---------------------------------------------------------
int KLFitter::LikelihoodTopLeptonJets_JetAngles::SaveResolutionFunctions() {
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
std::vector<double> KLFitter::LikelihoodTopLeptonJets_JetAngles::LogLikelihoodComponents(std::vector<double> parameters) {
  std::vector<double> vecci;

  // calculate 4-vectors
  CalculateLorentzVectors(parameters);

  // temporary flag for a safe use of the transfer functions
  bool TFgoodTmp(true);

  // jet energy resolution terms
  vecci.push_back(log((*fDetector)->ResEnergyBJet((*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kParton))->p(bhad_fit_e, bhad_meas_e, &TFgoodTmp)));  // comp0
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log((*fDetector)->ResEnergyBJet((*fParticlesPermuted)->DetEta(1, KLFitter::Particles::kParton))->p(blep_fit_e, blep_meas_e, &TFgoodTmp)));  // comp1
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log((*fDetector)->ResEnergyLightJet((*fParticlesPermuted)->DetEta(2, KLFitter::Particles::kParton))->p(lq1_fit_e, lq1_meas_e, &TFgoodTmp)));  // comp2
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log((*fDetector)->ResEnergyLightJet((*fParticlesPermuted)->DetEta(3, KLFitter::Particles::kParton))->p(lq2_fit_e, lq2_meas_e, &TFgoodTmp)));  // comp3
  if (!TFgoodTmp) fTFgood = false;

  // lepton energy resolution terms
  if (fTypeLepton == kElectron) {
    vecci.push_back(log(fResLepton->p(lep_fit_e, lep_meas_e, &TFgoodTmp)));  // comp4
  } else if (fTypeLepton == kMuon) {
    vecci.push_back(log(fResLepton->p(lep_fit_e* lep_meas_sintheta, lep_meas_pt, &TFgoodTmp)));  // comp4
  }
  if (!TFgoodTmp) fTFgood = false;

  // neutrino px and py
  vecci.push_back(log(fResMET->p(nu_fit_px, ETmiss_x, &TFgoodTmp, SumET)));  // comp5
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log(fResMET->p(nu_fit_py, ETmiss_y, &TFgoodTmp, SumET)));  // comp6
  if (!TFgoodTmp) fTFgood = false;

  // jet eta resolution terms
  vecci.push_back(log((*fDetector)->ResEtaBJet((*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kParton))->p(parameters[parBhadEta], (*fParticlesPermuted)->Parton(0)->Eta(), &TFgoodTmp)));
  if (!TFgoodTmp) fTFgood = false;
  vecci.push_back(log((*fDetector)->ResEtaBJet((*fParticlesPermuted)->DetEta(1, KLFitter::Particles::kParton))->p(parameters[parBlepEta], (*fParticlesPermuted)->Parton(1)->Eta(), &TFgoodTmp)));
  if (!TFgoodTmp) fTFgood = false;
  vecci.push_back(log((*fDetector)->ResEtaLightJet((*fParticlesPermuted)->DetEta(2, KLFitter::Particles::kParton))->p(parameters[parLQ1Eta], (*fParticlesPermuted)->Parton(2)->Eta(), &TFgoodTmp)));
  if (!TFgoodTmp) fTFgood = false;
  vecci.push_back(log((*fDetector)->ResEtaLightJet((*fParticlesPermuted)->DetEta(3, KLFitter::Particles::kParton))->p(parameters[parLQ2Eta], (*fParticlesPermuted)->Parton(3)->Eta(), &TFgoodTmp)));
  if (!TFgoodTmp) fTFgood = false;

  // jet phi resolution terms
  vecci.push_back(log((*fDetector)->ResPhiBJet((*fParticlesPermuted)->DetEta(0, KLFitter::Particles::kParton))->p(diffPhi(parameters[parBhadPhi], (*fParticlesPermuted)->Parton(0)->Phi()), 0., &TFgoodTmp)));
  if (!TFgoodTmp) fTFgood = false;
  vecci.push_back(log((*fDetector)->ResPhiBJet((*fParticlesPermuted)->DetEta(1, KLFitter::Particles::kParton))->p(diffPhi(parameters[parBlepPhi], (*fParticlesPermuted)->Parton(1)->Phi()), 0., &TFgoodTmp)));
  if (!TFgoodTmp) fTFgood = false;
  vecci.push_back(log((*fDetector)->ResPhiLightJet((*fParticlesPermuted)->DetEta(2, KLFitter::Particles::kParton))->p(diffPhi(parameters[parLQ1Phi], (*fParticlesPermuted)->Parton(2)->Phi()), 0., &TFgoodTmp)));
  if (!TFgoodTmp) fTFgood = false;
  vecci.push_back(log((*fDetector)->ResPhiLightJet((*fParticlesPermuted)->DetEta(3, KLFitter::Particles::kParton))->p(diffPhi(parameters[parLQ2Phi], (*fParticlesPermuted)->Parton(3)->Phi()), 0., &TFgoodTmp)));
  if (!TFgoodTmp) fTFgood = false;

  // physics constants
  double massW = fPhysicsConstants.MassW();
  double gammaW = fPhysicsConstants.GammaW();
  // note: top mass width should be made DEPENDENT on the top mass at a certain point
  //    fPhysicsConstants.SetMassTop(parameters[parTopM]);
  // (this will also set the correct width for the top)
  double gammaTop = fPhysicsConstants.GammaTop();

  // Breit-Wigner of hadronically decaying W-boson
  vecci.push_back(BCMath::LogBreitWignerRel(whad_fit_m, massW, gammaW));  // comp7

  // Breit-Wigner of leptonically decaying W-boson
  vecci.push_back(BCMath::LogBreitWignerRel(wlep_fit_m, massW, gammaW));  // comp8

  // Breit-Wigner of hadronically decaying top quark
  vecci.push_back(BCMath::LogBreitWignerRel(thad_fit_m, parameters[parTopM], gammaTop));  // comp9

  // Breit-Wigner of leptonically decaying top quark
  vecci.push_back(BCMath::LogBreitWignerRel(tlep_fit_m, parameters[parTopM], gammaTop));  // comp10

  // return log of likelihood
  return vecci;
}

// ---------------------------------------------------------
double KLFitter::LikelihoodTopLeptonJets_JetAngles::diffPhi(double phi1, double phi2) {
  double delta = phi1 - phi2;
  if (delta > TMath::Pi()) {
    delta -= TMath::TwoPi();
  } else if (delta < -TMath::Pi()) {
    delta += TMath::TwoPi();
  }
  return delta;
}
