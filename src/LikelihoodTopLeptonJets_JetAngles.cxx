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

namespace {
// A helper for calculating delta(phi1, phi2)
double diffPhi(double phi1, double phi2) {
  double delta = phi1 - phi2;
  if (delta > TMath::Pi()) {
    delta -= TMath::TwoPi();
  } else if (delta < -TMath::Pi()) {
    delta += TMath::TwoPi();
  }
  return delta;
}
}  // namespace


namespace KLFitter {
// ---------------------------------------------------------
LikelihoodTopLeptonJets_JetAngles::LikelihoodTopLeptonJets_JetAngles() = default;

// ---------------------------------------------------------
LikelihoodTopLeptonJets_JetAngles::~LikelihoodTopLeptonJets_JetAngles() = default;

// ---------------------------------------------------------
void LikelihoodTopLeptonJets_JetAngles::DefineParameters() {
  // add parameters of model
  AddParameter("energy hadronic b", fPhysicsConstants.MassBottom(), 1000.0);  // parBhadE
  AddParameter("energy leptonic b", fPhysicsConstants.MassBottom(), 1000.0);  // parBlepE
  AddParameter("energy light quark 1", 0.0, 1000.0);                          // parLQ1E
  AddParameter("energy light quark 2", 0.0, 1000.0);                          // parLQ2E
  AddParameter("energy lepton", 0.0, 1000.0);                                 // parLepE
  AddParameter("p_x neutrino", -1000.0, 1000.0);                              // parNuPx
  AddParameter("p_y neutrino", -1000.0, 1000.0);                              // parNuPy
  AddParameter("p_z neutrino", -1000.0, 1000.0);                              // parNuPz
  AddParameter("eta hadronic b", -2.5, 2.5);                                  // parBhadEta
  AddParameter("eta leptonic b", -2.5, 2.5);                                  // parBlepEta
  AddParameter("eta light quark 1", -2.5, 2.5);                               // parLQ1Eta
  AddParameter("eta light quark 2", -2.5, 2.5);                               // parLQ2Eta
  AddParameter("phi hadronic b", -TMath::Pi(), TMath::Pi());                  // parBhadPhi
  AddParameter("phi leptonic b", -TMath::Pi(), TMath::Pi());                  // parBlepPhi
  AddParameter("phi light quark 1", -TMath::Pi(), TMath::Pi());               // parLQ1Phi
  AddParameter("phi light quark 2", -TMath::Pi(), TMath::Pi());               // parLQ2Phi
  AddParameter("top mass", 100.0, 1000.0);                                    // parTopM
}

// ---------------------------------------------------------
int LikelihoodTopLeptonJets_JetAngles::CalculateLorentzVectors(std::vector <double> const& parameters) {
  double scale;
  double thad_fit_e;
  double thad_fit_px;
  double thad_fit_py;
  double thad_fit_pz;
  double tlep_fit_e;
  double tlep_fit_px;
  double tlep_fit_py;
  double tlep_fit_pz;

  TLorentzVector v;

  // hadronic b quark
  v.SetPtEtaPhiE(sqrt(parameters[parBhadE] * parameters[parBhadE] - m_bhad_meas_m * m_bhad_meas_m) / cosh(parameters[parBhadEta]),
                 parameters[parBhadEta],
                 parameters[parBhadPhi],
                 parameters[parBhadE]);
  m_bhad_fit_e = v.E();
  m_bhad_fit_px = v.Px();
  m_bhad_fit_py = v.Py();
  m_bhad_fit_pz = v.Pz();

  // leptonic b quark
  v.SetPtEtaPhiE(sqrt(parameters[parBlepE] * parameters[parBlepE] - m_blep_meas_m * m_blep_meas_m) / cosh(parameters[parBlepEta]),
                 parameters[parBlepEta],
                 parameters[parBlepPhi],
                 parameters[parBlepE]);
  m_blep_fit_e = v.E();
  m_blep_fit_px = v.Px();
  m_blep_fit_py = v.Py();
  m_blep_fit_pz = v.Pz();

  // light quark 1
  v.SetPtEtaPhiE(sqrt(parameters[parLQ1E] * parameters[parLQ1E] - m_lq1_meas_m * m_lq1_meas_m) / cosh(parameters[parLQ1Eta]),
                 parameters[parLQ1Eta],
                 parameters[parLQ1Phi],
                 parameters[parLQ1E]);
  m_lq1_fit_e = v.E();
  m_lq1_fit_px = v.Px();
  m_lq1_fit_py = v.Py();
  m_lq1_fit_pz = v.Pz();

  // light quark 2
  v.SetPtEtaPhiE(sqrt(parameters[parLQ2E] * parameters[parLQ2E] - m_lq2_meas_m * m_lq2_meas_m) / cosh(parameters[parLQ2Eta]),
                 parameters[parLQ2Eta],
                 parameters[parLQ2Phi],
                 parameters[parLQ2E]);
  m_lq2_fit_e = v.E();
  m_lq2_fit_px = v.Px();
  m_lq2_fit_py = v.Py();
  m_lq2_fit_pz = v.Pz();

  // lepton
  m_lep_fit_e = parameters[parLepE];
  scale = m_lep_fit_e / m_lep_meas_e;
  m_lep_fit_px = scale * m_lep_meas_px;
  m_lep_fit_py = scale * m_lep_meas_py;
  m_lep_fit_pz = scale * m_lep_meas_pz;

  // neutrino
  m_nu_fit_px = parameters[parNuPx];
  m_nu_fit_py = parameters[parNuPy];
  m_nu_fit_pz = parameters[parNuPz];
  m_nu_fit_e  = sqrt(m_nu_fit_px * m_nu_fit_px + m_nu_fit_py * m_nu_fit_py + m_nu_fit_pz * m_nu_fit_pz);

  // hadronic W
  m_whad_fit_e  = m_lq1_fit_e + m_lq2_fit_e;
  m_whad_fit_px = m_lq1_fit_px + m_lq2_fit_px;
  m_whad_fit_py = m_lq1_fit_py + m_lq2_fit_py;
  m_whad_fit_pz = m_lq1_fit_pz + m_lq2_fit_pz;
  m_whad_fit_m = sqrt(m_whad_fit_e * m_whad_fit_e - (m_whad_fit_px * m_whad_fit_px + m_whad_fit_py * m_whad_fit_py + m_whad_fit_pz * m_whad_fit_pz));

  // leptonic W
  m_wlep_fit_e  = m_lep_fit_e + m_nu_fit_e;
  m_wlep_fit_px = m_lep_fit_px + m_nu_fit_px;
  m_wlep_fit_py = m_lep_fit_py + m_nu_fit_py;
  m_wlep_fit_pz = m_lep_fit_pz + m_nu_fit_pz;
  m_wlep_fit_m = sqrt(m_wlep_fit_e * m_wlep_fit_e - (m_wlep_fit_px * m_wlep_fit_px + m_wlep_fit_py * m_wlep_fit_py + m_wlep_fit_pz * m_wlep_fit_pz));

  // hadronic top
  thad_fit_e = m_whad_fit_e + m_bhad_fit_e;
  thad_fit_px = m_whad_fit_px + m_bhad_fit_px;
  thad_fit_py = m_whad_fit_py + m_bhad_fit_py;
  thad_fit_pz = m_whad_fit_pz + m_bhad_fit_pz;
  m_thad_fit_m = sqrt(thad_fit_e * thad_fit_e - (thad_fit_px * thad_fit_px + thad_fit_py * thad_fit_py + thad_fit_pz * thad_fit_pz));

  // leptonic top
  tlep_fit_e = m_wlep_fit_e + m_blep_fit_e;
  tlep_fit_px = m_wlep_fit_px + m_blep_fit_px;
  tlep_fit_py = m_wlep_fit_py + m_blep_fit_py;
  tlep_fit_pz = m_wlep_fit_pz + m_blep_fit_pz;
  m_tlep_fit_m = sqrt(tlep_fit_e * tlep_fit_e - (tlep_fit_px * tlep_fit_px + tlep_fit_py * tlep_fit_py + tlep_fit_pz * tlep_fit_pz));

  // no error
  return 1;
}

// ---------------------------------------------------------
int LikelihoodTopLeptonJets_JetAngles::AdjustParameterRanges() {
  // adjust limits
  double nsigmas_jet    = m_flag_get_par_sigmas_from_TFs ? 10 : 7;
  double nsigmas_lepton = m_flag_get_par_sigmas_from_TFs ? 10 : 2;
  double nsigmas_met    = m_flag_get_par_sigmas_from_TFs ? 10 : 1;

  double E = (*fParticlesPermuted)->Parton(0)->E();
  double m = fPhysicsConstants.MassBottom();
  if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->Parton(0)->M());
  double sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_bhad->GetSigma(E) : sqrt(E);
  double Emin = std::max(m, E - nsigmas_jet* sigma);
  double Emax = E + nsigmas_jet* sigma;
  SetParameterRange(parBhadE, Emin, Emax);

  E = (*fParticlesPermuted)->Parton(1)->E();
  m = fPhysicsConstants.MassBottom();
  if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->Parton(1)->M());
  sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_blep->GetSigma(E) : sqrt(E);
  Emin = std::max(m, E - nsigmas_jet* sigma);
  Emax = E + nsigmas_jet* sigma;
  SetParameterRange(parBlepE, Emin, Emax);

  E = (*fParticlesPermuted)->Parton(2)->E();
  m = 0.001;
  if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->Parton(2)->M());
  sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_lq1->GetSigma(E) : sqrt(E);
  Emin = std::max(m, E - nsigmas_jet* sigma);
  Emax = E + nsigmas_jet* sigma;
  SetParameterRange(parLQ1E, Emin, Emax);

  E = (*fParticlesPermuted)->Parton(3)->E();
  m = 0.001;
  if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->Parton(3)->M());
  sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_lq2->GetSigma(E) : sqrt(E);
  Emin = std::max(m, E - nsigmas_jet* sigma);
  Emax = E + nsigmas_jet* sigma;
  SetParameterRange(parLQ2E, Emin, Emax);

  if (m_lepton_type == kElectron) {
    E = (*fParticlesPermuted)->Electron(0)->E();
    sigma = m_flag_get_par_sigmas_from_TFs ? m_res_lepton->GetSigma(E) : sqrt(E);
    Emin = std::max(0.001, E - nsigmas_lepton* sigma);
    Emax = E + nsigmas_lepton* sigma;
  } else if (m_lepton_type == kMuon) {
    E = (*fParticlesPermuted)->Muon(0)->E();
    double sintheta = sin((*fParticlesPermuted)->Muon(0)->Theta());
    sigma = m_flag_get_par_sigmas_from_TFs ? m_res_lepton->GetSigma(E*sintheta)/sintheta : E*E*sintheta;
    double sigrange = nsigmas_lepton* sigma;
    Emin = std::max(0.001, E -sigrange);
    Emax = E +sigrange;
  }
  SetParameterRange(parLepE, Emin, Emax);

  sigma = m_flag_get_par_sigmas_from_TFs ? m_res_met->GetSigma(m_et_miss_sum) : 100;
  double sigrange = nsigmas_met * sigma;
  SetParameterRange(parNuPx, m_et_miss_x - sigrange, m_et_miss_x + sigrange);
  SetParameterRange(parNuPy, m_et_miss_y - sigrange, m_et_miss_y + sigrange);

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

  if (m_flag_top_mass_fixed)
    SetParameterRange(parTopM, fPhysicsConstants.MassTop(), fPhysicsConstants.MassTop());

  // no error
  return 1;
}

// ---------------------------------------------------------
void LikelihoodTopLeptonJets_JetAngles::RequestResolutionFunctions() {
  (*fDetector)->RequestResolutionType(ResolutionType::EnergyLightJet);
  (*fDetector)->RequestResolutionType(ResolutionType::EnergyBJet);
  (*fDetector)->RequestResolutionType(ResolutionType::EnergyElectron);
  (*fDetector)->RequestResolutionType(ResolutionType::EnergyMuon);
  (*fDetector)->RequestResolutionType(ResolutionType::MissingET);
  (*fDetector)->RequestResolutionType(ResolutionType::EtaLightJet);
  (*fDetector)->RequestResolutionType(ResolutionType::EtaBJet);
  (*fDetector)->RequestResolutionType(ResolutionType::PhiLightJet);
  (*fDetector)->RequestResolutionType(ResolutionType::PhiBJet);
}

// ---------------------------------------------------------
double LikelihoodTopLeptonJets_JetAngles::LogLikelihood(const std::vector<double> & parameters) {
  // calculate 4-vectors
  CalculateLorentzVectors(parameters);

  // define log of likelihood
  double logprob(0.);

  // temporary flag for a safe use of the transfer functions
  bool TFgoodTmp(true);

  // jet energy resolution terms
  logprob += log(m_res_energy_bhad->p(m_bhad_fit_e, m_bhad_meas_e, &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;

  logprob += log(m_res_energy_blep->p(m_blep_fit_e, m_blep_meas_e, &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;

  logprob += log(m_res_energy_lq1->p(m_lq1_fit_e, m_lq1_meas_e, &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;

  logprob += log(m_res_energy_lq2->p(m_lq2_fit_e, m_lq2_meas_e, &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;

  // lepton energy resolution terms
  if (m_lepton_type == kElectron) {
    logprob += log(m_res_lepton->p(m_lep_fit_e, m_lep_meas_e, &TFgoodTmp));
  } else if (m_lepton_type == kMuon) {
    logprob += log(m_res_lepton->p(m_lep_fit_e* m_lep_meas_sintheta, m_lep_meas_pt, &TFgoodTmp));
  }
  if (!TFgoodTmp) fTFgood = false;

  // neutrino px and py
  logprob += log(m_res_met->p(m_nu_fit_px, m_et_miss_x, &TFgoodTmp, m_et_miss_sum));
  if (!TFgoodTmp) fTFgood = false;

  logprob += log(m_res_met->p(m_nu_fit_py, m_et_miss_y, &TFgoodTmp, m_et_miss_sum));
  if (!TFgoodTmp) fTFgood = false;

  // eta resolution
  logprob += log((*fDetector)->ResEtaBJet((*fParticlesPermuted)->DetEta(0, Particles::kParton))->p(parameters[parBhadEta], (*fParticlesPermuted)->Parton(0)->Eta(), &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;
  logprob += log((*fDetector)->ResEtaBJet((*fParticlesPermuted)->DetEta(1, Particles::kParton))->p(parameters[parBlepEta], (*fParticlesPermuted)->Parton(1)->Eta(), &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;
  logprob += log((*fDetector)->ResEtaLightJet((*fParticlesPermuted)->DetEta(2, Particles::kParton))->p(parameters[parLQ1Eta], (*fParticlesPermuted)->Parton(2)->Eta(), &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;
  logprob += log((*fDetector)->ResEtaLightJet((*fParticlesPermuted)->DetEta(3, Particles::kParton))->p(parameters[parLQ2Eta], (*fParticlesPermuted)->Parton(3)->Eta(), &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;

  // transform all phi values, so that they are centered around zero, and not around the measured phi

  // phi resolution
  logprob += log((*fDetector)->ResPhiBJet((*fParticlesPermuted)->DetEta(0, Particles::kParton))->p(diffPhi(parameters[parBhadPhi], (*fParticlesPermuted)->Parton(0)->Phi()), 0., &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;
  logprob += log((*fDetector)->ResPhiBJet((*fParticlesPermuted)->DetEta(1, Particles::kParton))->p(diffPhi(parameters[parBlepPhi], (*fParticlesPermuted)->Parton(1)->Phi()), 0., &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;
  logprob += log((*fDetector)->ResPhiLightJet((*fParticlesPermuted)->DetEta(2, Particles::kParton))->p(diffPhi(parameters[parLQ1Phi], (*fParticlesPermuted)->Parton(2)->Phi()), 0., &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;
  logprob += log((*fDetector)->ResPhiLightJet((*fParticlesPermuted)->DetEta(3, Particles::kParton))->p(diffPhi(parameters[parLQ2Phi], (*fParticlesPermuted)->Parton(3)->Phi()), 0., &TFgoodTmp));
  if (!TFgoodTmp) fTFgood = false;

  // physics constants
  double massW = fPhysicsConstants.MassW();
  double gammaW = fPhysicsConstants.GammaW();
  double gammaTop = fPhysicsConstants.GammaTop();

  // Breit-Wigner of hadronically decaying W-boson
  logprob += BCMath::LogBreitWignerRel(m_whad_fit_m, massW, gammaW);

  // Breit-Wigner of leptonically decaying W-boson
  logprob += BCMath::LogBreitWignerRel(m_wlep_fit_m, massW, gammaW);

  // Breit-Wigner of hadronically decaying top quark
  logprob += BCMath::LogBreitWignerRel(m_thad_fit_m, parameters[parTopM], gammaTop);

  // Breit-Wigner of leptonically decaying top quark
  logprob += BCMath::LogBreitWignerRel(m_tlep_fit_m, parameters[parTopM], gammaTop);

  // return log of likelihood
  return logprob;
}

// ---------------------------------------------------------
std::vector<double> LikelihoodTopLeptonJets_JetAngles::GetInitialParametersWoNeutrinoPz() {
  std::vector<double> values(GetNParameters());

  // energies of the quarks
  values[parBhadE] = m_bhad_meas_e;
  values[parBlepE] = m_blep_meas_e;
  values[parLQ1E]  = m_lq1_meas_e;
  values[parLQ2E]  = m_lq2_meas_e;

  // energy of the lepton
  if (m_lepton_type == kElectron) {
    values[parLepE] = (*fParticlesPermuted)->Electron(0)->E();
  } else if (m_lepton_type == kMuon) {
    values[parLepE] = (*fParticlesPermuted)->Muon(0)->E();
  }

  // missing px and py
  values[parNuPx] = m_et_miss_x;
  values[parNuPy] = m_et_miss_y;

  // eta of the quarks
  values[parBhadEta] = m_bhad_meas_eta;
  values[parBlepEta] = m_blep_meas_eta;
  values[parLQ1Eta]  = m_lq1_meas_eta;
  values[parLQ2Eta]  = m_lq2_meas_eta;

  // phi of the quarks
  values[parBhadPhi] = m_bhad_meas_phi;
  values[parBlepPhi] = m_blep_meas_phi;
  values[parLQ1Phi]  = m_lq1_meas_phi;
  values[parLQ2Phi]  = m_lq2_meas_phi;

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
std::vector<double> LikelihoodTopLeptonJets_JetAngles::LogLikelihoodComponents(std::vector<double> parameters) {
  std::vector<double> vecci;

  // calculate 4-vectors
  CalculateLorentzVectors(parameters);

  // temporary flag for a safe use of the transfer functions
  bool TFgoodTmp(true);

  // jet energy resolution terms
  vecci.push_back(log(m_res_energy_bhad->p(m_bhad_fit_e, m_bhad_meas_e, &TFgoodTmp)));  // comp0
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log(m_res_energy_blep->p(m_blep_fit_e, m_blep_meas_e, &TFgoodTmp)));  // comp1
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log(m_res_energy_lq1->p(m_lq1_fit_e, m_lq1_meas_e, &TFgoodTmp)));  // comp2
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log(m_res_energy_lq2->p(m_lq2_fit_e, m_lq2_meas_e, &TFgoodTmp)));  // comp3
  if (!TFgoodTmp) fTFgood = false;

  // lepton energy resolution terms
  if (m_lepton_type == kElectron) {
    vecci.push_back(log(m_res_lepton->p(m_lep_fit_e, m_lep_meas_e, &TFgoodTmp)));  // comp4
  } else if (m_lepton_type == kMuon) {
    vecci.push_back(log(m_res_lepton->p(m_lep_fit_e* m_lep_meas_sintheta, m_lep_meas_pt, &TFgoodTmp)));  // comp4
  }
  if (!TFgoodTmp) fTFgood = false;

  // neutrino px and py
  vecci.push_back(log(m_res_met->p(m_nu_fit_px, m_et_miss_x, &TFgoodTmp, m_et_miss_sum)));  // comp5
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(log(m_res_met->p(m_nu_fit_py, m_et_miss_y, &TFgoodTmp, m_et_miss_sum)));  // comp6
  if (!TFgoodTmp) fTFgood = false;

  // jet eta resolution terms
  vecci.push_back(log((*fDetector)->ResEtaBJet((*fParticlesPermuted)->DetEta(0, Particles::kParton))->p(parameters[parBhadEta], (*fParticlesPermuted)->Parton(0)->Eta(), &TFgoodTmp)));
  if (!TFgoodTmp) fTFgood = false;
  vecci.push_back(log((*fDetector)->ResEtaBJet((*fParticlesPermuted)->DetEta(1, Particles::kParton))->p(parameters[parBlepEta], (*fParticlesPermuted)->Parton(1)->Eta(), &TFgoodTmp)));
  if (!TFgoodTmp) fTFgood = false;
  vecci.push_back(log((*fDetector)->ResEtaLightJet((*fParticlesPermuted)->DetEta(2, Particles::kParton))->p(parameters[parLQ1Eta], (*fParticlesPermuted)->Parton(2)->Eta(), &TFgoodTmp)));
  if (!TFgoodTmp) fTFgood = false;
  vecci.push_back(log((*fDetector)->ResEtaLightJet((*fParticlesPermuted)->DetEta(3, Particles::kParton))->p(parameters[parLQ2Eta], (*fParticlesPermuted)->Parton(3)->Eta(), &TFgoodTmp)));
  if (!TFgoodTmp) fTFgood = false;

  // jet phi resolution terms
  vecci.push_back(log((*fDetector)->ResPhiBJet((*fParticlesPermuted)->DetEta(0, Particles::kParton))->p(diffPhi(parameters[parBhadPhi], (*fParticlesPermuted)->Parton(0)->Phi()), 0., &TFgoodTmp)));
  if (!TFgoodTmp) fTFgood = false;
  vecci.push_back(log((*fDetector)->ResPhiBJet((*fParticlesPermuted)->DetEta(1, Particles::kParton))->p(diffPhi(parameters[parBlepPhi], (*fParticlesPermuted)->Parton(1)->Phi()), 0., &TFgoodTmp)));
  if (!TFgoodTmp) fTFgood = false;
  vecci.push_back(log((*fDetector)->ResPhiLightJet((*fParticlesPermuted)->DetEta(2, Particles::kParton))->p(diffPhi(parameters[parLQ1Phi], (*fParticlesPermuted)->Parton(2)->Phi()), 0., &TFgoodTmp)));
  if (!TFgoodTmp) fTFgood = false;
  vecci.push_back(log((*fDetector)->ResPhiLightJet((*fParticlesPermuted)->DetEta(3, Particles::kParton))->p(diffPhi(parameters[parLQ2Phi], (*fParticlesPermuted)->Parton(3)->Phi()), 0., &TFgoodTmp)));
  if (!TFgoodTmp) fTFgood = false;

  // physics constants
  double massW = fPhysicsConstants.MassW();
  double gammaW = fPhysicsConstants.GammaW();
  double gammaTop = fPhysicsConstants.GammaTop();

  // Breit-Wigner of hadronically decaying W-boson
  vecci.push_back(BCMath::LogBreitWignerRel(m_whad_fit_m, massW, gammaW));  // comp7

  // Breit-Wigner of leptonically decaying W-boson
  vecci.push_back(BCMath::LogBreitWignerRel(m_wlep_fit_m, massW, gammaW));  // comp8

  // Breit-Wigner of hadronically decaying top quark
  vecci.push_back(BCMath::LogBreitWignerRel(m_thad_fit_m, parameters[parTopM], gammaTop));  // comp9

  // Breit-Wigner of leptonically decaying top quark
  vecci.push_back(BCMath::LogBreitWignerRel(m_tlep_fit_m, parameters[parTopM], gammaTop));  // comp10

  // return log of likelihood
  return vecci;
}
}  // namespace KLFitter
