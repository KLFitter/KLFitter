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

#include "KLFitter/LikelihoodTopLeptonJets.h"

#include <algorithm>
#include <iostream>

#include "BAT/BCMath.h"
#include "BAT/BCParameter.h"
#include "KLFitter/DetectorBase.h"
#include "KLFitter/ParticleCollection.h"
#include "KLFitter/PermutationHandler.h"
#include "KLFitter/PhysicsConstants.h"
#include "KLFitter/ResolutionBase.h"
#include "TLorentzVector.h"

namespace KLFitter {
// ---------------------------------------------------------
LikelihoodTopLeptonJets::LikelihoodTopLeptonJets()
  : LikelihoodBase::LikelihoodBase()
  , m_flag_top_mass_fixed(false)
  , m_flag_get_par_sigmas_from_TFs(false)
  , m_et_miss_x(0.)
  , m_et_miss_y(0.)
  , m_et_miss_sum(0.)
  , m_lepton_type(kElectron) {
  // define model particles
  this->DefineModelParticles();

  // define parameters
  this->DefineParameters();
}

// ---------------------------------------------------------
LikelihoodTopLeptonJets::~LikelihoodTopLeptonJets() = default;

// ---------------------------------------------------------
int LikelihoodTopLeptonJets::SetET_miss_XY_SumET(double etx, double ety, double sumet) {
  // set missing ET x and y component and the m_et_miss_sum
  m_et_miss_x = etx;
  m_et_miss_y = ety;
  m_et_miss_sum = sumet;

  // no error
  return 1;
}

// ---------------------------------------------------------
void LikelihoodTopLeptonJets::RequestResolutionFunctions() {
  (*fDetector)->RequestResolutionType(ResolutionType::EnergyLightJet);
  (*fDetector)->RequestResolutionType(ResolutionType::EnergyBJet);
  (*fDetector)->RequestResolutionType(ResolutionType::EnergyElectron);
  (*fDetector)->RequestResolutionType(ResolutionType::EnergyMuon);
  (*fDetector)->RequestResolutionType(ResolutionType::MissingET);
}

// ---------------------------------------------------------
void LikelihoodTopLeptonJets::SetLeptonType(LeptonType leptontype) {
  if (leptontype != kElectron && leptontype != kMuon) {
    std::cout << "KLFitter::SetLeptonType(). Warning: lepton type not defined. Set electron as lepton type." << std::endl;
    m_lepton_type = kElectron;
  } else {
    m_lepton_type = leptontype;
  }

  // define model particles
  DefineModelParticles();
}

// ---------------------------------------------------------
void LikelihoodTopLeptonJets::SetLeptonType(int leptontype) {
  if (leptontype != 1 && leptontype != 2) {
    std::cout << "KLFitter::SetLeptonType(). Warning: lepton type not defined. Set electron as lepton type." << std::endl;
    leptontype = 1;
  }

  if (leptontype == 1) {
    SetLeptonType(kElectron);
  } else if (leptontype == 2) {
    SetLeptonType(kMuon);
  }
}

// ---------------------------------------------------------
int LikelihoodTopLeptonJets::DefineModelParticles() {
  // create the particles of the model
  fParticlesModel.reset(new ParticleCollection{});

  // add model particles
  Particles::Parton parton0{"hadronic b quark", TLorentzVector{}};
  parton0.SetIdentifier(0);
  parton0.SetTrueFlavor(Particles::PartonTrueFlavor::kB);
  fParticlesModel->AddParticle(parton0);

  Particles::Parton parton1{"leptonic b quark", TLorentzVector{}};
  parton1.SetIdentifier(1);
  parton1.SetTrueFlavor(Particles::PartonTrueFlavor::kB);
  fParticlesModel->AddParticle(parton1);

  Particles::Parton parton2{"light quark 1", TLorentzVector{}};
  parton2.SetIdentifier(2);
  parton2.SetTrueFlavor(Particles::PartonTrueFlavor::kLight);
  fParticlesModel->AddParticle(parton2);

  Particles::Parton parton3{"light quark 2", TLorentzVector{}};
  parton3.SetIdentifier(3);
  parton3.SetTrueFlavor(Particles::PartonTrueFlavor::kLight);
  fParticlesModel->AddParticle(parton3);

  if (m_lepton_type == kElectron) {
    fParticlesModel->AddParticle(Particles::Electron{"electron", TLorentzVector{}});
  } else if (m_lepton_type == kMuon) {
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
void LikelihoodTopLeptonJets::DefineParameters() {
  // add parameters of model
  AddParameter("energy hadronic b", fPhysicsConstants.MassBottom(), 1000.0);  // parBhadE
  AddParameter("energy leptonic b", fPhysicsConstants.MassBottom(), 1000.0);  // parBlepE
  AddParameter("energy light quark 1", 0.0, 1000.0);                          // parLQ1E
  AddParameter("energy light quark 2", 0.0, 1000.0);                          // parLQ2E
  AddParameter("energy lepton", 0.0, 1000.0);                                 // parLepE
  AddParameter("p_x neutrino", -1000.0, 1000.0);                              // parNuPx
  AddParameter("p_y neutrino", -1000.0, 1000.0);                              // parNuPy
  AddParameter("p_z neutrino", -1000.0, 1000.0);                              // parNuPz
  AddParameter("top mass", 100.0, 1000.0);                                   // parTopM
}

// ---------------------------------------------------------
int LikelihoodTopLeptonJets::CalculateLorentzVectors(std::vector <double> const& parameters) {
  double scale;
  double thad_fit_e;
  double thad_fit_px;
  double thad_fit_py;
  double thad_fit_pz;
  double tlep_fit_e;
  double tlep_fit_px;
  double tlep_fit_py;
  double tlep_fit_pz;

  // hadronic b quark
  m_bhad_fit_e = parameters[parBhadE];
  scale = sqrt(m_bhad_fit_e * m_bhad_fit_e - m_bhad_meas_m * m_bhad_meas_m) / m_bhad_meas_p;
  m_bhad_fit_px = scale * m_bhad_meas_px;
  m_bhad_fit_py = scale * m_bhad_meas_py;
  m_bhad_fit_pz = scale * m_bhad_meas_pz;

  // leptonic b quark
  m_blep_fit_e = parameters[parBlepE];
  scale = sqrt(m_blep_fit_e * m_blep_fit_e - m_blep_meas_m * m_blep_meas_m) / m_blep_meas_p;
  m_blep_fit_px = scale * m_blep_meas_px;
  m_blep_fit_py = scale * m_blep_meas_py;
  m_blep_fit_pz = scale * m_blep_meas_pz;

  // light quark 1
  m_lq1_fit_e = parameters[parLQ1E];
  scale = sqrt(m_lq1_fit_e * m_lq1_fit_e - m_lq1_meas_m * m_lq1_meas_m) / m_lq1_meas_p;
  m_lq1_fit_px = scale * m_lq1_meas_px;
  m_lq1_fit_py = scale * m_lq1_meas_py;
  m_lq1_fit_pz = scale * m_lq1_meas_pz;

  // light quark 2
  m_lq2_fit_e = parameters[parLQ2E];
  scale = sqrt(m_lq2_fit_e * m_lq2_fit_e - m_lq2_meas_m * m_lq2_meas_m) / m_lq2_meas_p;
  m_lq2_fit_px  = scale * m_lq2_meas_px;
  m_lq2_fit_py  = scale * m_lq2_meas_py;
  m_lq2_fit_pz  = scale * m_lq2_meas_pz;

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
int LikelihoodTopLeptonJets::RemoveInvariantParticlePermutations() {
  // error code
  int err = 1;

  // remove the permutation from the second and the third jet
  Particles::Type ptype = Particles::Type::kParton;
  err *= (*fPermutations)->InvariantParticlePermutations(ptype, std::vector<int>{2, 3});

  // remove invariant jet permutations of all jets not considered
  const ParticleCollection* particles = (*fPermutations)->Particles();
  std::vector<int> indices;
  for (size_t iPartons = 4; iPartons < particles->partons.size(); iPartons++) {
    indices.emplace_back(iPartons);
  }
  err *= (*fPermutations)->InvariantParticlePermutations(ptype, indices);

  // remove the permutation from the other lepton
  if (m_lepton_type == kElectron) {
    ptype = Particles::Type::kMuon;
    indices.clear();
    for (size_t iMuon = 0; iMuon < particles->muons.size(); iMuon++) {
      indices.emplace_back(iMuon);
    }
    err *= (*fPermutations)->InvariantParticlePermutations(ptype, indices);
  } else if (m_lepton_type == kMuon) {
    ptype = Particles::Type::kElectron;
    indices.clear();
    for (size_t iElectron = 0; iElectron < particles->electrons.size(); iElectron++) {
      indices.emplace_back(iElectron);
    }
    err *= (*fPermutations)->InvariantParticlePermutations(ptype, indices);
  }

  // return error code
  return err;
}

// ---------------------------------------------------------
int LikelihoodTopLeptonJets::AdjustParameterRanges() {
  // adjust limits
  double nsigmas_jet    = m_flag_get_par_sigmas_from_TFs ? 10 : 7;
  double nsigmas_lepton = m_flag_get_par_sigmas_from_TFs ? 10 : 2;
  double nsigmas_met    = m_flag_get_par_sigmas_from_TFs ? 10 : 1;

  double E = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 0)->E();
  double m = fPhysicsConstants.MassBottom();
  if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 0)->M());
  double sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_bhad->GetSigma(E) : sqrt(E);
  double Emin = std::max(m, E - nsigmas_jet * sigma);
  double Emax = E + nsigmas_jet * sigma;
  SetParameterRange(parBhadE, Emin, Emax);

  E = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 1)->E();
  m = fPhysicsConstants.MassBottom();
  if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 1)->M());
  sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_blep->GetSigma(E) : sqrt(E);
  Emin = std::max(m, E - nsigmas_jet * sigma);
  Emax = E + nsigmas_jet * sigma;
  SetParameterRange(parBlepE, Emin, Emax);

  E = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 2)->E();
  m = 0.001;
  if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 2)->M());
  sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_lq1->GetSigma(E) : sqrt(E);
  Emin = std::max(m, E - nsigmas_jet * sigma);
  Emax = E + nsigmas_jet * sigma;
  SetParameterRange(parLQ1E, Emin, Emax);

  E = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 3)->E();
  m = 0.001;
  if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 3)->M());
  sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_lq2->GetSigma(E) : sqrt(E);
  Emin = std::max(m, E - nsigmas_jet * sigma);
  Emax = E + nsigmas_jet * sigma;
  SetParameterRange(parLQ2E, Emin, Emax);

  if (m_lepton_type == kElectron) {
    E = (*fParticlesPermuted)->GetP4(Particles::Type::kElectron, 0)->E();
    sigma = m_flag_get_par_sigmas_from_TFs ? m_res_lepton->GetSigma(E) : sqrt(E);
    Emin = std::max(0.001, E - nsigmas_lepton * sigma);
    Emax = E + nsigmas_lepton * sigma;
  } else if (m_lepton_type == kMuon) {
    E = (*fParticlesPermuted)->GetP4(Particles::Type::kMuon, 0)->E();
    double sintheta = sin((*fParticlesPermuted)->GetP4(Particles::Type::kMuon, 0)->Theta());
    sigma = m_flag_get_par_sigmas_from_TFs ? m_res_lepton->GetSigma(E * sintheta) / sintheta : E * E * sintheta;
    double sigrange = nsigmas_lepton * sigma;
    Emin = std::max(0.001, E - sigrange);
    Emax = E + sigrange;
  }
  SetParameterRange(parLepE, Emin, Emax);

  sigma = m_flag_get_par_sigmas_from_TFs ? m_res_met->GetSigma(m_et_miss_sum) : 100;
  double sigrange = nsigmas_met * sigma;
  SetParameterRange(parNuPx, m_et_miss_x - sigrange, m_et_miss_x + sigrange);
  SetParameterRange(parNuPy, m_et_miss_y - sigrange, m_et_miss_y + sigrange);

  if (m_flag_top_mass_fixed)
    SetParameterRange(parTopM, fPhysicsConstants.MassTop(), fPhysicsConstants.MassTop());

  // no error
  return 1;
}

// ---------------------------------------------------------
double LikelihoodTopLeptonJets::LogLikelihood(const std::vector<double> & parameters) {
  // calculate 4-vectors
  CalculateLorentzVectors(parameters);

  // define log of likelihood
  double logprob(0.);

  // temporary flag for a safe use of the transfer functions
  bool TFgoodTmp(true);

  // jet energy resolution terms
  logprob += m_res_energy_bhad->logp(m_bhad_fit_e, m_bhad_meas_e, &TFgoodTmp);
  if (!TFgoodTmp) fTFgood = false;

  logprob += m_res_energy_blep->logp(m_blep_fit_e, m_blep_meas_e, &TFgoodTmp);
  if (!TFgoodTmp) fTFgood = false;

  logprob += m_res_energy_lq1->logp(m_lq1_fit_e, m_lq1_meas_e, &TFgoodTmp);
  if (!TFgoodTmp) fTFgood = false;

  logprob += m_res_energy_lq2->logp(m_lq2_fit_e, m_lq2_meas_e, &TFgoodTmp);
  if (!TFgoodTmp) fTFgood = false;

  // lepton energy resolution terms
  if (m_lepton_type == kElectron) {
    logprob += m_res_lepton->logp(m_lep_fit_e, m_lep_meas_e, &TFgoodTmp);
  } else if (m_lepton_type == kMuon) {
    logprob += m_res_lepton->logp(m_lep_fit_e * m_lep_meas_sintheta, m_lep_meas_pt, &TFgoodTmp);
  }
  if (!TFgoodTmp) fTFgood = false;

  // neutrino px and py
  logprob += m_res_met->logp(m_nu_fit_px, m_et_miss_x, &TFgoodTmp, m_et_miss_sum);
  if (!TFgoodTmp) fTFgood = false;

  logprob += m_res_met->logp(m_nu_fit_py, m_et_miss_y, &TFgoodTmp, m_et_miss_sum);
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
std::vector<double> LikelihoodTopLeptonJets::GetInitialParameters() {
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

    if (sol1 > sol2) values[parNuPz] = neutrino_pz_solutions[0];
  }

  return values;
}

// ---------------------------------------------------------
std::vector<double> LikelihoodTopLeptonJets::GetInitialParametersWoNeutrinoPz() {
  std::vector<double> values(GetNParameters());

  // energies of the quarks
  values[parBhadE] = m_bhad_meas_e;
  values[parBlepE] = m_blep_meas_e;
  values[parLQ1E]  = m_lq1_meas_e;
  values[parLQ2E]  = m_lq2_meas_e;

  // energy of the lepton
  if (m_lepton_type == kElectron) {
    values[parLepE] = (*fParticlesPermuted)->GetP4(Particles::Type::kElectron, 0)->E();
  } else if (m_lepton_type == kMuon) {
    values[parLepE] = (*fParticlesPermuted)->GetP4(Particles::Type::kMuon, 0)->E();
  }

  // missing px and py
  values[parNuPx] = m_et_miss_x;
  values[parNuPy] = m_et_miss_y;

  // pz of the neutrino
  values[parNuPz] = 0.;

  // top mass
  double mtop = (*(*fParticlesPermuted)->GetP4(Particles::Type::kParton, 0) + *(*fParticlesPermuted)->GetP4(Particles::Type::kParton, 2) + *(*fParticlesPermuted)->GetP4(Particles::Type::kParton, 3)).M();
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
std::vector<double> LikelihoodTopLeptonJets::GetNeutrinoPzSolutions() {
  return CalculateNeutrinoPzSolutions();
}

// ---------------------------------------------------------
std::vector<double> LikelihoodTopLeptonJets::CalculateNeutrinoPzSolutions(TLorentzVector* additionalParticle) {
  std::vector<double> pz;

  class PhysicsConstants constants;
  // electron mass
  double mE = 0.;

  double px_c = 0.0;
  double py_c = 0.0;
  double pz_c = 0.0;
  double Ec = 0.0;

  if (m_lepton_type == kElectron) {
    px_c = (*fParticlesPermuted)->GetP4(Particles::Type::kElectron, 0)->Px();
    py_c = (*fParticlesPermuted)->GetP4(Particles::Type::kElectron, 0)->Py();
    pz_c = (*fParticlesPermuted)->GetP4(Particles::Type::kElectron, 0)->Pz();
    Ec = (*fParticlesPermuted)->GetP4(Particles::Type::kElectron, 0)->E();
  } else if (m_lepton_type == kMuon) {
    px_c = (*fParticlesPermuted)->GetP4(Particles::Type::kMuon, 0)->Px();
    py_c = (*fParticlesPermuted)->GetP4(Particles::Type::kMuon, 0)->Py();
    pz_c = (*fParticlesPermuted)->GetP4(Particles::Type::kMuon, 0)->Pz();
    Ec = (*fParticlesPermuted)->GetP4(Particles::Type::kMuon, 0)->E();
  }

  // add additional particle to "charged lepton" 4-vector
  if (additionalParticle) {
    px_c += additionalParticle->Px();
    py_c += additionalParticle->Py();
    pz_c += additionalParticle->Pz();
    Ec += additionalParticle->E();
  }

  double px_nu = m_et_miss_x;
  double py_nu = m_et_miss_y;
  double alpha = constants.MassW() * constants.MassW() - mE * mE + 2 * (px_c * px_nu + py_c * py_nu);

  double a = pz_c * pz_c - Ec * Ec;
  double b = alpha * pz_c;
  double c = - Ec * Ec * (px_nu * px_nu + py_nu * py_nu) + alpha * alpha / 4.;

  double discriminant = b * b - 4 * a * c;
  if (discriminant < 0.) return pz;

  double pz_offset = - b / (2 * a);

  double squareRoot = sqrt(discriminant);
  if (squareRoot < 1.e-6) {
    pz.push_back(pz_offset);
  } else {
    pz.push_back(pz_offset + squareRoot / (2 * a));
    pz.push_back(pz_offset - squareRoot / (2 * a));
  }

  return pz;
}

// ---------------------------------------------------------
int LikelihoodTopLeptonJets::SavePermutedParticles() {
  m_bhad_meas_e      = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 0)->E();
  m_bhad_meas_deteta = (*fParticlesPermuted)->partons.at(0).GetDetEta();
  m_bhad_meas_px     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 0)->Px();
  m_bhad_meas_py     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 0)->Py();
  m_bhad_meas_pz     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 0)->Pz();
  m_bhad_meas_m      = SetPartonMass((*fParticlesPermuted)->GetP4(Particles::Type::kParton, 0)->M(), fPhysicsConstants.MassBottom(), &m_bhad_meas_px, &m_bhad_meas_py, &m_bhad_meas_pz, m_bhad_meas_e);
  m_bhad_meas_p      = sqrt(m_bhad_meas_e*m_bhad_meas_e - m_bhad_meas_m*m_bhad_meas_m);

  m_blep_meas_e      = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 1)->E();
  m_blep_meas_deteta = (*fParticlesPermuted)->partons.at(1).GetDetEta();
  m_blep_meas_px     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 1)->Px();
  m_blep_meas_py     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 1)->Py();
  m_blep_meas_pz     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 1)->Pz();
  m_blep_meas_m      = SetPartonMass((*fParticlesPermuted)->GetP4(Particles::Type::kParton, 1)->M(), fPhysicsConstants.MassBottom(), &m_blep_meas_px, &m_blep_meas_py, &m_blep_meas_pz, m_blep_meas_e);
  m_blep_meas_p      = sqrt(m_blep_meas_e*m_blep_meas_e - m_blep_meas_m*m_blep_meas_m);

  m_lq1_meas_e      = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 2)->E();
  m_lq1_meas_deteta = (*fParticlesPermuted)->partons.at(2).GetDetEta();
  m_lq1_meas_px     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 2)->Px();
  m_lq1_meas_py     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 2)->Py();
  m_lq1_meas_pz     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 2)->Pz();
  m_lq1_meas_m      = SetPartonMass((*fParticlesPermuted)->GetP4(Particles::Type::kParton, 2)->M(), 0., &m_lq1_meas_px, &m_lq1_meas_py, &m_lq1_meas_pz, m_lq1_meas_e);
  m_lq1_meas_p      = sqrt(m_lq1_meas_e*m_lq1_meas_e - m_lq1_meas_m*m_lq1_meas_m);

  m_lq2_meas_e      = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 3)->E();
  m_lq2_meas_deteta = (*fParticlesPermuted)->partons.at(3).GetDetEta();
  m_lq2_meas_px     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 3)->Px();
  m_lq2_meas_py     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 3)->Py();
  m_lq2_meas_pz     = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 3)->Pz();
  m_lq2_meas_m      = SetPartonMass((*fParticlesPermuted)->GetP4(Particles::Type::kParton, 3)->M(), 0., &m_lq2_meas_px, &m_lq2_meas_py, &m_lq2_meas_pz, m_lq2_meas_e);
  m_lq2_meas_p      = sqrt(m_lq2_meas_e*m_lq2_meas_e - m_lq2_meas_m*m_lq2_meas_m);

  TLorentzVector* lepton(0);
  if (m_lepton_type == kElectron) {
    lepton = (*fParticlesPermuted)->GetP4(Particles::Type::kElectron, 0);
    m_lep_meas_deteta = (*fParticlesPermuted)->electrons.at(0).GetDetEta();
  } else {
    lepton = (*fParticlesPermuted)->GetP4(Particles::Type::kMuon, 0);
    m_lep_meas_deteta = (*fParticlesPermuted)->muons.at(0).GetDetEta();
  }
  m_lep_meas_e        = lepton->E();
  m_lep_meas_sintheta = sin(lepton->Theta());
  m_lep_meas_pt       = lepton->Pt();
  m_lep_meas_px       = lepton->Px();
  m_lep_meas_py       = lepton->Py();
  m_lep_meas_pz       = lepton->Pz();

  // no error
  return 1;
}

// ---------------------------------------------------------
int LikelihoodTopLeptonJets::SaveResolutionFunctions() {
  m_res_energy_bhad = (*fDetector)->ResEnergyBJet(m_bhad_meas_deteta);
  m_res_energy_blep = (*fDetector)->ResEnergyBJet(m_blep_meas_deteta);
  m_res_energy_lq1  = (*fDetector)->ResEnergyLightJet(m_lq1_meas_deteta);
  m_res_energy_lq2  = (*fDetector)->ResEnergyLightJet(m_lq2_meas_deteta);
  if (m_lepton_type == kElectron) {
    m_res_lepton = (*fDetector)->ResEnergyElectron(m_lep_meas_deteta);
  } else if (m_lepton_type == kMuon) {
    m_res_lepton = (*fDetector)->ResEnergyMuon(m_lep_meas_deteta);
  }
  m_res_met = (*fDetector)->ResMissingET();

  // no error
  return 1;
}

// ---------------------------------------------------------
int LikelihoodTopLeptonJets::BuildModelParticles() {
  if (GetBestFitParameters().size() > 0) CalculateLorentzVectors(GetBestFitParameters());

  TLorentzVector* bhad = fParticlesModel->GetP4(Particles::Type::kParton, 0);
  TLorentzVector* blep = fParticlesModel->GetP4(Particles::Type::kParton, 1);
  TLorentzVector* lq1  = fParticlesModel->GetP4(Particles::Type::kParton, 2);
  TLorentzVector* lq2  = fParticlesModel->GetP4(Particles::Type::kParton, 3);
  TLorentzVector* lep(0);
  if (m_lepton_type == kElectron) {
    lep  = fParticlesModel->GetP4(Particles::Type::kElectron, 0);
  } else if (m_lepton_type == kMuon) {
    lep  = fParticlesModel->GetP4(Particles::Type::kMuon, 0);
  }
  TLorentzVector* nu   = fParticlesModel->GetP4(Particles::Type::kNeutrino, 0);
  TLorentzVector* whad  = fParticlesModel->GetP4(Particles::Type::kBoson, 0);
  TLorentzVector* wlep  = fParticlesModel->GetP4(Particles::Type::kBoson, 1);
  TLorentzVector* thad  = fParticlesModel->GetP4(Particles::Type::kParton, 4);
  TLorentzVector* tlep  = fParticlesModel->GetP4(Particles::Type::kParton, 5);

  bhad->SetPxPyPzE(m_bhad_fit_px, m_bhad_fit_py, m_bhad_fit_pz, m_bhad_fit_e);
  blep->SetPxPyPzE(m_blep_fit_px, m_blep_fit_py, m_blep_fit_pz, m_blep_fit_e);
  lq1 ->SetPxPyPzE(m_lq1_fit_px,  m_lq1_fit_py,  m_lq1_fit_pz,  m_lq1_fit_e);
  lq2 ->SetPxPyPzE(m_lq2_fit_px,  m_lq2_fit_py,  m_lq2_fit_pz,  m_lq2_fit_e);
  lep ->SetPxPyPzE(m_lep_fit_px,  m_lep_fit_py,  m_lep_fit_pz,  m_lep_fit_e);
  nu  ->SetPxPyPzE(m_nu_fit_px,   m_nu_fit_py,   m_nu_fit_pz,   m_nu_fit_e);

  (*whad) = (*lq1)  + (*lq2);
  (*wlep) = (*lep)  + (*nu);
  (*thad) = (*whad) + (*bhad);
  (*tlep) = (*wlep) + (*blep);

  // no error
  return 1;
}

// ---------------------------------------------------------
std::vector<double> LikelihoodTopLeptonJets::LogLikelihoodComponents(std::vector<double> parameters) {
  std::vector<double> vecci;

  // calculate 4-vectors
  CalculateLorentzVectors(parameters);

  // temporary flag for a safe use of the transfer functions
  bool TFgoodTmp(true);

  // jet energy resolution terms
  vecci.push_back(m_res_energy_bhad->logp(m_bhad_fit_e, m_bhad_meas_e, &TFgoodTmp));  // comp0
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(m_res_energy_blep->logp(m_blep_fit_e, m_blep_meas_e, &TFgoodTmp));  // comp1
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(m_res_energy_lq1->logp(m_lq1_fit_e, m_lq1_meas_e, &TFgoodTmp));  // comp2
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(m_res_energy_lq2->logp(m_lq2_fit_e, m_lq2_meas_e, &TFgoodTmp));  // comp3
  if (!TFgoodTmp) fTFgood = false;

  // lepton energy resolution terms
  if (m_lepton_type == kElectron) {
    vecci.push_back(m_res_lepton->logp(m_lep_fit_e, m_lep_meas_e, &TFgoodTmp));  // comp4
  } else if (m_lepton_type == kMuon) {
    vecci.push_back(m_res_lepton->logp(m_lep_fit_e* m_lep_meas_sintheta, m_lep_meas_pt, &TFgoodTmp));  // comp4
  }
  if (!TFgoodTmp) fTFgood = false;

  // neutrino px and py
  vecci.push_back(m_res_met->logp(m_nu_fit_px, m_et_miss_x, &TFgoodTmp, m_et_miss_sum));  // comp5
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(m_res_met->logp(m_nu_fit_py, m_et_miss_y, &TFgoodTmp, m_et_miss_sum));  // comp6
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
