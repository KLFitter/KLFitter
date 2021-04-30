/*
 * Copyright (c) 2009--2021, the KLFitter developer team
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

#include "KLFitter/LikelihoodOneHadronicTop.h"

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

namespace KLFitter {
// ---------------------------------------------------------
LikelihoodOneHadronicTop::LikelihoodOneHadronicTop()
  : LikelihoodBase::LikelihoodBase()
  , m_flag_top_mass_fixed(false)
  , m_flag_get_par_sigmas_from_TFs(false) {
  // define model particles
  this->DefineModelParticles();

  // define parameters
  this->DefineParameters();
}

// ---------------------------------------------------------
LikelihoodOneHadronicTop::~LikelihoodOneHadronicTop() = default;

// ---------------------------------------------------------
void LikelihoodOneHadronicTop::RequestResolutionFunctions() {
  (*fDetector)->RequestResolutionType(ResolutionType::EnergyLightJet);
  (*fDetector)->RequestResolutionType(ResolutionType::EnergyBJet);
}

// ---------------------------------------------------------
int LikelihoodOneHadronicTop::DefineModelParticles() {
  // create the particles of the model
  fParticlesModel.reset(new Particles{});

  // add model particles
  TLorentzVector dummy{0, 0, 0, 0};
  fParticlesModel->AddParticle(&dummy,
                               Particles::kParton,  // type
                               "hadronic b quark",  // name
                               0,                   // index of corresponding particle
                               Particles::kB);      // b jet (truth)

  fParticlesModel->AddParticle(&dummy,
                               Particles::kParton,
                               "light quark 1",
                               1,                   // index of corresponding particle
                               Particles::kLight);  // light jet (truth)

  fParticlesModel->AddParticle(&dummy,
                               Particles::kParton,
                               "light quark 2",
                               2,                   // index of corresponding particle
                               Particles::kLight);  // light jet (truth)

  fParticlesModel->AddParticle(&dummy, Particles::kBoson, "hadronic W");

  fParticlesModel->AddParticle(&dummy, Particles::kParton, "hadronic top");

  // no error
  return 1;
}

// ---------------------------------------------------------
void LikelihoodOneHadronicTop::DefineParameters() {
  // add parameters of model
  AddParameter("energy hadronic b", fPhysicsConstants.MassBottom(), 1000.0);  // parBhadE
  AddParameter("energy light quark 1", 0.0, 1000.0);                          // parLQ1E
  AddParameter("energy light quark 2", 0.0, 1000.0);                          // parLQ2E
  AddParameter("top mass", 100.0, 1000.0);                                    // parTopM
}

// ---------------------------------------------------------
int LikelihoodOneHadronicTop::CalculateLorentzVectors(std::vector <double> const& parameters) {
  double scale;
  double thad_fit_e;
  double thad_fit_px;
  double thad_fit_py;
  double thad_fit_pz;

  // hadronic b quark
  m_bhad_fit_e = parameters[parBhadE];
  scale = sqrt(m_bhad_fit_e * m_bhad_fit_e - m_bhad_meas_m * m_bhad_meas_m) / m_bhad_meas_p;
  m_bhad_fit_px = scale * m_bhad_meas_px;
  m_bhad_fit_py = scale * m_bhad_meas_py;
  m_bhad_fit_pz = scale * m_bhad_meas_pz;

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

  // hadronic W
  m_whad_fit_e  = m_lq1_fit_e + m_lq2_fit_e;
  m_whad_fit_px = m_lq1_fit_px + m_lq2_fit_px;
  m_whad_fit_py = m_lq1_fit_py + m_lq2_fit_py;
  m_whad_fit_pz = m_lq1_fit_pz + m_lq2_fit_pz;
  m_whad_fit_m = sqrt(m_whad_fit_e * m_whad_fit_e - (m_whad_fit_px * m_whad_fit_px + m_whad_fit_py * m_whad_fit_py + m_whad_fit_pz * m_whad_fit_pz));

  // hadronic top
  thad_fit_e = m_whad_fit_e + m_bhad_fit_e;
  thad_fit_px = m_whad_fit_px + m_bhad_fit_px;
  thad_fit_py = m_whad_fit_py + m_bhad_fit_py;
  thad_fit_pz = m_whad_fit_pz + m_bhad_fit_pz;
  m_thad_fit_m = sqrt(thad_fit_e * thad_fit_e - (thad_fit_px * thad_fit_px + thad_fit_py * thad_fit_py + thad_fit_pz * thad_fit_pz));

  // no error
  return 1;
}

// ---------------------------------------------------------
int LikelihoodOneHadronicTop::RemoveInvariantParticlePermutations() {
  // error code
  int err = 1;

  // remove the permutation from the second and the third jet
  Particles::ParticleType ptype = Particles::kParton;
  std::vector<int> indexVector_Jets;
  indexVector_Jets.push_back(1);
  indexVector_Jets.push_back(2);
  err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Jets);

  // remove invariant jet permutations of all jets not considered
  Particles* particles = (*fPermutations)->Particles();
  indexVector_Jets.clear();
  for (int iPartons = 3; iPartons < particles->NPartons(); iPartons++) {
    indexVector_Jets.push_back(iPartons);
  }
  err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Jets);

  // return error code
  return err;
}

// ---------------------------------------------------------
int LikelihoodOneHadronicTop::AdjustParameterRanges() {
  // adjust limits
  double nsigmas_jet = m_flag_get_par_sigmas_from_TFs ? 10 : 7;

  double E = (*fParticlesPermuted)->Parton(0)->E();
  double m = fPhysicsConstants.MassBottom();
  if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->Parton(0)->M());
  double sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_bhad->GetSigma(E) : sqrt(E);
  double Emin = std::max(m, E - nsigmas_jet * sigma);
  double Emax = E + nsigmas_jet * sigma;
  SetParameterRange(parBhadE, Emin, Emax);

  E = (*fParticlesPermuted)->Parton(1)->E();
  m = 0.001;
  if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->Parton(1)->M());
  sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_lq1->GetSigma(E) : sqrt(E);
  Emin = std::max(m, E - nsigmas_jet * sigma);
  Emax = E + nsigmas_jet * sigma;
  SetParameterRange(parLQ1E, Emin, Emax);

  E = (*fParticlesPermuted)->Parton(2)->E();
  m = 0.001;
  if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->Parton(2)->M());
  sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_lq2->GetSigma(E) : sqrt(E);
  Emin = std::max(m, E - nsigmas_jet * sigma);
  Emax = E + nsigmas_jet * sigma;
  SetParameterRange(parLQ2E, Emin, Emax);

  if (m_flag_top_mass_fixed)
    SetParameterRange(parTopM, fPhysicsConstants.MassTop(), fPhysicsConstants.MassTop());

  // no error
  return 1;
}

// ---------------------------------------------------------
double LikelihoodOneHadronicTop::LogLikelihood(const std::vector<double> & parameters) {
  // calculate 4-vectors
  CalculateLorentzVectors(parameters);

  // define log of likelihood
  double logprob(0.);

  // temporary flag for a safe use of the transfer functions
  bool TFgoodTmp(true);

  // jet energy resolution terms
  logprob += m_res_energy_bhad->logp(m_bhad_fit_e, m_bhad_meas_e, &TFgoodTmp);
  if (!TFgoodTmp) fTFgood = false;

  logprob += m_res_energy_lq1->logp(m_lq1_fit_e, m_lq1_meas_e, &TFgoodTmp);
  if (!TFgoodTmp) fTFgood = false;

  logprob += m_res_energy_lq2->logp(m_lq2_fit_e, m_lq2_meas_e, &TFgoodTmp);
  if (!TFgoodTmp) fTFgood = false;

  // physics constants
  double massW = fPhysicsConstants.MassW();
  double gammaW = fPhysicsConstants.GammaW();
  double gammaTop = fPhysicsConstants.GammaTop();

  // Breit-Wigner of hadronically decaying W-boson
  logprob += BCMath::LogBreitWignerRel(m_whad_fit_m, massW, gammaW);

  // Breit-Wigner of hadronically decaying top quark
  logprob += BCMath::LogBreitWignerRel(m_thad_fit_m, parameters[parTopM], gammaTop);

  // return log of likelihood
  return logprob;
}

// ---------------------------------------------------------
std::vector<double> LikelihoodOneHadronicTop::GetInitialParameters() {
  std::vector<double> values(GetNParameters());

  // energies of the quarks
  values[parBhadE] = m_bhad_meas_e;
  values[parLQ1E]  = m_lq1_meas_e;
  values[parLQ2E]  = m_lq2_meas_e;

  // top mass
  double mtop = (*(*fParticlesPermuted)->Parton(0) + *(*fParticlesPermuted)->Parton(1) + *(*fParticlesPermuted)->Parton(2)).M();
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
int LikelihoodOneHadronicTop::SavePermutedParticles() {
  m_bhad_meas_e      = (*fParticlesPermuted)->Parton(0)->E();
  m_bhad_meas_deteta = (*fParticlesPermuted)->DetEta(0, Particles::kParton);
  m_bhad_meas_px     = (*fParticlesPermuted)->Parton(0)->Px();
  m_bhad_meas_py     = (*fParticlesPermuted)->Parton(0)->Py();
  m_bhad_meas_pz     = (*fParticlesPermuted)->Parton(0)->Pz();
  m_bhad_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(0)->M(), fPhysicsConstants.MassBottom(), &m_bhad_meas_px, &m_bhad_meas_py, &m_bhad_meas_pz, m_bhad_meas_e);
  m_bhad_meas_p      = sqrt(m_bhad_meas_e*m_bhad_meas_e - m_bhad_meas_m*m_bhad_meas_m);

  m_lq1_meas_e      = (*fParticlesPermuted)->Parton(1)->E();
  m_lq1_meas_deteta = (*fParticlesPermuted)->DetEta(1, Particles::kParton);
  m_lq1_meas_px     = (*fParticlesPermuted)->Parton(1)->Px();
  m_lq1_meas_py     = (*fParticlesPermuted)->Parton(1)->Py();
  m_lq1_meas_pz     = (*fParticlesPermuted)->Parton(1)->Pz();
  m_lq1_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(1)->M(), 0., &m_lq1_meas_px, &m_lq1_meas_py, &m_lq1_meas_pz, m_lq1_meas_e);
  m_lq1_meas_p      = sqrt(m_lq1_meas_e*m_lq1_meas_e - m_lq1_meas_m*m_lq1_meas_m);

  m_lq2_meas_e      = (*fParticlesPermuted)->Parton(2)->E();
  m_lq2_meas_deteta = (*fParticlesPermuted)->DetEta(2, Particles::kParton);
  m_lq2_meas_px     = (*fParticlesPermuted)->Parton(2)->Px();
  m_lq2_meas_py     = (*fParticlesPermuted)->Parton(2)->Py();
  m_lq2_meas_pz     = (*fParticlesPermuted)->Parton(2)->Pz();
  m_lq2_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(2)->M(), 0., &m_lq2_meas_px, &m_lq2_meas_py, &m_lq2_meas_pz, m_lq2_meas_e);
  m_lq2_meas_p      = sqrt(m_lq2_meas_e*m_lq2_meas_e - m_lq2_meas_m*m_lq2_meas_m);

  // no error
  return 1;
}

// ---------------------------------------------------------
int LikelihoodOneHadronicTop::SaveResolutionFunctions() {
  m_res_energy_bhad = (*fDetector)->ResEnergyBJet(m_bhad_meas_deteta);
  m_res_energy_lq1  = (*fDetector)->ResEnergyLightJet(m_lq1_meas_deteta);
  m_res_energy_lq2  = (*fDetector)->ResEnergyLightJet(m_lq2_meas_deteta);

  // no error
  return 1;
}

// ---------------------------------------------------------
int LikelihoodOneHadronicTop::BuildModelParticles() {
  if (GetBestFitParameters().size() > 0) CalculateLorentzVectors(GetBestFitParameters());

  TLorentzVector* bhad = fParticlesModel->Parton(0);
  TLorentzVector* lq1  = fParticlesModel->Parton(1);
  TLorentzVector* lq2  = fParticlesModel->Parton(2);

  TLorentzVector* whad = fParticlesModel->Boson(0);
  TLorentzVector* thad = fParticlesModel->Parton(3);

  bhad->SetPxPyPzE(m_bhad_fit_px, m_bhad_fit_py, m_bhad_fit_pz, m_bhad_fit_e);
  lq1 ->SetPxPyPzE(m_lq1_fit_px,  m_lq1_fit_py,  m_lq1_fit_pz,  m_lq1_fit_e);
  lq2 ->SetPxPyPzE(m_lq2_fit_px,  m_lq2_fit_py,  m_lq2_fit_pz,  m_lq2_fit_e);

  (*whad) = (*lq1)  + (*lq2);
  (*thad) = (*whad) + (*bhad);

  // no error
  return 1;
}

// ---------------------------------------------------------
std::vector<double> LikelihoodOneHadronicTop::LogLikelihoodComponents(std::vector<double> parameters) {
  std::vector<double> vecci;

  // calculate 4-vectors
  CalculateLorentzVectors(parameters);

  // temporary flag for a safe use of the transfer functions
  bool TFgoodTmp(true);

  // jet energy resolution terms
  vecci.push_back(m_res_energy_bhad->logp(m_bhad_fit_e, m_bhad_meas_e, &TFgoodTmp));  // comp0
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(m_res_energy_lq1->logp(m_lq1_fit_e, m_lq1_meas_e, &TFgoodTmp));  // comp2
  if (!TFgoodTmp) fTFgood = false;

  vecci.push_back(m_res_energy_lq2->logp(m_lq2_fit_e, m_lq2_meas_e, &TFgoodTmp));  // comp3
  if (!TFgoodTmp) fTFgood = false;

  // physics constants
  double massW = fPhysicsConstants.MassW();
  double gammaW = fPhysicsConstants.GammaW();
  double gammaTop = fPhysicsConstants.GammaTop();

  // Breit-Wigner of hadronically decaying W-boson
  vecci.push_back(BCMath::LogBreitWignerRel(m_whad_fit_m, massW, gammaW));  // comp6

  // Breit-Wigner of hadronically decaying top quark
  vecci.push_back(BCMath::LogBreitWignerRel(m_thad_fit_m, parameters[parTopM], gammaTop));  // comp8

  // return log of likelihood
  return vecci;
}
}  // namespace KLFitter
