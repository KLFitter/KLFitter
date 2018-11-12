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

#include "KLFitter/LikelihoodTopLeptonJets_Angular.h"

#include <algorithm>
#include <iostream>

#include "BAT/BCMath.h"
#include "KLFitter/ParticleCollection.h"
#include "KLFitter/Permutations.h"
#include "KLFitter/PhysicsConstants.h"
#include "KLFitter/ResolutionBase.h"
#include "TLorentzVector.h"

namespace KLFitter {
// ---------------------------------------------------------
LikelihoodTopLeptonJets_Angular::LikelihoodTopLeptonJets_Angular() = default;

// ---------------------------------------------------------
LikelihoodTopLeptonJets_Angular::~LikelihoodTopLeptonJets_Angular() = default;

// ---------------------------------------------------------
int LikelihoodTopLeptonJets_Angular::AdjustParameterRanges() {
  // adjust limits
  double nsigmas_jet = 7.0;
  double nsigmas_lepton = 2.0;

  double E = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 0)->E();
  double m = fPhysicsConstants.MassBottom();
  if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 0)->M());
  double Emin = std::max(m, E - nsigmas_jet * sqrt(E));
  double Emax = E + nsigmas_jet * sqrt(E);
  SetParameterRange(parBhadE, Emin, Emax);

  E = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 1)->E();
  m = fPhysicsConstants.MassBottom();
  if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 1)->M());
  Emin = std::max(m, E - nsigmas_jet * sqrt(E));
  Emax = E + nsigmas_jet * sqrt(E);
  SetParameterRange(parBlepE, Emin, Emax);

  E = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 2)->E();
  m = 0.001;
  if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 2)->M());
  Emin = std::max(m, E - nsigmas_jet * sqrt(E));
  Emax = E + nsigmas_jet * sqrt(E);
  SetParameterRange(parLQ1E, Emin, Emax);

  E = (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 3)->E();
  m = 0.001;
  if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->GetP4(Particles::Type::kParton, 3)->M());
  Emin = std::max(m, E - nsigmas_jet * sqrt(E));
  Emax = E + nsigmas_jet * sqrt(E);
  SetParameterRange(parLQ2E, Emin, Emax);

  if (m_lepton_type == kElectron) {
    E = (*fParticlesPermuted)->GetP4(Particles::Type::kElectron, 0)->E();
    Emin = std::max(0.001, E - nsigmas_lepton * sqrt(E));
    Emax = E + nsigmas_lepton * sqrt(E);
  } else if (m_lepton_type == kMuon) {
    E = (*fParticlesPermuted)->GetP4(Particles::Type::kMuon, 0)->E();
    double sintheta = sin((*fParticlesPermuted)->GetP4(Particles::Type::kMuon, 0)->Theta());
    double sigrange = nsigmas_lepton * (E * E * sintheta);
    Emin = std::max(0.001, E - sigrange);
    Emax = E + sigrange;
  }
  SetParameterRange(parLepE, Emin, Emax);

  SetParameterRange(parNuPx, m_et_miss_x - 100.0, m_et_miss_x + 100);
  SetParameterRange(parNuPy, m_et_miss_y - 100.0, m_et_miss_y + 100);

  if (m_flag_top_mass_fixed)
    SetParameterRange(parTopM, fPhysicsConstants.MassTop(), fPhysicsConstants.MassTop());

  // no error
  return 1;
}

// ---------------------------------------------------------
double LikelihoodTopLeptonJets_Angular::LogLikelihood(const std::vector<double> & parameters) {
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

  // angular information of leptonic decay

  // create 4-vector for leptonically decaying W boson, charge lepton and corresponding b quark
  TLorentzVector Wlep(m_wlep_fit_px, m_wlep_fit_py, m_wlep_fit_pz, m_wlep_fit_e);
  TLorentzVector Whad(m_whad_fit_px, m_whad_fit_py, m_whad_fit_pz, m_whad_fit_e);
  TLorentzVector lep(m_lep_fit_px,   m_lep_fit_py,  m_lep_fit_pz,  m_lep_fit_e);
  TLorentzVector blep(m_blep_fit_px, m_blep_fit_py, m_blep_fit_pz, m_blep_fit_e);
  TLorentzVector bhad(m_bhad_fit_px, m_bhad_fit_py, m_bhad_fit_pz, m_bhad_fit_e);
  TLorentzVector lq1(m_lq1_fit_px,   m_lq1_fit_py,  m_lq1_fit_pz,  m_lq1_fit_e);
  TLorentzVector lq2(m_lq2_fit_px,   m_lq2_fit_py,  m_lq2_fit_pz,  m_lq2_fit_e);

  // get boost vectors
  TVector3  Wlep_bo(0.0, 0.0, 0.0);
  Wlep_bo = Wlep.BoostVector();

  TVector3  Whad_bo(0.0, 0.0, 0.0);
  Whad_bo = Whad.BoostVector();

  // boost everything into W rest frames
  blep.Boost(-Wlep_bo);
  lep.Boost(-Wlep_bo);

  bhad.Boost(-Whad_bo);
  lq1.Boost(-Whad_bo);
  lq2.Boost(-Whad_bo);

  // calculate 3-vectors
  TVector3 lep3(0.0, 0.0, 0.0);
  TVector3 blep3(0.0, 0.0, 0.0);
  TVector3 bhad3(0.0, 0.0, 0.0);
  TVector3 lq13(0.0, 0.0, 0.0);
  TVector3 lq23(0.0, 0.0, 0.0);

  lep3.SetXYZ(lep.Px(),  lep.Py(),  lep.Pz());
  blep3.SetXYZ(blep.Px(), blep.Py(), blep.Pz());
  bhad3.SetXYZ(bhad.Px(), bhad.Py(), bhad.Pz());
  lq13.SetXYZ(lq1.Px(), lq1.Py(), lq1.Pz());
  lq23.SetXYZ(lq2.Px(), lq2.Py(), lq2.Pz());

  // calculate cos theta *
  double cos_theta     = cos(lep3.Angle(-blep3));
  double cos_theta_had = cos(lq23.Angle(-bhad3));

  // fix helicity fractions
  double F0 = 0.687;
  double FL = 0.311;
  double FR = 0.0017;

  // calculate probability
  double p_angular_lep = (3. / 4. * (1. - cos_theta * cos_theta) * F0 +
                          3. / 8. * (1. - cos_theta) * (1. - cos_theta) * FL +
                          3. / 8. * (1. + cos_theta) * (1. + cos_theta) * FR);

  double p_angular_had = (3. / 4. * (1. - cos_theta_had * cos_theta_had) * F0 +
                          3. / 8. * (1. + cos_theta_had * cos_theta_had) * (FL + FR));

  logprob += log(p_angular_lep);
  logprob += log(p_angular_had);

  // return log of likelihood
  return logprob;
}
}  // namespace KLFitter
