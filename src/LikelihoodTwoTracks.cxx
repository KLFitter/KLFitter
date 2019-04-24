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

#include "KLFitter/LikelihoodTwoTracks.h"

#include <iostream>
#include <algorithm>

#include "BAT/BCMath.h"
#include "BAT/BCParameter.h"
#include "KLFitter/DetectorBase.h"
#include "KLFitter/ParticleCollection.h"
#include "KLFitter/Permutations.h"
#include "KLFitter/PhysicsConstants.h"
#include "KLFitter/ResolutionBase.h"

namespace KLFitter {
// ---------------------------------------------------------
LikelihoodTwoTracks::LikelihoodTwoTracks()
    : LikelihoodBase::LikelihoodBase()
    , m_pion_mass(139.57)
    , m_kshort_mass(497.6)
    , m_kshort_width(7.351e-12) {
  // define model particles
  this->DefineModelParticles();

  // define parameters
  this->DefineParameters();

  fFlagIntegrate = 0;
}

// ---------------------------------------------------------
LikelihoodTwoTracks::~LikelihoodTwoTracks() = default;

// ---------------------------------------------------------
double LikelihoodTwoTracks::Log3DGaus(double x0, double x1, double x2,
                                      double mean0, double mean1, double mean2,
                                      double sigma00, double sigma10, double sigma11,
                                      double sigma20, double sigma21, double sigma22) {
  // the likelihood of a multivariate gaussian is \ln L = -\frac{1}{2} \left( \ln (|\boldsymbol\Sigma|\,) + (\mathbf{x}-\boldsymbol\mu)^{\rm T}\boldsymbol\Sigma^{-1}(\mathbf{x}-\boldsymbol\mu) + k\ln(2\pi) \right)

  // determinant of the correlation matrix, where sigmaAB is the
  // element in the Ath row and Bth column
  double det = sigma00 * sigma11 * sigma22 + 2 * (sigma10 * sigma21 * sigma20) - sigma11 * sigma20 * sigma20 - sigma00 * sigma21 * sigma21 - sigma22 * sigma10 * sigma10;

  // the matrix product of the vector (x-mean)^T*Sigma^-1*(x-mean)
  double inverse = (1. / det) * (
                                 (x0 - mean0) * (x0 - mean0) * (sigma11 * sigma22 - sigma21 * sigma21) +
                                 2 * (x0 - mean0) * (x1 - mean1) * (sigma20 * sigma21 - sigma10 * sigma22) +
                                 2 * (x0 - mean0) * (x2 - mean2) * (sigma10 * sigma21 - sigma11 * sigma20) +
                                 (x1 - mean1) * (x1 - mean1) * (sigma00 * sigma22 - sigma20 * sigma20) +
                                 2 * (x1 - mean1) * (x2 - mean2) * (sigma10 * sigma20 - sigma00 * sigma21) +
                                 (x2 - mean2) * (x2 - mean2) * (sigma00 * sigma11 - sigma10 * sigma10)
                                 );

  // return the entire likelihood
  return -0.5 * (log(det) + inverse + 3 + log(2 * M_PI));
}

// ---------------------------------------------------------
int LikelihoodTwoTracks::DefineModelParticles() {
  // create the particles of the model
  fParticlesModel.reset(new ParticleCollection{});

  // add model particles
  Particles::Boson kshort{"Kshort", TLorentzVector{}};
  kshort.SetIdentifier(0);
  fParticlesModel->AddParticle(kshort);

  Particles::Track pi_plus{"pi plus", TLorentzVector{}};
  pi_plus.SetIdentifier(0);
  fParticlesModel->AddParticle(pi_plus);

  Particles::Track pi_minus{"pi minus", TLorentzVector{}};
  pi_minus.SetIdentifier(0);
  fParticlesModel->AddParticle(pi_minus);

  // no error
  return 1;
}

// ---------------------------------------------------------
void LikelihoodTwoTracks::DefineParameters() {
  // add parameters of model
  AddParameter("phi pi plus",    -M_PI, M_PI);          // parPiPlusPhi
  AddParameter("theta pi plus",  0.0, M_PI);            // parPiPlusTheta
  AddParameter("p pi plus",      0.0, 100000.0);        // parPiPlusP
  AddParameter("phi pi minus",   -M_PI, M_PI);          // parPiMinusPhi
  AddParameter("theta pi minus", 0.0, M_PI);            // parPiMinusTheta
  AddParameter("p pi minus",     0.0, 100000.0);        // parPiMinusP
  AddParameter("m Kshort",       m_kshort_mass-250., m_kshort_mass+250.);  // parKShortM
}

// ---------------------------------------------------------
int LikelihoodTwoTracks::CalculateLorentzVectors(const std::vector<double>& parameters) {
  TLorentzVector t1 = TLorentzVector();
  TLorentzVector t2 = TLorentzVector();
  TLorentzVector Ks = TLorentzVector();

  t1.SetPtEtaPhiM(sin(parameters[parPiPlusTheta]) * parameters[parPiPlusP], -log(tan(parameters[parPiPlusTheta] / 2.)), parameters[parPiPlusPhi], m_pion_mass);

  t2.SetPtEtaPhiM(sin(parameters[parPiMinusTheta]) * parameters[parPiMinusP], -log(tan(parameters[parPiMinusTheta] / 2.)), parameters[parPiMinusPhi], m_pion_mass);

  Ks = t1 + t2;

  m_t1_fit_phi = parameters[parPiPlusPhi];
  m_t1_fit_theta = parameters[parPiPlusTheta];
  m_t1_fit_p = parameters[parPiPlusP];
  m_t1_fit_m = m_pion_mass;
  m_t2_fit_phi = parameters[parPiMinusPhi];
  m_t2_fit_theta = parameters[parPiMinusTheta];
  m_t2_fit_p = parameters[parPiMinusP];
  m_t2_fit_m = m_pion_mass;
  m_ks_fit_m = Ks.M();

  // no error
  return 1;
}

// ---------------------------------------------------------
int LikelihoodTwoTracks::Initialize() {
  // error code
  int err = 1;

  // save the current permuted particles
  err *= SavePermutedParticles();

  // set initial values
  // (only for Markov chains - initial parameters for other minimisation methods are set in Fitter.cxx)
  SetInitialParameters( GetInitialParameters() );

  // return error code
  return err;
}

// ---------------------------------------------------------
int LikelihoodTwoTracks::RemoveInvariantParticlePermutations() {
  // error code
  int err = 1;

  // remove the permutation from the second and the third jet
  Particles::Type ptype = Particles::Type::kTrack;
  err *= (*fPermutations)->InvariantParticlePermutations(ptype, std::vector<int>{0, 1});

  // return error code
  return err;
}

// ---------------------------------------------------------
double LikelihoodTwoTracks::LogLikelihood(const std::vector<double>& parameters) {
  // calculate 4-vectors
  CalculateLorentzVectors(parameters);

  // define log of likelihood
  double logprob(0.);

  logprob += Log3DGaus(m_t1_fit_phi, m_t1_fit_theta, 1/m_t1_fit_p, m_t1_meas_phi, m_t1_meas_theta, 1/m_t1_meas_p, m_t1_meas_sigma00, m_t1_meas_sigma10, m_t1_meas_sigma11, m_t1_meas_sigma20, m_t1_meas_sigma21, m_t1_meas_sigma22);

  logprob += Log3DGaus(m_t2_fit_phi, m_t2_fit_theta, 1/m_t2_fit_p, m_t2_meas_phi, m_t2_meas_theta, 1/m_t2_meas_p, m_t2_meas_sigma00, m_t2_meas_sigma10, m_t2_meas_sigma11, m_t2_meas_sigma20, m_t2_meas_sigma21, m_t2_meas_sigma22);

  logprob += BCMath::LogBreitWignerRel(m_ks_fit_m, m_kshort_mass, m_kshort_width);

  // return log of likelihood
  return logprob;
}

// ---------------------------------------------------------
std::vector<double> LikelihoodTwoTracks::GetInitialParameters() {
  std::vector<double> values(GetNParameters());

  TLorentzVector t1 = TLorentzVector();
  TLorentzVector t2 = TLorentzVector();
  TLorentzVector Ks = TLorentzVector();

  t1.SetPtEtaPhiM(sin(m_t1_meas_theta)*m_t1_fit_p, -log(tan(m_t1_meas_theta/2.)), m_t1_meas_phi, m_pion_mass);

  t2.SetPtEtaPhiM(sin(m_t2_meas_theta)*m_t1_fit_p, -log(tan(m_t2_meas_theta/2.)), m_t2_meas_phi, m_pion_mass);

  Ks = t1 + t2;

  values[parPiPlusPhi] = m_t1_meas_phi;
  values[parPiPlusTheta] = m_t1_meas_theta;
  values[parPiPlusP]  = m_t1_meas_p;
  values[parPiMinusPhi]  = m_t2_meas_phi;
  values[parPiMinusTheta] = m_t2_meas_theta;
  values[parPiMinusP]  = m_t2_meas_p;
  values[parKShortM] = Ks.M();

  // return the vector
  return values;
}

// ---------------------------------------------------------
int LikelihoodTwoTracks::SavePermutedParticles() {
  m_t1_meas_phi = (*fParticlesPermuted)->GetP4(Particles::Type::kTrack, 0)->Phi();
  m_t1_meas_theta = (*fParticlesPermuted)->GetP4(Particles::Type::kTrack, 0)->Theta();
  m_t1_meas_p = (*fParticlesPermuted)->GetP4(Particles::Type::kTrack, 0)->P();

  m_t1_meas_sigma00 = (*fParticlesPermuted)->tracks.at(0).GetUncertainties().at(0);
  m_t1_meas_sigma10 = (*fParticlesPermuted)->tracks.at(0).GetUncertainties().at(1);
  m_t1_meas_sigma11 = (*fParticlesPermuted)->tracks.at(0).GetUncertainties().at(2);
  m_t1_meas_sigma20 = (*fParticlesPermuted)->tracks.at(0).GetUncertainties().at(3);
  m_t1_meas_sigma21 = (*fParticlesPermuted)->tracks.at(0).GetUncertainties().at(4);
  m_t1_meas_sigma22 = (*fParticlesPermuted)->tracks.at(0).GetUncertainties().at(5);

  m_t2_meas_phi = (*fParticlesPermuted)->GetP4(Particles::Type::kTrack, 1)->Phi();
  m_t2_meas_theta = (*fParticlesPermuted)->GetP4(Particles::Type::kTrack, 1)->Theta();
  m_t2_meas_p = (*fParticlesPermuted)->GetP4(Particles::Type::kTrack, 1)->P();

  m_t2_meas_sigma00 = (*fParticlesPermuted)->tracks.at(1).GetUncertainties().at(0);
  m_t2_meas_sigma10 = (*fParticlesPermuted)->tracks.at(1).GetUncertainties().at(1);
  m_t2_meas_sigma11 = (*fParticlesPermuted)->tracks.at(1).GetUncertainties().at(2);
  m_t2_meas_sigma20 = (*fParticlesPermuted)->tracks.at(1).GetUncertainties().at(3);
  m_t2_meas_sigma21 = (*fParticlesPermuted)->tracks.at(1).GetUncertainties().at(4);
  m_t2_meas_sigma22 = (*fParticlesPermuted)->tracks.at(1).GetUncertainties().at(5);


  // no error
  return 1;
}

// ---------------------------------------------------------
int LikelihoodTwoTracks::BuildModelParticles() {
  if (GetBestFitParameters().size() > 0) CalculateLorentzVectors(GetBestFitParameters());

  TLorentzVector * t1 = fParticlesModel->GetP4(Particles::Type::kTrack, 0);
  TLorentzVector * t2 = fParticlesModel->GetP4(Particles::Type::kTrack, 1);
  TLorentzVector * Ks = fParticlesModel->GetP4(Particles::Type::kBoson, 0);


  t1->SetPtEtaPhiM(sin(m_t1_fit_theta)*m_t1_fit_p, -log(tan(m_t1_fit_theta/2.)), m_t1_fit_phi, m_t1_fit_m);

  t2->SetPtEtaPhiM(sin(m_t2_fit_theta)*m_t1_fit_p, -log(tan(m_t2_fit_theta/2.)), m_t2_fit_phi, m_t2_fit_m);

  (*Ks) = (*t1) + (*t2);

  // no error
  return 1;
}

// ---------------------------------------------------------
std::vector<double> LikelihoodTwoTracks::LogLikelihoodComponents(std::vector<double> parameters) {
  std::vector<double> vecci;

  // calculate 4-vectors
  CalculateLorentzVectors(parameters);

  vecci.push_back(Log3DGaus(m_t1_fit_phi, m_t1_fit_theta, m_t1_fit_p, m_t1_meas_phi, m_t1_meas_theta, m_t1_meas_p, m_t1_meas_sigma00, m_t1_meas_sigma10, m_t1_meas_sigma11, m_t1_meas_sigma20, m_t1_meas_sigma21, m_t1_meas_sigma22)); //comp0

  vecci.push_back(Log3DGaus(m_t2_fit_phi, m_t2_fit_theta, m_t2_fit_p, m_t2_meas_phi, m_t2_meas_theta, m_t2_meas_p, m_t2_meas_sigma00, m_t2_meas_sigma10, m_t2_meas_sigma11, m_t2_meas_sigma20, m_t2_meas_sigma21, m_t2_meas_sigma22)); //comp1

  vecci.push_back(BCMath::LogBreitWignerRel(m_ks_fit_m, m_kshort_mass, m_kshort_width));  //comp2

  // return log of likelihood
  return vecci;
}
}  // namespace KLFitter
