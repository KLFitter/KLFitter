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

#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "KLFitter/DetectorSnowmass.h"
#include "KLFitter/Fitter.h"
#include "KLFitter/LikelihoodTopLeptonJets.h"
#include "KLFitter/Permutations.h"
#include "TLorentzVector.h"

namespace {
std::unique_ptr<KLFitter::ParticleCollection> getExampleParticles(float tag_eff, float tag_ineff) {
  std::unique_ptr<KLFitter::ParticleCollection> particles{new KLFitter::ParticleCollection};

  KLFitter::Particles::Parton parton1{"parton1", TLorentzVector{}};
  parton1.GetP4().SetPtEtaPhiE(133.56953, 0.2231264, 1.7798618, 137.56292);
  parton1.SetDetEta(parton1.GetP4().Eta());
  parton1.SetIdentifier(1);
  parton1.SetBTagWeight(0.6868029);
  parton1.SetIsBTagged(false);
  parton1.SetBTagEfficiency(tag_eff);
  parton1.SetBTagRejection(tag_ineff);
  particles->AddParticle(parton1);

  KLFitter::Particles::Parton parton2{"parton2", TLorentzVector{}};
  parton2.GetP4().SetPtEtaPhiE(77.834281, 0.8158330, -1.533635, 105.72334);
  parton2.SetDetEta(parton2.GetP4().Eta());
  parton2.SetIdentifier(2);
  parton2.SetBTagWeight(-0.869940);
  parton2.SetIsBTagged(false);
  parton2.SetBTagEfficiency(tag_eff);
  parton2.SetBTagRejection(tag_ineff);
  particles->AddParticle(parton2);

  KLFitter::Particles::Parton parton3{"parton3", TLorentzVector{}};
  parton3.GetP4().SetPtEtaPhiE(49.327293, 1.9828589, -1.878274, 182.64006);
  parton3.SetDetEta(parton3.GetP4().Eta());
  parton3.SetIdentifier(3);
  parton3.SetBTagWeight(0.9999086);
  parton3.SetIsBTagged(true);
  parton3.SetBTagEfficiency(tag_eff);
  parton3.SetBTagRejection(tag_ineff);
  particles->AddParticle(parton3);

  KLFitter::Particles::Parton parton4{"parton4", TLorentzVector{}};
  parton4.GetP4().SetPtEtaPhiE(43.140816, 0.4029131, -0.472721, 47.186804);
  parton4.SetDetEta(parton4.GetP4().Eta());
  parton4.SetIdentifier(4);
  parton4.SetBTagWeight(-0.223728);
  parton4.SetIsBTagged(false);
  parton4.SetBTagEfficiency(tag_eff);
  parton4.SetBTagRejection(tag_ineff);
  particles->AddParticle(parton4);

  KLFitter::Particles::Muon muon{"muon", TLorentzVector{}};
  muon.GetP4().SetPtEtaPhiE(30.501886, 0.4483959, 2.9649317, 33.620113);
  muon.SetIdentifier(0);
  particles->AddParticle(muon);

  return particles;
}

void normalizeValues(std::vector<float>* vector) {
  float scale{0};
  for (const auto& i : *vector) { scale += i; }
  for (auto& i : *vector) { i *= 1./scale; }
}
}  // namespace


// ---------------------------------------------------------
// ---------------------------------------------------------

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Wrong number of arguments." << std::endl;
    std::cerr << "Usage: test-ljets-lh [base directory]" << std::endl;
    return -1;
  }
  const auto base_dir = std::string(argv[1]);

  KLFitter::Fitter fitter{};

  // Get one set of example particles. Assume the following
  // efficiencies for the jet b-tagging algorithm:
  //   - 0.7 tag rate
  //   - 1/125 type-II error (false positives)
  const auto particles = getExampleParticles(0.7, 125);
  fitter.SetParticles(particles.get());

  const float met{26.125748};
  const float met_phi{0.3639200};
  fitter.SetET_miss_XY_SumET(met * std::cos(met_phi), met * std::sin(met_phi), met);

  KLFitter::LikelihoodTopLeptonJets lh{};
  lh.SetLeptonType(KLFitter::LikelihoodTopLeptonJets::LeptonType::kMuon);
  lh.SetBTagging(KLFitter::LikelihoodBase::BtaggingMethod::kWorkingPoint);

  if (!fitter.SetLikelihood(&lh)) {
    std::cerr << "Setting up the likelihood failed" << std::endl;
    return -1;
  }

  KLFitter::DetectorSnowmass detector{base_dir + "/data/transferfunctions/snowmass"};
  if (!fitter.SetDetector(&detector)) {
    std::cerr << "Setting up the detector failed" << std::endl;
    return -1;
  }

  const auto nperm = fitter.Permutations()->NPermutations();
  std::vector<float> lh_values{};
  std::vector<float> evt_probs{};
  for (int perm = 0; perm < nperm; ++perm) {
    fitter.Fit(perm);
    lh_values.emplace_back(fitter.Likelihood()->LogLikelihood(fitter.Likelihood()->GetBestFitParameters()));
    evt_probs.emplace_back(std::exp(fitter.Likelihood()->LogEventProbability()));
  }

  normalizeValues(&evt_probs);

  std::cout << std::fixed;  // enforce fixed precision for output
  for (int perm = 0; perm < nperm; ++perm) {
    std::cout << "Permutation: " << perm + 1;
    std::cout << std::setprecision(2);
    std::cout << "  \tLogLikelihood: " << lh_values.at(perm);
    std::cout << std::setprecision(5);
    std::cout << "  \tEvtProbability: " << evt_probs.at(perm);
    std::cout << std::endl;
  }
  return 0;
}
