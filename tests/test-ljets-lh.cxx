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

#include "KLFitter/DetectorAtlas_8TeV.h"
#include "KLFitter/Fitter.h"
#include "KLFitter/LikelihoodTopLeptonJets.h"
#include "KLFitter/Permutations.h"
#include "TLorentzVector.h"

namespace {
std::unique_ptr<KLFitter::Particles> getExampleParticles(float tag_eff, float tag_ineff) {
  TLorentzVector jet1{};
  jet1.SetPtEtaPhiE(133.56953, 0.2231264, 1.7798618, 137.56292);
  const float jet1_btag_weight{0.6868029};
  const bool jet1_has_btag{false};

  TLorentzVector jet2{};
  jet2.SetPtEtaPhiE(77.834281, 0.8158330, -1.533635, 105.72334);
  const float jet2_btag_weight{-0.869940};
  const bool jet2_has_btag{false};

  TLorentzVector jet3{};
  jet3.SetPtEtaPhiE(49.327293, 1.9828589, -1.878274, 182.64006);
  const float jet3_btag_weight{0.9999086};
  const bool jet3_has_btag{true};

  TLorentzVector jet4{};
  jet4.SetPtEtaPhiE(43.140816, 0.4029131, -0.472721, 47.186804);
  const float jet4_btag_weight{-0.223728};
  const bool jet4_has_btag{false};

  TLorentzVector lep{};
  lep.SetPtEtaPhiE(30.501886, 0.4483959, 2.9649317, 33.620113);

  std::unique_ptr<KLFitter::Particles> particles{new KLFitter::Particles};
  particles->AddParticle(&jet1, jet1.Eta(),
      KLFitter::Particles::kParton, "", 0,
      jet1_has_btag, tag_eff, tag_ineff,
      KLFitter::Particles::kNone, jet1_btag_weight);
  particles->AddParticle(&jet2, jet2.Eta(),
      KLFitter::Particles::kParton, "", 1,
      jet2_has_btag, tag_eff, tag_ineff,
      KLFitter::Particles::kNone, jet2_btag_weight);
  particles->AddParticle(&jet3, jet3.Eta(),
      KLFitter::Particles::kParton, "", 2,
      jet3_has_btag, tag_eff, tag_ineff,
      KLFitter::Particles::kNone, jet3_btag_weight);
  particles->AddParticle(&jet4, jet4.Eta(),
      KLFitter::Particles::kParton, "", 3,
      jet4_has_btag, tag_eff, tag_ineff,
      KLFitter::Particles::kNone, jet4_btag_weight);
  particles->AddParticle(&lep, lep.Eta(), KLFitter::Particles::kMuon, "", 0);
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
  fitter.SetLikelihood(&lh);

  KLFitter::DetectorAtlas_8TeV detector{base_dir + "/data/transferfunctions/8TeV/ttbar/mc12_LCJets_v1"};
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
    std::cout << std::setprecision(8);
    std::cout << "  \tEvtProbability: " << evt_probs.at(perm);
    std::cout << std::endl;
  }
  return 0;
}
