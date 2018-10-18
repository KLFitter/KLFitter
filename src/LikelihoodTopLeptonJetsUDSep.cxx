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

#include "KLFitter/LikelihoodTopLeptonJetsUDSep.h"

#include <iostream>
#include <algorithm>

#include "BAT/BCMath.h"
#include "BAT/BCParameter.h"
#include "KLFitter/DetectorBase.h"
#include "KLFitter/ParticleCollection.h"
#include "KLFitter/Permutations.h"
#include "KLFitter/PhysicsConstants.h"
#include "KLFitter/ResolutionBase.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"

namespace KLFitter {
// ---------------------------------------------------------
LikelihoodTopLeptonJetsUDSep::LikelihoodTopLeptonJetsUDSep()
  : LikelihoodTopLeptonJets::LikelihoodTopLeptonJets()
  , m_ljet_separation_method(LikelihoodTopLeptonJetsUDSep::kNone) {
  // define model particles
  this->DefineModelParticles();

  // define parameters
  this->DefineParameters();
  }

// ---------------------------------------------------------
LikelihoodTopLeptonJetsUDSep::~LikelihoodTopLeptonJetsUDSep() = default;

// ---------------------------------------------------------
int LikelihoodTopLeptonJetsUDSep::DefineModelParticles() {
  // create the particles of the model
  fParticlesModel.reset(new ParticleCollection{});

  // add model particles
  TLorentzVector dummy{0, 0, 0, 0};
  fParticlesModel->AddParticle(&dummy,
                               ParticleCollection::kParton,  // type
                               "hadronic b quark",  // name
                               0,                   // index of corresponding particle
                               Particle::JetTrueFlavor::kB);      // b jet (truth)

  fParticlesModel->AddParticle(&dummy,
                               ParticleCollection::kParton,
                               "leptonic b quark",
                               1,                   // index of corresponding particle
                               Particle::JetTrueFlavor::kB);      // b jet (truth)

  fParticlesModel->AddParticle(&dummy,
                               ParticleCollection::kParton,
                               "light up type quark",
                               2,                      // index of corresponding particle
                               Particle::JetTrueFlavor::kLightUp);   // light up type jet (truth)

  fParticlesModel->AddParticle(&dummy,
                               ParticleCollection::kParton,
                               "light down type quark",
                               3,                        // index of corresponding particle
                               Particle::JetTrueFlavor::kLightDown);   // light down type jet (truth)

  if (m_lepton_type == kElectron) {
    fParticlesModel->AddParticle(&dummy, ParticleCollection::kElectron, "electron");
  } else if (m_lepton_type == kMuon) {
    fParticlesModel->AddParticle(&dummy, ParticleCollection::kMuon, "muon");
  }

  fParticlesModel->AddParticle(&dummy, ParticleCollection::kNeutrino, "neutrino");

  fParticlesModel->AddParticle(&dummy, ParticleCollection::kBoson, "hadronic W");
  fParticlesModel->AddParticle(&dummy, ParticleCollection::kBoson, "leptonic W");

  fParticlesModel->AddParticle(&dummy, ParticleCollection::kParton, "hadronic top");
  fParticlesModel->AddParticle(&dummy, ParticleCollection::kParton, "leptonic top");

  // no error
  return 1;
}

// ---------------------------------------------------------
void LikelihoodTopLeptonJetsUDSep::DefineParameters() {
  // rename light quark parameters
  this->GetParameter("energy light quark 1")->SetName("energy light up type quark");
  this->GetParameter("energy light quark 2")->SetName("energy light down type quark");
}

// ---------------------------------------------------------
int LikelihoodTopLeptonJetsUDSep::RemoveInvariantParticlePermutations() {
  // error code
  int err = 1;

  ParticleCollection::ParticleType ptype = ParticleCollection::kParton;
  std::vector<int> indexVector_Jets;
  // remove invariant jet permutations of all jets not considered
  ParticleCollection* particles = (*fPermutations)->Particles();
  indexVector_Jets.clear();
  for (int iPartons = 4; iPartons < particles->NPartons(); iPartons++) {
    indexVector_Jets.push_back(iPartons);
  }
  err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Jets);

  // remove the permutation from the other lepton
  if (m_lepton_type == kElectron) {
    ptype = ParticleCollection::kMuon;
    std::vector<int> indexVector_Muons;
    for (int iMuon = 0; iMuon < particles->NMuons(); iMuon++) {
      indexVector_Muons.push_back(iMuon);
    }
    err *= (*fPermutations)->InvariantParticlePermutations(ptype, indexVector_Muons);
  } else if (m_lepton_type == kMuon) {
    ptype = ParticleCollection::kElectron;
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
double LikelihoodTopLeptonJetsUDSep::LogEventProbability() {
  double logprob = 0;

  if (fBTagMethod != kNotag) {
    double logprobbtag = LogEventProbabilityBTag();
    if (logprobbtag <= -1e99) return -1e99;
    logprob += logprobbtag;
  }
  if (m_ljet_separation_method != kNone) {
    double logprobljetweight = LogEventProbabilityLJetReweight();
    if (logprobljetweight <= -1e99) return -1e99;
    logprob += logprobljetweight;
  }

  // use integrated value of LogLikelihood (default)
  if (fFlagIntegrate) {
    logprob += log(GetIntegral());
  } else {
    logprob += LogLikelihood(GetBestFitParameters());
  }

  return logprob;
}

// ---------------------------------------------------------
double LikelihoodTopLeptonJetsUDSep::LogEventProbabilityLJetReweight() {
  double logprob = 0;
  switch (m_ljet_separation_method) {
  case kPermReweight:

    if (!(m_up_jet_pt_histo && m_down_jet_pt_histo&& m_bjet_pt_histo && m_up_jet_tag_weight_histo && m_down_jet_tag_weight_histo && m_bjet_tag_weight_histo)) {
      std::cout <<  " KLFitter::LikelihoodTopLeptonJetsUDSep::LogEventProbabilityLJetReweight() : Histograms were not set properly! " << std::endl;
      return -1e99;
    }

    for (int i = 0; i < fParticlesModel->NPartons(); ++i) {
      // get index of corresponding measured particle.

      int index = fParticlesModel->JetIndex(i);

      if (index < 0) {
        continue;
      }
      if (!((*fParticlesPermuted)->BTagWeightSet(index))) {
        std::cout <<  " KLFitter::LikelihoodTopLeptonJetsUDSep::LogEventProbabilityLJetReweight() : bTag weight for particle was not set ! " << std::endl;
        return -1e99;
      }
      Particle::JetTrueFlavor trueFlavor = fParticlesModel->TrueFlavor(i);
      if (trueFlavor == Particle::JetTrueFlavor::kB) {
        logprob += log(BJetPt((*fParticlesPermuted)->Parton(index)->Pt()));
        logprob += log(BJetTagWeight((*fParticlesPermuted)->BTagWeight(index)));
      }
      if (trueFlavor == Particle::JetTrueFlavor::kLightUp) {
        logprob += log(UpJetPt((*fParticlesPermuted)->Parton(index)->Pt()));
        logprob += log(UpJetTagWeight((*fParticlesPermuted)->BTagWeight(index)));
      }
      if (trueFlavor == Particle::JetTrueFlavor::kLightDown) {
        logprob += log(DownJetPt((*fParticlesPermuted)->Parton(index)->Pt()));
        logprob += log(DownJetTagWeight((*fParticlesPermuted)->BTagWeight(index)));
      }
    }
    return logprob;
    break;

  case kPermReweight2D:
    if (!(m_up_jet_2d_weight_histo && m_down_jet_2d_weight_histo && m_bjet_2d_weight_histo)) {
      std::cout <<  " KLFitter::LikelihoodTopLeptonJetsUDSep::LogEventProbabilityLJetReweight() : 2D Histograms were not set properly! " << std::endl;
      return -1e99;
    }

    for (int i = 0; i < fParticlesModel->NPartons(); ++i) {
      // get index of corresponding measured particle.

      int index = fParticlesModel->JetIndex(i);

      if (index < 0) {
        continue;
      }
      if (!((*fParticlesPermuted)->BTagWeightSet(index))) {
        std::cout <<  " KLFitter::LikelihoodTopLeptonJetsUDSep::LogEventProbabilityLJetReweight() : bTag weight for particle was not set ! " << std::endl;
        return -1e99;
      }
      Particle::JetTrueFlavor trueFlavor = fParticlesModel->TrueFlavor(i);
      if (trueFlavor == Particle::JetTrueFlavor::kB) {
        logprob += log(BJetProb((*fParticlesPermuted)->BTagWeight(index), (*fParticlesPermuted)->Parton(index)->Pt()));
      }
      if (trueFlavor == Particle::JetTrueFlavor::kLightUp) {
        logprob += log(UpJetProb((*fParticlesPermuted)->BTagWeight(index), (*fParticlesPermuted)->Parton(index)->Pt()));
      }
      if (trueFlavor == Particle::JetTrueFlavor::kLightDown) {
        logprob += log(DownJetProb((*fParticlesPermuted)->BTagWeight(index), (*fParticlesPermuted)->Parton(index)->Pt()));
      }
    }
    return logprob;
    break;

  default:
    return logprob;
    break;
  }
}

// ---------------------------------------------------------
double LikelihoodTopLeptonJetsUDSep::LogEventProbabilityBTag() {
  double logprob = 0;

  double probbtag = 1;

  if (fBTagMethod == kVeto) {
    // loop over all model particles.  calculate the overall b-tagging
    // probability which is the product of all probabilities.
    for (int i = 0; i < fParticlesModel->NPartons(); ++i) {
      // get index of corresponding measured particle.
      int index = fParticlesModel->JetIndex(i);
      if (index < 0)
        continue;

      Particle::JetTrueFlavor trueFlavor = fParticlesModel->TrueFlavor(i);
      bool isBTagged = fParticlesModel->IsBTagged(i);
      if (((trueFlavor == Particle::JetTrueFlavor::kLightUp) || (trueFlavor == Particle::JetTrueFlavor::kLightDown)) && isBTagged == true)
        probbtag = 0.;
    }

    if (probbtag > 0) {
      logprob += log(probbtag);
    } else {
      return -1e99;
    }
  } else if (fBTagMethod == kWorkingPoint) {
    // loop over all model particles.  calculate the overall b-tagging
    // probability which is the product of all probabilities.
    for (int i = 0; i < fParticlesModel->NPartons(); ++i) {
      // get index of corresponding measured particle.
      int index = fParticlesModel->JetIndex(i);
      if (index < 0)
        continue;

      Particle::JetTrueFlavor trueFlavor = fParticlesModel->TrueFlavor(i);
      bool isBTagged = fParticlesModel->IsBTagged(i);
      double efficiency = fParticlesModel->BTaggingEfficiency(i);
      double rejection = fParticlesModel->BTaggingRejection(i);
      if (rejection < 0 || efficiency < 0) {
        std::cout <<  " KLFitter::LikelihoodTopLeptonJetsUDSep::LogEventProbability() : Your working points are not set properly! Returning 0 probability " << std::endl;
        return -1e99;
      }

      if (((trueFlavor == Particle::JetTrueFlavor::kLightUp) || (trueFlavor == Particle::JetTrueFlavor::kLightDown)) && isBTagged) {
        logprob += log(1./rejection);
      } else if (((trueFlavor == Particle::JetTrueFlavor::kLightUp) || (trueFlavor == Particle::JetTrueFlavor::kLightDown)) && !isBTagged) {
        logprob += log(1 - 1./rejection);
      } else if (trueFlavor == Particle::JetTrueFlavor::kB && isBTagged) {
        logprob += log(efficiency);
      } else if (trueFlavor == Particle::JetTrueFlavor::kB && !isBTagged) {
        logprob += log(1 - efficiency);
      } else {
        std::cout << " KLFitter::LikelihoodTopLeptonJetsUDSep::LogEventProbability() : b-tagging association failed! " << std::endl;
      }
    }
  }

  return logprob;
}

// ---------------------------------------------------------
double LikelihoodTopLeptonJetsUDSep::UpJetPt(double pt) {
  return m_up_jet_pt_histo->GetBinContent(m_up_jet_pt_histo->GetXaxis()->FindBin(pt));
}

// ---------------------------------------------------------
double LikelihoodTopLeptonJetsUDSep::DownJetPt(double pt) {
  return m_down_jet_pt_histo->GetBinContent(m_down_jet_pt_histo->GetXaxis()->FindBin(pt));
}

// ---------------------------------------------------------
double LikelihoodTopLeptonJetsUDSep::BJetPt(double pt) {
  return m_bjet_pt_histo->GetBinContent(m_bjet_pt_histo->GetXaxis()->FindBin(pt));
}

// ---------------------------------------------------------
double LikelihoodTopLeptonJetsUDSep::UpJetTagWeight(double tagweight) {
  return m_up_jet_tag_weight_histo->GetBinContent(m_up_jet_tag_weight_histo->GetXaxis()->FindBin(tagweight));
}

// ---------------------------------------------------------
double LikelihoodTopLeptonJetsUDSep::DownJetTagWeight(double tagweight) {
  return m_down_jet_tag_weight_histo->GetBinContent(m_down_jet_tag_weight_histo->GetXaxis()->FindBin(tagweight));
}

// ---------------------------------------------------------
double LikelihoodTopLeptonJetsUDSep::BJetTagWeight(double tagweight) {
  return m_bjet_tag_weight_histo->GetBinContent(m_bjet_tag_weight_histo->GetXaxis()->FindBin(tagweight));
}

// ---------------------------------------------------------
double LikelihoodTopLeptonJetsUDSep::UpJetProb(double tagweight, double pt) {
  return m_up_jet_2d_weight_histo->GetBinContent(m_up_jet_2d_weight_histo->GetXaxis()->FindBin(tagweight), m_up_jet_2d_weight_histo->GetYaxis()->FindBin(pt));
}

// ---------------------------------------------------------
double LikelihoodTopLeptonJetsUDSep::DownJetProb(double tagweight, double pt) {
  return m_down_jet_2d_weight_histo->GetBinContent(m_down_jet_2d_weight_histo->GetXaxis()->FindBin(tagweight), m_down_jet_2d_weight_histo->GetYaxis()->FindBin(pt));
}

// ---------------------------------------------------------
double LikelihoodTopLeptonJetsUDSep::BJetProb(double tagweight, double pt) {
  return m_bjet_2d_weight_histo->GetBinContent(m_bjet_2d_weight_histo->GetXaxis()->FindBin(tagweight), m_bjet_2d_weight_histo->GetYaxis()->FindBin(pt));
}

// ---------------------------------------------------------
int LikelihoodTopLeptonJetsUDSep::LHInvariantPermutationPartner(int iperm, int nperms, int *switchpar1, int *switchpar2) {
  int partnerid = -1;
  int cache = iperm % 6;
  switch (nperms) {
  case 24:
    if (iperm % 2) {
      partnerid = iperm - 1;
    } else {
      partnerid = iperm + 1;
    }
    break;

  case 120:
    if (cache > 2) {
      partnerid = iperm - 3;
    } else {
      partnerid = iperm + 3;
    }
    break;

  default: partnerid = -1;
  }
  *switchpar1 = 2;
  *switchpar2 = 3;
  return partnerid;
}
}  // namespace KLFitter
