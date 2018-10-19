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

#include "KLFitter/ParticleCollection.h"

#include <iostream>

#include "TLorentzVector.h"

namespace KLFitter {
// ---------------------------------------------------------
ParticleCollection::ParticleCollection() = default;

// ---------------------------------------------------------
ParticleCollection::ParticleCollection(const ParticleCollection& o)
    : jets(o.jets)
    , electrons(o.electrons)
    , muons(o.muons)
    , taus(o.taus)
    , neutrinos(o.neutrinos)
    , bosons(o.bosons)
    , photons(o.photons)
    , tracks(o.tracks) {
  // empty
}

// ---------------------------------------------------------
ParticleCollection::~ParticleCollection() = default;

// ---------------------------------------------------------
ParticleCollection& ParticleCollection::operator=(const ParticleCollection& o) {
  jets = o.jets;
  electrons = o.electrons;
  muons = o.muons;
  taus = o.taus;
  neutrinos = o.neutrinos;
  bosons = o.bosons;
  photons = o.photons;
  tracks = o.tracks;

  return *this;
}

// ---------------------------------------------------------
void ParticleCollection::AddParticle(const Particles::Jet& p) {
  if (FindParticle(Particles::Type::kParton, p.GetName())) {
    throw std::invalid_argument("Particle with name " + p.GetName() + " exists already");
  } else {
    jets.emplace_back(p);
  }
}

// ---------------------------------------------------------
void ParticleCollection::AddParticle(const Particles::Electron& p) {
  if (FindParticle(Particles::Type::kElectron, p.GetName())) {
    throw std::invalid_argument("Particle with name " + p.GetName() + " exists already");
  } else {
    electrons.emplace_back(p);
  }
}

// ---------------------------------------------------------
void ParticleCollection::AddParticle(const Particles::Muon& p) {
  if (FindParticle(Particles::Type::kMuon, p.GetName())) {
    throw std::invalid_argument("Particle with name " + p.GetName() + " exists already");
  } else {
    muons.emplace_back(p);
  }
}

// ---------------------------------------------------------
void ParticleCollection::AddParticle(const Particles::Photon& p) {
  if (FindParticle(Particles::Type::kPhoton, p.GetName())) {
    throw std::invalid_argument("Particle with name " + p.GetName() + " exists already");
  } else {
    photons.emplace_back(p);
  }
}

// ---------------------------------------------------------
void ParticleCollection::AddParticle(const Particles::Tau& p) {
  if (FindParticle(Particles::Type::kTau, p.GetName())) {
    throw std::invalid_argument("Particle with name " + p.GetName() + " exists already");
  } else {
    taus.emplace_back(p);
  }
}

// ---------------------------------------------------------
void ParticleCollection::AddParticle(const Particles::Neutrino& p) {
  if (FindParticle(Particles::Type::kNeutrino, p.GetName())) {
    throw std::invalid_argument("Particle with name " + p.GetName() + " exists already");
  } else {
    neutrinos.emplace_back(p);
  }
}

// ---------------------------------------------------------
void ParticleCollection::AddParticle(const Particles::Boson& p) {
  if (FindParticle(Particles::Type::kBoson, p.GetName())) {
    throw std::invalid_argument("Particle with name " + p.GetName() + " exists already");
  } else {
    bosons.emplace_back(p);
  }
}

// ---------------------------------------------------------
void ParticleCollection::AddParticle(const Particles::Track& p) {
  if (FindParticle(Particles::Type::kTrack, p.GetName())) {
    throw std::invalid_argument("Particle with name " + p.GetName() + " exists already");
  } else {
    tracks.emplace_back(p);
  }
}

// ---------------------------------------------------------
int ParticleCollection::AddParticle(const TLorentzVector& particle, double DetEta, float LepCharge, Particles::Type ptype, std::string name, int measuredindex) {
  // check name
  if (name == "") name = Form("particle_%i", NParticles());

  // check if particle with name exists already
  if (!FindParticle(name)) {
    if (ptype == Particles::Type::kParton) {
      Particles::Jet jet{name, particle};
      jet.SetIdentifier(measuredindex);
      jet.SetDetEta(DetEta);
      jets.emplace_back(std::move(jet));
    } else if (ptype == Particles::Type::kElectron) {
      Particles::Electron electron{name, particle};
      electron.SetIdentifier(measuredindex);
      electron.SetDetEta(DetEta);
      electron.SetCharge(LepCharge);
      electrons.emplace_back(std::move(electron));
    } else if (ptype == Particles::Type::kMuon) {
      Particles::Muon muon{name, particle};
      muon.SetIdentifier(measuredindex);
      muon.SetDetEta(DetEta);
      muon.SetCharge(LepCharge);
      muons.emplace_back(std::move(muon));
    } else if (ptype == Particles::Type::kPhoton) {
      Particles::Photon photon{name, particle};
      photon.SetIdentifier(measuredindex);
      photon.SetDetEta(DetEta);
      photons.emplace_back(std::move(photon));
    } else if (ptype == Particles::Type::kTau) {
      taus.emplace_back(name, particle);
    } else if (ptype == Particles::Type::kNeutrino) {
      neutrinos.emplace_back(name, particle);
    } else if (ptype == Particles::Type::kBoson) {
      bosons.emplace_back(name, particle);
    } else if (ptype == Particles::Type::kTrack) {
      tracks.emplace_back(name, particle);
    } else {
      std::cout << "KLFitter::ParticleCollection::AddParticle(). Particle type " << static_cast<std::underlying_type<Particles::Type>::type>(ptype) << " does not exist." << std::endl;
      return 0;
    }
  } else {
    std::cout << "KLFitter::ParticleCollection::AddParticle(). Particle with the name " << name << " exists already." << std::endl;
    return 0;
  }

  if (fabs(particle.P()/particle.E()-1) > 1.e-6 && particle.M() < 0) {  // No Warning if P differs less than 1e-6 from E
    std::cout << "KLFitter::ParticleCollection::AddParticle(). WARNING : A particle with negative mass " << particle.M() << " of type " << static_cast<std::underlying_type<Particles::Type>::type>(ptype) << " was added." << std::endl;
    return 1;
  }

  // no error
  return 1;
}

// --------------------------------------------------------- // THIS FUNCTION IS TO BE REMOVED IN THE NEXT MAJOR RELEASE
int ParticleCollection::AddParticle(const TLorentzVector* const particle, double DetEta, float LepCharge, Particles::Type ptype, std::string name, int measuredindex) {
  return AddParticle(*particle, DetEta, LepCharge, ptype, name, measuredindex);
}

// ---------------------------------------------------------
  int ParticleCollection::AddParticle(const TLorentzVector& particle, double DetEta, Particles::Type ptype, std::string name, int measuredindex, bool isBtagged, double bTagEff, double bTagRej, Particles::JetTrueFlavor trueflav, double btagweight) {
  // check name
  if (name == "") name = Form("particle_%i", NParticles());

  // check if particle with name exists already
  if (!FindParticle(name)) {
    if (ptype == Particles::Type::kParton) {
      Particles::Jet jet{name, particle, isBtagged};
      jet.SetIdentifier(measuredindex);
      jet.SetDetEta(DetEta);
      jet.SetBTagEfficiency(bTagEff);
      jet.SetBTagRejection(bTagRej);
      jet.SetTrueFlavor(trueflav);
      jet.SetBTagWeight(btagweight);
      jets.emplace_back(std::move(jet));
    } else if (ptype == Particles::Type::kElectron) {
      Particles::Electron electron{name, particle};
      electron.SetIdentifier(measuredindex);
      electron.SetDetEta(DetEta);
      electrons.emplace_back(std::move(electron));
    } else if (ptype == Particles::Type::kMuon) {
      Particles::Muon muon{name, particle};
      muon.SetIdentifier(measuredindex);
      muon.SetDetEta(DetEta);
      muons.emplace_back(std::move(muon));
    } else if (ptype == Particles::Type::kPhoton) {
      Particles::Photon photon{name, particle};
      photon.SetIdentifier(measuredindex);
      photon.SetDetEta(DetEta);
      photons.emplace_back(std::move(photon));
    } else if (ptype == Particles::Type::kTau) {
      taus.emplace_back(name, particle);
    } else if (ptype == Particles::Type::kNeutrino) {
      neutrinos.emplace_back(name, particle);
    } else if (ptype == Particles::Type::kBoson) {
      bosons.emplace_back(name, particle);
    } else if (ptype == Particles::Type::kTrack) {
      Particles::Track track{name, particle};
      track.SetIdentifier(measuredindex);
      tracks.emplace_back(std::move(track));
    } else {
      std::cout << "KLFitter::ParticleCollection::AddParticle(). Particle type " << static_cast<std::underlying_type<Particles::Type>::type>(ptype) << " does not exist." << std::endl;
      return 0;
    }
  } else {
    std::cout << "KLFitter::ParticleCollection::AddParticle(). Particle with the name " << name << " exists already." << std::endl;
    return 0;
  }

  if (fabs(particle.P()/particle.E()-1) > 1.e-6 && particle.M() < 0) {  // No Warning if P differs less than 1e-6 from E
    std::cout << "KLFitter::ParticleCollection::AddParticle(). WARNING : A particle with negative mass " << particle.M() << " of type " << static_cast<std::underlying_type<Particles::Type>::type>(ptype) << " was added." << std::endl;
    return 1;
  }

  // no error
  return 1;
}

// --------------------------------------------------------- // THIS FUNCTION IS TO BE REMOVED IN THE NEXT MAJOR RELEASE
int ParticleCollection::AddParticle(const TLorentzVector* const particle, double DetEta, Particles::Type ptype, std::string name, int measuredindex, bool isBtagged, double bTagEff, double bTagRej, Particles::JetTrueFlavor trueflav, double btagweight) {
  return AddParticle(*particle, DetEta, ptype, name, measuredindex, isBtagged, bTagEff, bTagRej, trueflav, btagweight);
}

// ---------------------------------------------------------
int ParticleCollection::AddParticle(const TLorentzVector& particle, Particles::Type ptype, std::string name, int measuredindex, bool isBtagged, double bTagEff, double bTagRej, Particles::JetTrueFlavor trueflav, double btagweight) {
  return AddParticle(particle, -999, ptype, name, measuredindex, isBtagged, bTagEff, bTagRej, trueflav, btagweight);
}

// --------------------------------------------------------- // THIS FUNCTION IS TO BE REMOVED IN THE NEXT MAJOR RELEASE
int ParticleCollection::AddParticle(const TLorentzVector* const particle, Particles::Type ptype, std::string name, int measuredindex, bool isBtagged, double bTagEff, double bTagRej, Particles::JetTrueFlavor trueflav, double btagweight) {
  return AddParticle(*particle, ptype, name, measuredindex, isBtagged, bTagEff, bTagRej, trueflav, btagweight);
}

// ---------------------------------------------------------
int ParticleCollection::AddParticle(const TLorentzVector& particle, Particles::Type ptype, std::string name, int measuredindex, Particles::JetTrueFlavor trueflav, double btagweight) {
  return AddParticle(particle, -999, ptype, name, measuredindex, false, -1., -1., trueflav, btagweight);
}

// --------------------------------------------------------- // THIS FUNCTION IS TO BE REMOVED IN THE NEXT MAJOR RELEASE
int ParticleCollection::AddParticle(const TLorentzVector* const particle, Particles::Type ptype, std::string name, int measuredindex, Particles::JetTrueFlavor trueflav, double btagweight) {
  return AddParticle(*particle, ptype, name, measuredindex, trueflav, btagweight);
}

// ---------------------------------------------------------
int ParticleCollection::AddParticle(const TLorentzVector& particle, Particles::Type ptype, std::string name, int measuredindex, const std::vector<double>& uncertainies) {
  if (ptype == Particles::Type::kTrack) {
    Particles::Track track{name, particle};
    track.SetIdentifier(measuredindex);
    track.SetUncertainties(uncertainies);
    tracks.emplace_back(std::move(track));
    return 1;
  }

  // Return with error
  return 0;
}

// ---------------------------------------------------------
void ParticleCollection::RemoveParticle(Particles::Type ptype, size_t index) {
  if (ptype == Particles::Type::kParton) {
    jets.erase(jets.begin() + index);
  } else if (ptype == Particles::Type::kElectron) {
    electrons.erase(electrons.begin() + index);
  } else if (ptype == Particles::Type::kMuon) {
    muons.erase(muons.begin() + index);
  } else if (ptype == Particles::Type::kPhoton) {
    photons.erase(photons.begin() + index);
  } else if (ptype == Particles::Type::kTau) {
    taus.erase(taus.begin() + index);
  } else if (ptype == Particles::Type::kNeutrino) {
    neutrinos.erase(neutrinos.begin() + index);
  } else if (ptype == Particles::Type::kBoson) {
    bosons.erase(bosons.begin() + index);
  } else if (ptype == Particles::Type::kTrack) {
    tracks.erase(tracks.begin() + index);
  }
}

// ---------------------------------------------------------
const Particles::Base* ParticleCollection::FindParticle(const std::string& name) const {
  // loop over all jets
  for (auto jet = jets.begin(); jet != jets.end(); ++jet) {
    if (name != jet->GetName()) continue;
    return &*jet;
  }

  // loop over all electrons
  for (auto el = electrons.begin(); el != electrons.end(); ++el) {
    if (name != el->GetName()) continue;
    return &*el;
  }

  // loop over all muons
  for (auto mu = muons.begin(); mu != muons.end(); ++mu) {
    if (name != mu->GetName()) continue;
    return &*mu;
  }

  // loop over all taus
  for (auto tau = taus.begin(); tau != taus.end(); ++tau) {
    if (name != tau->GetName()) continue;
    return &*tau;
  }

  // loop over all neutrinos
  for (auto neutrino = neutrinos.begin(); neutrino != neutrinos.end(); ++neutrino) {
    if (name != neutrino->GetName()) continue;
    return &*neutrino;
  }

  // loop over all bosons
  for (auto boson = bosons.begin(); boson != bosons.end(); ++boson) {
    if (name != boson->GetName()) continue;
    return &*boson;
  }

  // loop over all photons
  for (auto ph = photons.begin(); ph != photons.end(); ++ph) {
    if (name != ph->GetName()) continue;
    return &*ph;
  }

  // loop over all tracks
  for (auto track = tracks.begin(); track != tracks.end(); ++track) {
    if (name != track->GetName()) continue;
    return &*track;
  }

  // particle not found
  return nullptr;
}

// ---------------------------------------------------------
const Particles::Base* ParticleCollection::FindParticle(Particles::Type ptype, const std::string& name) const {
  if (ptype == Particles::Type::kParton) {
    for (auto jet = jets.begin(); jet != jets.end(); ++jet) {
      if (name != jet->GetName()) continue;
      return &*jet;
    }
  } else if (ptype == Particles::Type::kElectron) {
    for (auto el = electrons.begin(); el != electrons.end(); ++el) {
      if (name != el->GetName()) continue;
      return &*el;
    }
  } else if (ptype == Particles::Type::kMuon) {
    for (auto mu = muons.begin(); mu != muons.end(); ++mu) {
      if (name != mu->GetName()) continue;
      return &*mu;
    }
  } else if (ptype == Particles::Type::kTau) {
    for (auto tau = taus.begin(); tau != taus.end(); ++tau) {
      if (name != tau->GetName()) continue;
      return &*tau;
    }
  } else if (ptype == Particles::Type::kNeutrino) {
    for (auto neutrino = neutrinos.begin(); neutrino != neutrinos.end(); ++neutrino) {
      if (name != neutrino->GetName()) continue;
      return &*neutrino;
    }
  } else if (ptype == Particles::Type::kBoson) {
    for (auto boson = bosons.begin(); boson != bosons.end(); ++boson) {
      if (name != boson->GetName()) continue;
      return &*boson;
    }
  } else if (ptype == Particles::Type::kPhoton) {
    for (auto ph = photons.begin(); ph != photons.end(); ++ph) {
      if (name != ph->GetName()) continue;
      return &*ph;
    }
  } else if (ptype == Particles::Type::kTrack) {
    for (auto track = tracks.begin(); track != tracks.end(); ++track) {
      if (name != track->GetName()) continue;
      return &*track;
    }
  }

  // particle not found
  return nullptr;
}

// ---------------------------------------------------------
const TLorentzVector* ParticleCollection::GetP4(Particles::Type ptype, size_t index) const {
  if (ptype == Particles::Type::kBoson) {
    return &bosons.at(index).GetP4();
  } else if (ptype == Particles::Type::kElectron) {
    return &electrons.at(index).GetP4();
  } else if (ptype == Particles::Type::kMuon) {
    return &muons.at(index).GetP4();
  } else if (ptype == Particles::Type::kNeutrino) {
    return &neutrinos.at(index).GetP4();
  } else if (ptype == Particles::Type::kPhoton) {
    return &photons.at(index).GetP4();
  } else if (ptype == Particles::Type::kTau) {
    return &taus.at(index).GetP4();
  } else if (ptype == Particles::Type::kTrack) {
    return &tracks.at(index).GetP4();
  } else {
    return nullptr;
  }
}

// ---------------------------------------------------------
TLorentzVector* ParticleCollection::GetP4(Particles::Type ptype, size_t index) {
  if (ptype == Particles::Type::kBoson) {
    return &bosons.at(index).GetP4();
  } else if (ptype == Particles::Type::kElectron) {
    return &electrons.at(index).GetP4();
  } else if (ptype == Particles::Type::kMuon) {
    return &muons.at(index).GetP4();
  } else if (ptype == Particles::Type::kNeutrino) {
    return &neutrinos.at(index).GetP4();
  } else if (ptype == Particles::Type::kPhoton) {
      return &photons.at(index).GetP4();
  } else if (ptype == Particles::Type::kTau) {
    return &taus.at(index).GetP4();
  } else if (ptype == Particles::Type::kTrack) {
    return &tracks.at(index).GetP4();
  } else {
    return nullptr;
  }
}

// ---------------------------------------------------------
int ParticleCollection::NParticles(Particles::Type ptype) const {
  if (ptype == Particles::Type::kParton) {
    return static_cast<int>(jets.size());
  } else if (ptype == Particles::Type::kElectron) {
    return static_cast<int>(electrons.size());
  } else if (ptype == Particles::Type::kMuon) {
    return static_cast<int>(muons.size());
  } else if (ptype == Particles::Type::kPhoton) {
    return static_cast<int>(photons.size());
  } else if (ptype == Particles::Type::kTau) {
    return static_cast<int>(taus.size());
  } else if (ptype == Particles::Type::kNeutrino) {
    return static_cast<int>(neutrinos.size());
  } else if (ptype == Particles::Type::kBoson) {
    return static_cast<int>(bosons.size());
  } else if (ptype == Particles::Type::kTrack) {
    return static_cast<int>(tracks.size());
  }
  return 0;
}

// ---------------------------------------------------------
int ParticleCollection::NBTags() const {
  int sum{0};
  for (const auto& jet : jets) {
    if (jet.GetIsBTagged()) sum++;
  }

  return sum;
}
}  // namespace KLFitter
