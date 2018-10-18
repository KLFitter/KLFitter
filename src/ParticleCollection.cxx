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
int ParticleCollection::AddParticle(const TLorentzVector& particle, double DetEta, float LepCharge, Particle::Type ptype, std::string name, int measuredindex) {
  // check name
  if (name == "") name = Form("particle_%i", NParticles());

  // get index and type
  TLorentzVector* vect = 0;
  int index = 0;
  Particle::Type temptype = Particle::Type::kParton;

  // check if particle with name exists already
  if (!FindParticle(name, vect, &index, &temptype)) {
    if (ptype == Particle::Type::kParton) {
      Particle::Jet jet{name, particle};
      jet.SetIdentifier(measuredindex);
      jet.SetDetEta(DetEta);
      jets.emplace_back(std::move(jet));
    } else if (ptype == Particle::Type::kElectron) {
      Particle::Electron electron{name, particle};
      electron.SetIdentifier(measuredindex);
      electron.SetDetEta(DetEta);
      electron.SetCharge(LepCharge);
      electrons.emplace_back(std::move(electron));
    } else if (ptype == Particle::Type::kMuon) {
      Particle::Muon muon{name, particle};
      muon.SetIdentifier(measuredindex);
      muon.SetDetEta(DetEta);
      muon.SetCharge(LepCharge);
      muons.emplace_back(std::move(muon));
    } else if (ptype == Particle::Type::kPhoton) {
      Particle::Photon photon{name, particle};
      photon.SetIdentifier(measuredindex);
      photon.SetDetEta(DetEta);
      photons.emplace_back(std::move(photon));
    } else if (ptype == Particle::Type::kTau) {
      taus.emplace_back(name, particle);
    } else if (ptype == Particle::Type::kNeutrino) {
      neutrinos.emplace_back(name, particle);
    } else if (ptype == Particle::Type::kBoson) {
      bosons.emplace_back(name, particle);
    } else if (ptype == Particle::Type::kTrack) {
      tracks.emplace_back(name, particle);
    } else {
      std::cout << "KLFitter::ParticleCollection::AddParticle(). Particle type " << static_cast<std::underlying_type<Particle::Type>::type>(ptype) << " does not exist." << std::endl;
      return 0;
    }
  } else {
    std::cout << "KLFitter::ParticleCollection::AddParticle(). Particle with the name " << name << " exists already." << std::endl;
    return 0;
  }

  if (fabs(particle.P()/particle.E()-1) > 1.e-6 && particle.M() < 0) {  // No Warning if P differs less than 1e-6 from E
    std::cout << "KLFitter::ParticleCollection::AddParticle(). WARNING : A particle with negative mass " << particle.M() << " of type " << static_cast<std::underlying_type<Particle::Type>::type>(ptype) << " was added." << std::endl;
    return 1;
  }

  // no error
  return 1;
}

// --------------------------------------------------------- // THIS FUNCTION IS TO BE REMOVED IN THE NEXT MAJOR RELEASE
int ParticleCollection::AddParticle(const TLorentzVector* const particle, double DetEta, float LepCharge, Particle::Type ptype, std::string name, int measuredindex) {
  return AddParticle(*particle, DetEta, LepCharge, ptype, name, measuredindex);
}

// ---------------------------------------------------------
  int ParticleCollection::AddParticle(const TLorentzVector& particle, double DetEta, Particle::Type ptype, std::string name, int measuredindex, bool isBtagged, double bTagEff, double bTagRej, Particle::JetTrueFlavor trueflav, double btagweight) {
  // check name
  if (name == "") name = Form("particle_%i", NParticles());

  // get index and type
  TLorentzVector* vect = 0;
  int index = 0;
  Particle::Type temptype = Particle::Type::kParton;

  // check if particle with name exists already
  if (!FindParticle(name, vect, &index, &temptype)) {
    if (ptype == Particle::Type::kParton) {
      Particle::Jet jet{name, particle, isBtagged};
      jet.SetIdentifier(measuredindex);
      jet.SetDetEta(DetEta);
      jet.SetBTagEfficiency(bTagEff);
      jet.SetBTagRejection(bTagRej);
      jet.SetTrueFlavor(static_cast<Particle::JetTrueFlavor>(trueflav));
      jet.SetBTagWeight(btagweight);
      jets.emplace_back(std::move(jet));
    } else if (ptype == Particle::Type::kElectron) {
      Particle::Electron electron{name, particle};
      electron.SetIdentifier(measuredindex);
      electron.SetDetEta(DetEta);
      electrons.emplace_back(std::move(electron));
    } else if (ptype == Particle::Type::kMuon) {
      Particle::Muon muon{name, particle};
      muon.SetIdentifier(measuredindex);
      muon.SetDetEta(DetEta);
      muons.emplace_back(std::move(muon));
    } else if (ptype == Particle::Type::kPhoton) {
      Particle::Photon photon{name, particle};
      photon.SetIdentifier(measuredindex);
      photon.SetDetEta(DetEta);
      photons.emplace_back(std::move(photon));
    } else if (ptype == Particle::Type::kTau) {
      taus.emplace_back(name, particle);
    } else if (ptype == Particle::Type::kNeutrino) {
      neutrinos.emplace_back(name, particle);
    } else if (ptype == Particle::Type::kBoson) {
      bosons.emplace_back(name, particle);
    } else if (ptype == Particle::Type::kTrack) {
      Particle::Track track{name, particle};
      track.SetIdentifier(measuredindex);
      tracks.emplace_back(std::move(track));
    } else {
      std::cout << "KLFitter::ParticleCollection::AddParticle(). Particle type " << static_cast<std::underlying_type<Particle::Type>::type>(ptype) << " does not exist." << std::endl;
      return 0;
    }
  } else {
    std::cout << "KLFitter::ParticleCollection::AddParticle(). Particle with the name " << name << " exists already." << std::endl;
    return 0;
  }

  if (fabs(particle.P()/particle.E()-1) > 1.e-6 && particle.M() < 0) {  // No Warning if P differs less than 1e-6 from E
    std::cout << "KLFitter::ParticleCollection::AddParticle(). WARNING : A particle with negative mass " << particle.M() << " of type " << static_cast<std::underlying_type<Particle::Type>::type>(ptype) << " was added." << std::endl;
    return 1;
  }

  // no error
  return 1;
}

// --------------------------------------------------------- // THIS FUNCTION IS TO BE REMOVED IN THE NEXT MAJOR RELEASE
int ParticleCollection::AddParticle(const TLorentzVector* const particle, double DetEta, Particle::Type ptype, std::string name, int measuredindex, bool isBtagged, double bTagEff, double bTagRej, Particle::JetTrueFlavor trueflav, double btagweight) {
  return AddParticle(*particle, DetEta, ptype, name, measuredindex, isBtagged, bTagEff, bTagRej, trueflav, btagweight);
}

// ---------------------------------------------------------
int ParticleCollection::AddParticle(const TLorentzVector& particle, Particle::Type ptype, std::string name, int measuredindex, bool isBtagged, double bTagEff, double bTagRej, Particle::JetTrueFlavor trueflav, double btagweight) {
  return AddParticle(particle, -999, ptype, name, measuredindex, isBtagged, bTagEff, bTagRej, trueflav, btagweight);
}

// --------------------------------------------------------- // THIS FUNCTION IS TO BE REMOVED IN THE NEXT MAJOR RELEASE
int ParticleCollection::AddParticle(const TLorentzVector* const particle, Particle::Type ptype, std::string name, int measuredindex, bool isBtagged, double bTagEff, double bTagRej, Particle::JetTrueFlavor trueflav, double btagweight) {
  return AddParticle(*particle, ptype, name, measuredindex, isBtagged, bTagEff, bTagRej, trueflav, btagweight);
}

// ---------------------------------------------------------
int ParticleCollection::AddParticle(const TLorentzVector& particle, Particle::Type ptype, std::string name, int measuredindex, Particle::JetTrueFlavor trueflav, double btagweight) {
  return AddParticle(particle, -999, ptype, name, measuredindex, false, -1., -1., trueflav, btagweight);
}

// --------------------------------------------------------- // THIS FUNCTION IS TO BE REMOVED IN THE NEXT MAJOR RELEASE
int ParticleCollection::AddParticle(const TLorentzVector* const particle, Particle::Type ptype, std::string name, int measuredindex, Particle::JetTrueFlavor trueflav, double btagweight) {
  return AddParticle(*particle, ptype, name, measuredindex, trueflav, btagweight);
}

// ---------------------------------------------------------
int ParticleCollection::AddParticle(const TLorentzVector& particle, Particle::Type ptype, std::string name, int measuredindex, const std::vector<double>& uncertainies) {
  if (ptype == Particle::Type::kTrack) {
    Particle::Track track{name, particle};
    track.SetIdentifier(measuredindex);
    track.SetUncertainties(uncertainies);
    tracks.emplace_back(std::move(track));
    return 1;
  }

  // Return with error
  return 0;
}

// ---------------------------------------------------------
int ParticleCollection::RemoveParticle(int index, Particle::Type ptype) {
  if (ptype == Particle::Type::kParton) {
    jets.erase(jets.begin() + index);
  } else if (ptype == Particle::Type::kElectron) {
    electrons.erase(electrons.begin() + index);
  } else if (ptype == Particle::Type::kMuon) {
    muons.erase(muons.begin() + index);
  } else if (ptype == Particle::Type::kPhoton) {
    photons.erase(photons.begin() + index);
  } else if (ptype == Particle::Type::kTau) {
    taus.erase(taus.begin() + index);
  } else if (ptype == Particle::Type::kNeutrino) {
    neutrinos.erase(neutrinos.begin() + index);
  } else if (ptype == Particle::Type::kBoson) {
    bosons.erase(bosons.begin() + index);
  } else if (ptype == Particle::Type::kTrack) {
    tracks.erase(tracks.begin() + index);
  }

  // no error
  std::cout << "KLFitter::ParticleCollection::RemoveParticle(). Particle type " << static_cast<std::underlying_type<Particle::Type>::type>(ptype) << " does not exist." << std::endl;
  return 1;
}

// ---------------------------------------------------------
int ParticleCollection::RemoveParticle(const std::string& name) {
  // get index and type
  TLorentzVector* vect = 0;
  int index = 0;
  Particle::Type ptype = Particle::Type::kParton;

  // remove particle
  if (FindParticle(name, vect, &index, &ptype)) {
    return RemoveParticle(index, ptype);
  } else {
    std::cout << "KLFitter::ParticleCollection::RemoveParticles(). Could not find particle with name " << name << "." << std::endl;
    return 0;
  }
}

// ---------------------------------------------------------
TLorentzVector* ParticleCollection::Particle(const std::string& name) {
  TLorentzVector* particle = 0;
  int index = 0;
  Particle::Type ptype = Particle::Type::kParton;

  // find particle
  if (!FindParticle(name, particle, &index, &ptype)) {
    std::cout << "KLFitter::ParticleCollection::Particle(). Could not find particles." << std::endl;
    return 0;
  }

  // return 4-vector
  return particle;
}

// ---------------------------------------------------------
TLorentzVector* ParticleCollection::Particle(int index, Particle::Type ptype) {
  if (ptype == Particle::Type::kParton) {
    return &jets.at(index).GetP4();
  } else if (ptype == Particle::Type::kElectron) {
    return &electrons.at(index).GetP4();
  } else if (ptype == Particle::Type::kMuon) {
    return &muons.at(index).GetP4();
  } else if (ptype == Particle::Type::kPhoton) {
    return &photons.at(index).GetP4();
  } else if (ptype == Particle::Type::kTau) {
    return &taus.at(index).GetP4();
  } else if (ptype == Particle::Type::kNeutrino) {
    return &neutrinos.at(index).GetP4();
  } else if (ptype == Particle::Type::kBoson) {
    return &bosons.at(index).GetP4();
  } else if (ptype == Particle::Type::kTrack) {
    return &tracks.at(index).GetP4();
  }

  // Return nullptr
  std::cout << "KLFitter::ParticleCollection::Particle(). Particle type " << static_cast<std::underlying_type<Particle::Type>::type>(ptype) << " does not exist." << std::endl;
  return nullptr;
}

// ---------------------------------------------------------
int ParticleCollection::FindParticle(const std::string& name, TLorentzVector* &particle, int *index, Particle::Type *ptype) {
  // loop over all jets
  for (auto jet = jets.begin(); jet != jets.end(); ++jet) {
    if (name != jet->GetName()) continue;
    particle = &jet->GetP4();
    *index = jet - jets.begin();
    *ptype = Particle::Type::kParton;
    return 1;
  }

  // loop over all electrons
  for (auto el = electrons.begin(); el != electrons.end(); ++el) {
    if (name != el->GetName()) continue;
    particle = &el->GetP4();
    *index = el - electrons.begin();
    *ptype = Particle::Type::kElectron;
    return 1;
  }

  // loop over all muons
  for (auto mu = muons.begin(); mu != muons.end(); ++mu) {
    if (name != mu->GetName()) continue;
    particle = &mu->GetP4();
    *index = mu - muons.begin();
    *ptype = Particle::Type::kMuon;
    return 1;
  }

  // loop over all taus
  for (auto tau = taus.begin(); tau != taus.end(); ++tau) {
    if (name != tau->GetName()) continue;
    particle = &tau->GetP4();
    *index = tau - taus.begin();
    *ptype = Particle::Type::kTau;
    return 1;
  }

  // loop over all neutrinos
  for (auto neutrino = neutrinos.begin(); neutrino != neutrinos.end(); ++neutrino) {
    if (name != neutrino->GetName()) continue;
    particle = &neutrino->GetP4();
    *index = neutrino - neutrinos.begin();
    *ptype = Particle::Type::kNeutrino;
    return 1;
  }

  // loop over all bosons
  for (auto boson = bosons.begin(); boson != bosons.end(); ++boson) {
    if (name != boson->GetName()) continue;
    particle = &boson->GetP4();
    *index = boson - bosons.begin();
    *ptype = Particle::Type::kBoson;
    return 1;
  }

  // loop over all photons
  for (auto ph = photons.begin(); ph != photons.end(); ++ph) {
    if (name != ph->GetName()) continue;
    particle = &ph->GetP4();
    *index = ph - photons.begin();
    *ptype = Particle::Type::kPhoton;
    return 1;
  }

  // loop over all tracks
  for (auto track = tracks.begin(); track != tracks.end(); ++track) {
    if (name != track->GetName()) continue;
    particle = &track->GetP4();
    *index = track - tracks.begin();
    *ptype = Particle::Type::kTrack;
    return 1;
  }

  // particle not found
  return 0;
}

// ---------------------------------------------------------
TLorentzVector* ParticleCollection::Parton(int index) {
  return &jets.at(index).GetP4();
}

// ---------------------------------------------------------
TLorentzVector* ParticleCollection::Electron(int index) {
  return &electrons.at(index).GetP4();
}

// ---------------------------------------------------------
TLorentzVector* ParticleCollection::Muon(int index) {
  return &muons.at(index).GetP4();
}

// ---------------------------------------------------------
TLorentzVector* ParticleCollection::Tau(int index) {
  return &taus.at(index).GetP4();
}

// ---------------------------------------------------------
TLorentzVector* ParticleCollection::Boson(int index) {
  return &bosons.at(index).GetP4();
}

// ---------------------------------------------------------
TLorentzVector* ParticleCollection::Neutrino(int index) {
  return &neutrinos.at(index).GetP4();
}

// ---------------------------------------------------------
TLorentzVector* ParticleCollection::Photon(int index) {
  return &photons.at(index).GetP4();
}

// ---------------------------------------------------------
TLorentzVector* ParticleCollection::Track(int index) {
  return &tracks.at(index).GetP4();
}

// ---------------------------------------------------------
int ParticleCollection::NParticles(Particle::Type ptype) const {
  if (ptype == Particle::Type::kParton) {
    return static_cast<int>(jets.size());
  } else if (ptype == Particle::Type::kElectron) {
    return static_cast<int>(electrons.size());
  } else if (ptype == Particle::Type::kMuon) {
    return static_cast<int>(muons.size());
  } else if (ptype == Particle::Type::kPhoton) {
    return static_cast<int>(photons.size());
  } else if (ptype == Particle::Type::kTau) {
    return static_cast<int>(taus.size());
  } else if (ptype == Particle::Type::kNeutrino) {
    return static_cast<int>(neutrinos.size());
  } else if (ptype == Particle::Type::kBoson) {
    return static_cast<int>(bosons.size());
  } else if (ptype == Particle::Type::kTrack) {
    return static_cast<int>(tracks.size());
  }
  return 0;
}

// ---------------------------------------------------------
const std::vector<double>* ParticleCollection::Uncertainties(int index, Particle::Type ptype) const {
  if (ptype == Particle::Type::kTrack) return &tracks.at(index).GetUncertainties();

  // return error value
  std::cout << "KLFitter::ParticleCollection::Uncertainties(). Particle type " << static_cast<std::underlying_type<Particle::Type>::type>(ptype) << " does not store uncertainties." << std::endl;
  return nullptr;
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
