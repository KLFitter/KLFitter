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

#include "KLFitter/Particles.h"

#include <iostream>

#include "TLorentzVector.h"

namespace KLFitter {
// ---------------------------------------------------------
Particles::Particles() = default;

// ---------------------------------------------------------
Particles::Particles(const Particles& o) = default;

// ---------------------------------------------------------
Particles::~Particles() = default;

// ---------------------------------------------------------
Particles& Particles::operator=(const Particles& o) = default;

// ---------------------------------------------------------
int Particles::AddParticle(const TLorentzVector& particle, double DetEta, float LepCharge, Particles::ParticleType ptype, std::string name, int measuredindex) {
  // check name
  if (name == "") name = Form("particle_%i", NParticles());

  // get index and type
  TLorentzVector* vect = 0;
  int index = 0;
  Particles::ParticleType temptype = kParton;

  // check if particle with name exists already
  if (!FindParticle(name, vect, &index, &temptype)) {
    if (ptype == Particles::kParton) {
      Particle::Jet jet{name, particle};
      jet.SetIdentifier(measuredindex);
      jet.SetDetEta(DetEta);
      m_jets.emplace_back(std::move(jet));
    } else if (ptype == Particles::kElectron) {
      Particle::Electron electron{name, particle};
      electron.SetIdentifier(measuredindex);
      electron.SetDetEta(DetEta);
      electron.SetCharge(LepCharge);
      m_electrons.emplace_back(std::move(electron));
    } else if (ptype == Particles::kMuon) {
      Particle::Muon muon{name, particle};
      muon.SetIdentifier(measuredindex);
      muon.SetDetEta(DetEta);
      muon.SetCharge(LepCharge);
      m_muons.emplace_back(std::move(muon));
    } else if (ptype == Particles::kPhoton) {
      Particle::Photon photon{name, particle};
      photon.SetIdentifier(measuredindex);
      photon.SetDetEta(DetEta);
      m_photons.emplace_back(std::move(photon));
    } else if (ptype == Particles::kTau) {
      m_taus.emplace_back(name, particle);
    } else if (ptype == Particles::kNeutrino) {
      m_neutrinos.emplace_back(name, particle);
    } else if (ptype == Particles::kBoson) {
      m_bosons.emplace_back(name, particle);
    } else if (ptype == Particles::kTrack) {
      m_tracks.emplace_back(name, particle);
    } else {
      std::cout << "KLFitter::Particles::AddParticle(). Particle type " << ptype << " does not exist." << std::endl;
      return 0;
    }
  } else {
    std::cout << "KLFitter::Particles::AddParticle(). Particle with the name " << name << " exists already." << std::endl;
    return 0;
  }

  if (fabs(particle.P()/particle.E()-1) > 1.e-6 && particle.M() < 0) {  // No Warning if P differs less than 1e-6 from E
    std::cout << "KLFitter::Particles::AddParticle(). WARNING : A particle with negative mass " << particle.M() << " of type " << ptype << " was added." << std::endl;
    return 1;
  }

  // no error
  return 1;
}

// --------------------------------------------------------- // THIS FUNCTION IS TO BE REMOVED IN THE NEXT MAJOR RELEASE
int Particles::AddParticle(const TLorentzVector* const particle, double DetEta, float LepCharge, Particles::ParticleType ptype, std::string name, int measuredindex) {
  return AddParticle(*particle, DetEta, LepCharge, ptype, name, measuredindex);
}

// ---------------------------------------------------------
int Particles::AddParticle(const TLorentzVector& particle, double DetEta, Particles::ParticleType ptype, std::string name, int measuredindex, bool isBtagged, double bTagEff, double bTagRej, TrueFlavorType trueflav, double btagweight) {
  // check name
  if (name == "") name = Form("particle_%i", NParticles());

  // get index and type
  TLorentzVector* vect = 0;
  int index = 0;
  Particles::ParticleType temptype = kParton;

  // check if particle with name exists already
  if (!FindParticle(name, vect, &index, &temptype)) {
    if (ptype == Particles::kParton) {
      Particle::Jet jet{name, particle, isBtagged};
      jet.SetIdentifier(measuredindex);
      jet.SetDetEta(DetEta);
      jet.SetBTagEfficiency(bTagEff);
      jet.SetBTagRejection(bTagRej);
      jet.SetTrueFlavor(static_cast<Particle::JetTrueFlavor>(trueflav));
      jet.SetBTagWeight(btagweight);
      m_jets.emplace_back(std::move(jet));
    } else if (ptype == Particles::kElectron) {
      Particle::Electron electron{name, particle};
      electron.SetIdentifier(measuredindex);
      electron.SetDetEta(DetEta);
      m_electrons.emplace_back(std::move(electron));
    } else if (ptype == Particles::kMuon) {
      Particle::Muon muon{name, particle};
      muon.SetIdentifier(measuredindex);
      muon.SetDetEta(DetEta);
      m_muons.emplace_back(std::move(muon));
    } else if (ptype == Particles::kPhoton) {
      Particle::Photon photon{name, particle};
      photon.SetIdentifier(measuredindex);
      photon.SetDetEta(DetEta);
      m_photons.emplace_back(std::move(photon));
    } else if (ptype == Particles::kTau) {
      m_taus.emplace_back(name, particle);
    } else if (ptype == Particles::kNeutrino) {
      m_neutrinos.emplace_back(name, particle);
    } else if (ptype == Particles::kBoson) {
      m_bosons.emplace_back(name, particle);
    } else if (ptype == Particles::kTrack) {
      Particle::Track track{name, particle};
      track.SetIdentifier(measuredindex);
      m_tracks.emplace_back(std::move(track));
    } else {
      std::cout << "KLFitter::Particles::AddParticle(). Particle type " << ptype << " does not exist." << std::endl;
      return 0;
    }
  } else {
    std::cout << "KLFitter::Particles::AddParticle(). Particle with the name " << name << " exists already." << std::endl;
    return 0;
  }

  if (fabs(particle.P()/particle.E()-1) > 1.e-6 && particle.M() < 0) {  // No Warning if P differs less than 1e-6 from E
    std::cout << "KLFitter::Particles::AddParticle(). WARNING : A particle with negative mass " << particle.M() << " of type " << ptype << " was added." << std::endl;
    return 1;
  }

  // no error
  return 1;
}

// --------------------------------------------------------- // THIS FUNCTION IS TO BE REMOVED IN THE NEXT MAJOR RELEASE
int Particles::AddParticle(const TLorentzVector* const particle, double DetEta, Particles::ParticleType ptype, std::string name, int measuredindex, bool isBtagged, double bTagEff, double bTagRej, TrueFlavorType trueflav, double btagweight) {
  return AddParticle(*particle, DetEta, ptype, name, measuredindex, isBtagged, bTagEff, bTagRej, trueflav, btagweight);
}

// ---------------------------------------------------------
int Particles::AddParticle(const TLorentzVector& particle, Particles::ParticleType ptype, std::string name, int measuredindex, bool isBtagged, double bTagEff, double bTagRej, TrueFlavorType trueflav, double btagweight) {
  return AddParticle(particle, -999, ptype, name, measuredindex, isBtagged, bTagEff, bTagRej, trueflav, btagweight);
}

// --------------------------------------------------------- // THIS FUNCTION IS TO BE REMOVED IN THE NEXT MAJOR RELEASE
int Particles::AddParticle(const TLorentzVector* const particle, Particles::ParticleType ptype, std::string name, int measuredindex, bool isBtagged, double bTagEff, double bTagRej, TrueFlavorType trueflav, double btagweight) {
  return AddParticle(*particle, ptype, name, measuredindex, isBtagged, bTagEff, bTagRej, trueflav, btagweight);
}

// ---------------------------------------------------------
int Particles::AddParticle(const TLorentzVector& particle, Particles::ParticleType ptype, std::string name, int measuredindex, TrueFlavorType trueflav, double btagweight) {
  return AddParticle(particle, -999, ptype, name, measuredindex, false, -1., -1., trueflav, btagweight);
}

// --------------------------------------------------------- // THIS FUNCTION IS TO BE REMOVED IN THE NEXT MAJOR RELEASE
int Particles::AddParticle(const TLorentzVector* const particle, Particles::ParticleType ptype, std::string name, int measuredindex, TrueFlavorType trueflav, double btagweight) {
  return AddParticle(*particle, ptype, name, measuredindex, trueflav, btagweight);
}

// ---------------------------------------------------------
int Particles::AddParticle(const TLorentzVector& particle, Particles::ParticleType ptype, std::string name, int measuredindex, const std::vector<double>& uncertainies) {
  if (ptype == Particles::kTrack) {
    Particle::Track track{name, particle};
    track.SetIdentifier(measuredindex);
    track.SetUncertainties(uncertainies);
    m_tracks.emplace_back(std::move(track));
    return 1;
  }

  // Return with error
  return 0;
}

// ---------------------------------------------------------
int Particles::RemoveParticle(int index, Particles::ParticleType ptype) {
  if (ptype == ParticleType::kParton) {
    m_jets.erase(m_jets.begin() + index);
  } else if (ptype == ParticleType::kElectron) {
    m_electrons.erase(m_electrons.begin() + index);
  } else if (ptype == ParticleType::kMuon) {
    m_muons.erase(m_muons.begin() + index);
  } else if (ptype == ParticleType::kPhoton) {
    m_photons.erase(m_photons.begin() + index);
  } else if (ptype == ParticleType::kTau) {
    m_taus.erase(m_taus.begin() + index);
  } else if (ptype == ParticleType::kNeutrino) {
    m_neutrinos.erase(m_neutrinos.begin() + index);
  } else if (ptype == ParticleType::kBoson) {
    m_bosons.erase(m_bosons.begin() + index);
  } else if (ptype == ParticleType::kTrack) {
    m_tracks.erase(m_tracks.begin() + index);
  }

  // no error
  std::cout << "KLFitter::Particles::RemoveParticle(). Particle type " << ptype << " does not exist." << std::endl;
  return 1;
}

// ---------------------------------------------------------
int Particles::RemoveParticle(const std::string& name) {
  // get index and type
  TLorentzVector* vect = 0;
  int index = 0;
  Particles::ParticleType ptype = kParton;

  // remove particle
  if (FindParticle(name, vect, &index, &ptype)) {
    return RemoveParticle(index, ptype);
  } else {
    std::cout << "KLFitter::Particles::RemoveParticles(). Could not find particle with name " << name << "." << std::endl;
    return 0;
  }
}

// ---------------------------------------------------------
TLorentzVector* Particles::Particle(const std::string& name) {
  TLorentzVector* particle = 0;
  int index = 0;
  Particles::ParticleType ptype = kParton;

  // find particle
  if (!FindParticle(name, particle, &index, &ptype)) {
    std::cout << "KLFitter::Particles::Particle(). Could not find particles." << std::endl;
    return 0;
  }

  // return 4-vector
  return particle;
}

// ---------------------------------------------------------
TLorentzVector* Particles::Particle(int index, Particles::ParticleType ptype) {
  if (ptype == ParticleType::kParton) {
    return &m_jets.at(index).GetP4();
  } else if (ptype == ParticleType::kElectron) {
    return &m_electrons.at(index).GetP4();
  } else if (ptype == ParticleType::kMuon) {
    return &m_muons.at(index).GetP4();
  } else if (ptype == ParticleType::kPhoton) {
    return &m_photons.at(index).GetP4();
  } else if (ptype == ParticleType::kTau) {
    return &m_taus.at(index).GetP4();
  } else if (ptype == ParticleType::kNeutrino) {
    return &m_neutrinos.at(index).GetP4();
  } else if (ptype == ParticleType::kBoson) {
    return &m_bosons.at(index).GetP4();
  } else if (ptype == ParticleType::kTrack) {
    return &m_tracks.at(index).GetP4();
  }

  // Return nullptr
  std::cout << "KLFitter::Particles::Particle(). Particle type " << ptype << " does not exist." << std::endl;
  return nullptr;
}

// ---------------------------------------------------------
int Particles::FindParticle(const std::string& name, TLorentzVector* &particle, int *index, Particles::ParticleType *ptype) {
  // loop over all jets
  for (auto jet = m_jets.begin(); jet != m_jets.end(); ++jet) {
    if (name != jet->GetName()) continue;
    particle = &jet->GetP4();
    *index = jet - m_jets.begin();
    *ptype = Particles::kParton;
    return 1;
  }

  // loop over all electrons
  for (auto el = m_electrons.begin(); el != m_electrons.end(); ++el) {
    if (name != el->GetName()) continue;
    particle = &el->GetP4();
    *index = el - m_electrons.begin();
    *ptype = Particles::kElectron;
    return 1;
  }

  // loop over all muons
  for (auto mu = m_muons.begin(); mu != m_muons.end(); ++mu) {
    if (name != mu->GetName()) continue;
    particle = &mu->GetP4();
    *index = mu - m_muons.begin();
    *ptype = Particles::kMuon;
    return 1;
  }

  // loop over all taus
  for (auto tau = m_taus.begin(); tau != m_taus.end(); ++tau) {
    if (name != tau->GetName()) continue;
    particle = &tau->GetP4();
    *index = tau - m_taus.begin();
    *ptype = Particles::kTau;
    return 1;
  }

  // loop over all neutrinos
  for (auto neutrino = m_neutrinos.begin(); neutrino != m_neutrinos.end(); ++neutrino) {
    if (name != neutrino->GetName()) continue;
    particle = &neutrino->GetP4();
    *index = neutrino - m_neutrinos.begin();
    *ptype = Particles::kNeutrino;
    return 1;
  }

  // loop over all bosons
  for (auto boson = m_bosons.begin(); boson != m_bosons.end(); ++boson) {
    if (name != boson->GetName()) continue;
    particle = &boson->GetP4();
    *index = boson - m_bosons.begin();
    *ptype = Particles::kBoson;
    return 1;
  }

  // loop over all photons
  for (auto ph = m_photons.begin(); ph != m_photons.end(); ++ph) {
    if (name != ph->GetName()) continue;
    particle = &ph->GetP4();
    *index = ph - m_photons.begin();
    *ptype = Particles::kPhoton;
    return 1;
  }

  // loop over all tracks
  for (auto track = m_tracks.begin(); track != m_tracks.end(); ++track) {
    if (name != track->GetName()) continue;
    particle = &track->GetP4();
    *index = track - m_tracks.begin();
    *ptype = Particles::kTrack;
    return 1;
  }

  // particle not found
  return 0;
}

// ---------------------------------------------------------
TLorentzVector* Particles::Parton(int index) {
  return &m_jets.at(index).GetP4();
}

// ---------------------------------------------------------
TLorentzVector* Particles::Electron(int index) {
  return &m_electrons.at(index).GetP4();
}

// ---------------------------------------------------------
TLorentzVector* Particles::Muon(int index) {
  return &m_muons.at(index).GetP4();
}

// ---------------------------------------------------------
TLorentzVector* Particles::Tau(int index) {
  return &m_taus.at(index).GetP4();
}

// ---------------------------------------------------------
TLorentzVector* Particles::Boson(int index) {
  return &m_bosons.at(index).GetP4();
}

// ---------------------------------------------------------
TLorentzVector* Particles::Neutrino(int index) {
  return &m_neutrinos.at(index).GetP4();
}

// ---------------------------------------------------------
TLorentzVector* Particles::Photon(int index) {
  return &m_photons.at(index).GetP4();
}

// ---------------------------------------------------------
TLorentzVector* Particles::Track(int index) {
  return &m_tracks.at(index).GetP4();
}

// ---------------------------------------------------------
int Particles::NParticles(KLFitter::Particles::ParticleType ptype) const {
  if (ptype == ParticleType::kParton) {
    return static_cast<int>(m_jets.size());
  } else if (ptype == ParticleType::kElectron) {
    return static_cast<int>(m_electrons.size());
  } else if (ptype == ParticleType::kMuon) {
    return static_cast<int>(m_muons.size());
  } else if (ptype == ParticleType::kPhoton) {
    return static_cast<int>(m_photons.size());
  } else if (ptype == ParticleType::kTau) {
    return static_cast<int>(m_taus.size());
  } else if (ptype == ParticleType::kNeutrino) {
    return static_cast<int>(m_neutrinos.size());
  } else if (ptype == ParticleType::kBoson) {
    return static_cast<int>(m_bosons.size());
  } else if (ptype == ParticleType::kTrack) {
    return static_cast<int>(m_tracks.size());
  }
  return 0;
}

// ---------------------------------------------------------
std::string Particles::NameParticle(int index, Particles::ParticleType ptype) const {
  if (ptype == ParticleType::kParton) {
    return m_jets.at(index).GetName();
  } else if (ptype == ParticleType::kElectron) {
    return m_electrons.at(index).GetName();
  } else if (ptype == ParticleType::kMuon) {
    return m_muons.at(index).GetName();
  } else if (ptype == ParticleType::kPhoton) {
    return m_photons.at(index).GetName();
  } else if (ptype == ParticleType::kTau) {
    return m_taus.at(index).GetName();
  } else if (ptype == ParticleType::kNeutrino) {
    return m_neutrinos.at(index).GetName();
  } else if (ptype == ParticleType::kBoson) {
    return m_bosons.at(index).GetName();
  } else if (ptype == ParticleType::kTrack) {
    return m_tracks.at(index).GetName();
  }

  // return name
  std::cout << "KLFitter::Particles::NameParticle(). Particle type " << ptype << " does not exist." << std::endl;
  return "";
}

// ---------------------------------------------------------
double Particles::DetEta(int index, Particles::ParticleType ptype) const {
  if (ptype == Particles::kParton) {
    return m_jets.at(index).GetDetEta();
  } else if (ptype == Particles::kElectron) {
    return m_electrons.at(index).GetDetEta();
  } else if (ptype == Particles::kMuon) {
    return m_muons.at(index).GetDetEta();
  } else if (ptype == Particles::kPhoton) {
    return m_photons.at(index).GetDetEta();
  }

  // return error value
  std::cout << "KLFitter::Particles::DetEta(). Particle type " << ptype << " does not store detector eta." << std::endl;
  return -100;
}

// ---------------------------------------------------------
float Particles::LeptonCharge(int index, Particles::ParticleType ptype) const {
  if (ptype == Particles::kElectron) {
    return m_electrons.at(index).GetCharge();
  } else if (ptype == Particles::kMuon) {
    return m_muons.at(index).GetCharge();
  }

  // return error value
  std::cout << "KLFitter::Particles::LepCharge NO LEPTON TYPE!" << std::endl;
  return -9;
}

// ---------------------------------------------------------
const std::vector<double>* Particles::Uncertainties(int index, Particles::ParticleType ptype) const {
  if (ptype == Particles::kTrack) return &m_tracks.at(index).GetUncertainties();

  // return error value
  std::cout << "KLFitter::Particles::Uncertainties(). Particle type " << ptype << " does not store uncertainties." << std::endl;
  return nullptr;
}

// ---------------------------------------------------------
int Particles::JetIndex(int index) const {
  return m_jets.at(index).GetIdentifier();
}

// ---------------------------------------------------------
int Particles::ElectronIndex(int index) const {
  return m_electrons.at(index).GetIdentifier();
}

// ---------------------------------------------------------
int Particles::MuonIndex(int index) const {
 return m_muons.at(index).GetIdentifier();
}

// ---------------------------------------------------------
int Particles::PhotonIndex(int index) const {
  return m_photons.at(index).GetIdentifier();
}

// ---------------------------------------------------------
int Particles::TrackIndex(int index) const {
  return m_tracks.at(index).GetIdentifier();
}

// ---------------------------------------------------------
int Particles::SetIsBTagged(int index, bool isBTagged) {
  m_jets.at(index).SetIsBTagged(isBTagged);
  return 1;
}

// ---------------------------------------------------------
int Particles::SetBTagWeight(int index, double btagweight) {
  m_jets.at(index).SetBTagWeight(btagweight);
  return 1;
}

// ---------------------------------------------------------
int Particles::SetBTagWeightSet(int index, bool btagweightset) {
  m_jets.at(index).SetBTagWeightIsSet(btagweightset);
  return 1;
}

// ---------------------------------------------------------
int Particles::SetBTaggingEfficiency(int index, double btagEff) {
  m_jets.at(index).SetBTagEfficiency(btagEff);
  return 1;
}

// ---------------------------------------------------------
int Particles::SetBTaggingRejection(int index, double btagRej) {
  m_jets.at(index).SetBTagRejection(btagRej);
  return 1;
}

// ---------------------------------------------------------
int Particles::NBTags() const {
  int sum{0};
  for (const auto& jet : m_jets) {
    if (jet.GetIsBTagged()) sum++;
  }

  return sum;
}
}  // namespace KLFitter
