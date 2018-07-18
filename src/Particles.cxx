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
  // get particle container
  auto container = ParticleContainer(ptype);

  // check if container exists
  if (!container) {
    std::cout << "KLFitter::Particles::AddParticle(). Container does not exist." << std::endl;
    return 0;
  }

  // check name
  if (name == "")
    name = Form("particle_%i", NParticles());

  // get index and type
  TLorentzVector* vect = 0;
  int index = 0;
  Particles::ParticleType temptype = kParton;

  // check if particle with name exists already
  if (!FindParticle(name, vect, &index, &temptype)) {
    container->emplace_back(particle);
    ParticleNameContainer(ptype)->push_back(name);
    if (ptype == Particles::kElectron) {
      m_electron_indices.push_back(measuredindex);
      m_electron_det_etas.push_back(DetEta);
      m_electron_charges.push_back(LepCharge);
    } else if (ptype == Particles::kMuon) {
      m_muon_indices.push_back(measuredindex);
      m_muon_det_etas.push_back(DetEta);
      m_muon_charges.push_back(LepCharge);
    } else if (ptype == Particles::kPhoton) {
      m_photon_indices.push_back(measuredindex);
      m_photon_det_etas.push_back(DetEta);
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
  // get particle container
  auto container = ParticleContainer(ptype);

  // check if container exists
  if (!container) {
    std::cout << "KLFitter::Particles::AddParticle(). Container does not exist." << std::endl;
    return 0;
  }

  // check name
  if (name == "")
    name = Form("particle_%i", NParticles());

  // get index and type
  TLorentzVector* vect = 0;
  int index = 0;
  Particles::ParticleType temptype = kParton;

  // check if particle with name exists already
  if (!FindParticle(name, vect, &index, &temptype)) {
    container->emplace_back(particle);
    ParticleNameContainer(ptype)->push_back(name);
    if (ptype == Particles::kParton) {
      m_true_flavors.push_back(trueflav);
      m_jet_btagged_bools.push_back(isBtagged);
      m_btag_efficiencies.push_back(bTagEff);
      m_btag_rejections.push_back(bTagRej);
      m_jet_indices.push_back(measuredindex);
      m_jet_det_etas.push_back(DetEta);
      m_btag_weights.push_back(btagweight);
      if (btagweight != 999) {
        m_btag_weights_set.push_back(true);
      } else {
        m_btag_weights_set.push_back(false);
      }
    } else if (ptype == Particles::kElectron) {
      m_electron_indices.push_back(measuredindex);
      m_electron_det_etas.push_back(DetEta);
    } else if (ptype == Particles::kMuon) {
      m_muon_indices.push_back(measuredindex);
      m_muon_det_etas.push_back(DetEta);
    } else if (ptype == Particles::kPhoton) {
      m_photon_indices.push_back(measuredindex);
      m_photon_det_etas.push_back(DetEta);
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
int Particles::RemoveParticle(int index, Particles::ParticleType ptype) {
  // check container and index
  if (!CheckIndex(*ParticleContainer(ptype), index))
    return 0;

  // remove particle
  ParticleContainer(ptype)->erase(ParticleContainer(ptype)->begin() + index);
  ParticleNameContainer(ptype)->erase(ParticleNameContainer(ptype)->begin() + index);

  // no error
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
  // get particle container
  auto container = ParticleContainer(ptype);

  if (index < 0 || index > NParticles(ptype)) {
    std::cout << "KLFitter::Particles::Particle(). Index out of range." << std::endl;
    return 0;
  }

  // return pointer
  return &container->at(index);
}

// ---------------------------------------------------------
int Particles::FindParticle(const std::string& name, TLorentzVector* &particle, int *index, Particles::ParticleType *ptype) {
  // loop over all partons
  unsigned int npartons = m_parton_names.size();
  for (unsigned int i = 0; i < npartons; ++i)
    if (name == m_parton_names[i]) {
      particle = &m_partons[i];
      *index = i;
      *ptype = Particles::kParton;
      return 1;
    }

  // loop over all electrons
  unsigned int nelectrons = m_electron_names.size();
  for (unsigned int i = 0; i < nelectrons; ++i)
    if (name == m_electron_names[i]) {
      particle = &m_electrons[i];
      *index = i;
      *ptype = Particles::kElectron;
      return 1;
    }

  // loop over all muons
  unsigned int nmuons = m_muon_names.size();
  for (unsigned int i = 0; i < nmuons; ++i)
    if (name == m_muon_names[i]) {
      particle = &m_muons[i];
      *index = i;
      *ptype = Particles::kMuon;
      return 1;
    }

  // loop over all taus
  unsigned int ntaus = m_tau_names.size();
  for (unsigned int i = 0; i < ntaus; ++i)
    if (name == m_tau_names[i]) {
      particle = &m_taus[i];
      *index = i;
      *ptype = Particles::kTau;
      return 1;
    }

  // loop over all neutrinos
  unsigned int nneutrinos = m_neutrino_names.size();
  for (unsigned int i = 0; i < nneutrinos; ++i)
    if (name == m_neutrino_names[i]) {
      particle = &m_neutrinos[i];
      *index = i;
      *ptype = Particles::kNeutrino;
      return 1;
    }

  // loop over all bosons
  unsigned int nbosons = m_boson_names.size();
  for (unsigned int i = 0; i < nbosons; ++i)
    if (name == m_boson_names[i]) {
      particle = &m_bosons[i];
      *index = i;
      *ptype = Particles::kBoson;
      return 1;
    }

  // loop over all photons
  unsigned int nphotons = m_photon_names.size();
  for (unsigned int i = 0; i < nphotons; ++i)
    if (name == m_photon_names[i]) {
      particle = &m_photons[i];
      *index = i;
      *ptype = Particles::kPhoton;
      return 1;
    }

  // particle not found
  return 0;
}

// ---------------------------------------------------------
TLorentzVector* Particles::Parton(int index) {
  // no check on index range for CPU-time reasons
  return &m_partons[index];
}

// ---------------------------------------------------------
TLorentzVector* Particles::Electron(int index) {
  // no check on index range for CPU-time reasons
  return &m_electrons[index];
}

// ---------------------------------------------------------
TLorentzVector* Particles::Muon(int index) {
  // no check on index range for CPU-time reasons
  return &m_muons[index];
}

// ---------------------------------------------------------
TLorentzVector* Particles::Tau(int index) {
  // no check on index range for CPU-time reasons
  return &m_taus[index];
}

// ---------------------------------------------------------
TLorentzVector* Particles::Boson(int index) {
  // no check on index range for CPU-time reasons
  return &m_bosons[index];
}

// ---------------------------------------------------------
TLorentzVector* Particles::Neutrino(int index) {
  // no check on index range for CPU-time reasons
  return &m_neutrinos[index];
}

// ---------------------------------------------------------
TLorentzVector* Particles::Photon(int index) {
  // no check on index range for CPU-time reasons
  return &m_photons[index];
}

// ---------------------------------------------------------
int Particles::NParticles(KLFitter::Particles::ParticleType ptype) {
  return static_cast<int>(ParticleContainer(ptype)->size());
}

// ---------------------------------------------------------
std::string Particles::NameParticle(int index, Particles::ParticleType ptype) {
  // get particle container
  auto container = ParticleContainer(ptype);

  // check container and index
  if (!CheckIndex(*container, index))
    return "";

  // return name
  return ParticleNameContainer(ptype)->at(index);
}

// ---------------------------------------------------------
int Particles::CheckIndex(const std::vector<TLorentzVector>& container, int index) const {
  // check index
  if (index < 0 || index >= static_cast<int>(container.size())) {
    std::cout << "KLFitter::Particles::CheckIndex(). Index out of range." << std::endl;
    return 0;
  }

  // no error
  return 1;
}

// ---------------------------------------------------------
std::vector<TLorentzVector>* Particles::ParticleContainer(KLFitter::Particles::ParticleType ptype) {
  // return particle container
  switch (ptype) {
  case Particles::kParton:
    return &m_partons;
    break;
  case Particles::kElectron:
    return &m_electrons;
    break;
  case Particles::kMuon:
    return &m_muons;
    break;
  case Particles::kTau:
    return &m_taus;
    break;
  case Particles::kNeutrino:
    return &m_neutrinos;
    break;
  case Particles::kBoson:
    return &m_bosons;
    break;
  case Particles::kPhoton:
    return &m_photons;
    break;
  }

  // or null pointer
  std::cout << "KLFitter::Particles::ParticleContainer(). Particle type unknown." << std::endl;
  return 0;
}

// ---------------------------------------------------------
std::vector <std::string>* Particles::ParticleNameContainer(KLFitter::Particles::ParticleType ptype) {
  // return container
  if (ptype == Particles::kParton) {
    return &m_parton_names;
  } else if (ptype == Particles::kElectron) {
    return &m_electron_names;
  } else if (ptype == Particles::kMuon) {
    return &m_muon_names;
  } else if (ptype == Particles::kTau) {
    return &m_tau_names;
  } else if (ptype == Particles::kBoson) {
    return &m_boson_names;
  } else if (ptype == Particles::kNeutrino) {
    return &m_neutrino_names;
  } else if (ptype == Particles::kPhoton) {
    return &m_photon_names;
  } else {
    // or null pointer
    std::cout << "KLFitter::Particles::ParticleNameContainer(). Particle type not known." << std::endl;
    return 0;
  }
}

// ---------------------------------------------------------
double Particles::DetEta(int index, Particles::ParticleType ptype) {
  if (index < 0 || index > NParticles(ptype)) {
    std::cout << "KLFitter::Particles::DetEta(). Index out of range." << std::endl;
    return 0;
  }

  if (ptype == Particles::kParton) {
    return m_jet_det_etas[index];
  } else if (ptype == Particles::kElectron) {
    return m_electron_det_etas[index];
  } else if (ptype == Particles::kMuon) {
    return m_muon_det_etas[index];
  } else if (ptype == Particles::kPhoton) {
    return m_photon_det_etas[index];
  }

  // return error value
  return -100;
}

// ---------------------------------------------------------
float Particles::LeptonCharge(int index, Particles::ParticleType ptype) {
  if (index < 0 || index > NParticles(ptype)) {
    std::cout << "KLFitter::Particles::LepCharge(). Index out of range." << std::endl;
    return 0;
  }

  if (ptype == Particles::kElectron) {
    if (m_electron_charges.size()== 0) {
      return -9;
    } else {
      return m_electron_charges[index];
    }
  } else if (ptype == Particles::kMuon) {
    if (m_muon_charges.size()== 0) {
      return -9;
    } else {
      return m_muon_charges[index];
    }
  } else {
    std::cout << "KLFitter::Particles::LepCharge NO LEPTON TYPE!" << std::endl;
  }

  // return error value
  return -9;
}

// ---------------------------------------------------------
int Particles::JetIndex(int index) const {
  // no check on index range for CPU-time reasons
  return m_jet_indices[index];
}

// ---------------------------------------------------------
int Particles::ElectronIndex(int index) const {
  // no check on index range for CPU-time reasons
  return m_electron_indices[index];
}

// ---------------------------------------------------------
int Particles::MuonIndex(int index) const {
  // no check on index range for CPU-time reasons
  return m_muon_indices[index];
}

// ---------------------------------------------------------
int Particles::PhotonIndex(int index) const {
  // no check on index range for CPU-time reasons
  return m_photon_indices[index];
}

// ---------------------------------------------------------
int Particles::SetIsBTagged(int index, bool isBTagged) {
  // check index
  if (index < 0 || index >= static_cast<int>(m_jet_btagged_bools.size())) {
    std::cout << "KLFitter::SetIsBTagged(). Index out of range." << std::endl;
    return 0;
  }

  m_jet_btagged_bools[index] = isBTagged;

  return 1;
}

// ---------------------------------------------------------
int Particles::SetBTagWeight(int index, double btagweight) {
  // check index
  if (index < 0 || index >= static_cast<int>(m_btag_weights.size())) {
    std::cout << "KLFitter::SetBTagWeight(). Index out of range." << std::endl;
    return 0;
  }

  m_btag_weights[index] = btagweight;
  SetBTagWeightSet(index, true);

  return 1;
}

// ---------------------------------------------------------
int Particles::SetBTagWeightSet(int index, bool btagweightset) {
  // check index
  if (index < 0 || index >= static_cast<int>(m_btag_weights_set.size())) {
    std::cout << "KLFitter::SetBTagWeightSet(). Index out of range." << std::endl;
    return 0;
  }

  m_btag_weights_set[index] = btagweightset;

  return 1;
}

// ---------------------------------------------------------
int Particles::SetBTaggingEfficiency(int index, double btagEff) {
  // check index
  if (index < 0 || index >= static_cast<int>(m_btag_efficiencies.size())) {
    std::cout << "KLFitter::SetBTaggingEfficiency(). Index out of range." << std::endl;
    return 0;
  }

  m_btag_efficiencies[index] = btagEff;

  return 1;
}

// ---------------------------------------------------------
int Particles::SetBTaggingRejection(int index, double btagRej) {
  // check index
  if (index < 0 || index >= static_cast<int>(m_btag_rejections.size())) {
    std::cout << "KLFitter::SetBTaggingRejection(). Index out of range." << std::endl;
    return 0;
  }

  m_btag_rejections[index] = btagRej;

  return 1;
}

// ---------------------------------------------------------
int Particles::NBTags() const {
  unsigned int n = m_jet_btagged_bools.size();
  int sum = 0;

  for (unsigned int i = 0; i < n; ++i) {
    if (m_jet_btagged_bools[i]) {
      sum++;
    }
  }

  return sum;
}
}  // namespace KLFitter
