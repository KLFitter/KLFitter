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

// ---------------------------------------------------------
KLFitter::Particles::Particles() = default;

// ---------------------------------------------------------
KLFitter::Particles::Particles(const KLFitter::Particles& o) :
    m_parton_names(std::vector<std::string>{o.m_parton_names}),
    m_electron_names(std::vector<std::string>{o.m_electron_names}),
    m_muon_names(std::vector<std::string>{o.m_muon_names}),
    m_tau_names(std::vector<std::string>{o.m_tau_names}),
    m_neutrino_names(std::vector<std::string>{o.m_neutrino_names}),
    m_boson_names(std::vector<std::string>{o.m_boson_names}),
    m_photon_names(std::vector<std::string>{o.m_photon_names}),
    fJetIndex(std::vector<int>{o.fJetIndex}),
    fElectronIndex(std::vector<int>{o.fElectronIndex}),
    fMuonIndex(std::vector<int>{o.fMuonIndex}),
    fPhotonIndex(std::vector<int>{o.fPhotonIndex}),
    fTrueFlavor(std::vector<TrueFlavorType>{o.fTrueFlavor}),
    fIsBTagged(std::vector<bool>{o.fIsBTagged}),
    fBTaggingEfficiency(std::vector<double>{o.fBTaggingEfficiency}),
    fBTaggingRejection(std::vector<double>{o.fBTaggingRejection}),
    fBTagWeight(std::vector<double>{o.fBTagWeight}),
    fBTagWeightSet(std::vector<bool>{o.fBTagWeightSet}),
    fElectronDetEta(std::vector<double>{o.fElectronDetEta}),
    fMuonDetEta(std::vector<double>{o.fMuonDetEta}),
    fJetDetEta(std::vector<double>{o.fJetDetEta}),
    fPhotonDetEta(std::vector<double>{o.fPhotonDetEta}),
    fElectronCharge(std::vector<float>{o.fElectronCharge}),
    fMuonCharge(std::vector<float>{o.fMuonCharge}) {

  // Make deep copies of the vectors of unique pointers.
  m_partons.reserve(o.m_partons.size());
  for (const auto& i : o.m_partons) {
    m_partons.emplace_back(new TLorentzVector{*i});
  }

  m_electrons.reserve(o.m_electrons.size());
  for (const auto& i : o.m_electrons) {
    m_electrons.emplace_back(new TLorentzVector{*i});
  }

  m_muons.reserve(o.m_muons.size());
  for (const auto& i : o.m_muons) {
    m_muons.emplace_back(new TLorentzVector{*i});
  }

  m_taus.reserve(o.m_taus.size());
  for (const auto& i : o.m_taus) {
    m_taus.emplace_back(new TLorentzVector{*i});
  }

  m_neutrinos.reserve(o.m_neutrinos.size());
  for (const auto& i : o.m_neutrinos) {
    m_neutrinos.emplace_back(new TLorentzVector{*i});
  }

  m_bosons.reserve(o.m_bosons.size());
  for (const auto& i : o.m_bosons) {
    m_bosons.emplace_back(new TLorentzVector{*i});
  }

  m_photons.reserve(o.m_photons.size());
  for (const auto& i : o.m_photons) {
    m_photons.emplace_back(new TLorentzVector{*i});
  }
}

// ---------------------------------------------------------
KLFitter::Particles::~Particles() = default;

// ---------------------------------------------------------
KLFitter::Particles& KLFitter::Particles::operator=(const KLFitter::Particles& o) {
  m_parton_names = o.m_parton_names;
  m_electron_names = o.m_electron_names;
  m_muon_names = o.m_muon_names;
  m_tau_names = o.m_tau_names;
  m_neutrino_names = o.m_neutrino_names;
  m_boson_names = o.m_boson_names;
  m_photon_names = o.m_photon_names;

  fJetIndex = o.fJetIndex;
  fElectronIndex = o.fElectronIndex;
  fMuonIndex = o.fMuonIndex;
  fPhotonIndex = o.fPhotonIndex;
  fTrueFlavor = o.fTrueFlavor;
  fIsBTagged = o.fIsBTagged;
  fBTaggingEfficiency = o.fBTaggingEfficiency;
  fBTaggingRejection = o.fBTaggingRejection;
  fBTagWeight = o.fBTagWeight;
  fBTagWeightSet = o.fBTagWeightSet;
  fElectronDetEta = o.fElectronDetEta;
  fMuonDetEta = o.fMuonDetEta;
  fJetDetEta = o.fJetDetEta;
  fPhotonDetEta = o.fPhotonDetEta;
  fElectronCharge = o.fElectronCharge;
  fMuonCharge = o.fMuonCharge;

  // Make deep copies of the vectors of unique pointers.
  m_partons = std::vector<std::unique_ptr<TLorentzVector> >{};
  m_partons.reserve(o.m_partons.size());
  for (const auto& i : o.m_partons) {
    m_partons.emplace_back(new TLorentzVector{*i});
  }

  m_electrons = std::vector<std::unique_ptr<TLorentzVector> >{};
  m_electrons.reserve(o.m_electrons.size());
  for (const auto& i : o.m_electrons) {
    m_electrons.emplace_back(new TLorentzVector{*i});
  }

  m_muons = std::vector<std::unique_ptr<TLorentzVector> >{};
  m_muons.reserve(o.m_muons.size());
  for (const auto& i : o.m_muons) {
    m_muons.emplace_back(new TLorentzVector{*i});
  }

  m_taus = std::vector<std::unique_ptr<TLorentzVector> >{};
  m_taus.reserve(o.m_taus.size());
  for (const auto& i : o.m_taus) {
    m_taus.emplace_back(new TLorentzVector{*i});
  }

  m_neutrinos = std::vector<std::unique_ptr<TLorentzVector> >{};
  m_neutrinos.reserve(o.m_neutrinos.size());
  for (const auto& i : o.m_neutrinos) {
    m_neutrinos.emplace_back(new TLorentzVector{*i});
  }

  m_bosons = std::vector<std::unique_ptr<TLorentzVector> >{};
  m_bosons.reserve(o.m_bosons.size());
  for (const auto& i : o.m_bosons) {
    m_bosons.emplace_back(new TLorentzVector{*i});
  }

  m_photons = std::vector<std::unique_ptr<TLorentzVector> >{};
  m_photons.reserve(o.m_photons.size());
  for (const auto& i : o.m_photons) {
    m_photons.emplace_back(new TLorentzVector{*i});
  }

  return *this;
}

// ---------------------------------------------------------
int KLFitter::Particles::AddParticle(const TLorentzVector& particle, double DetEta, float LepCharge, KLFitter::Particles::ParticleType ptype, std::string name, int measuredindex) {
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
  KLFitter::Particles::ParticleType temptype = kParton;

  // check if particle with name exists already
  if (!FindParticle(name, vect, &index, &temptype)) {
    // add particle
    // create pointer copy of particle content which is owend by Particles
    std::unique_ptr<TLorentzVector> cparticle{new TLorentzVector{particle.Px(), particle.Py(), particle.Pz(), particle.E()}};
    container->emplace_back(std::move(cparticle));
    ParticleNameContainer(ptype)->push_back(name);
    if (ptype == KLFitter::Particles::kElectron) {
      fElectronIndex.push_back(measuredindex);
      fElectronDetEta.push_back(DetEta);
      fElectronCharge.push_back(LepCharge);
    } else if (ptype == KLFitter::Particles::kMuon) {
      fMuonIndex.push_back(measuredindex);
      fMuonDetEta.push_back(DetEta);
      fMuonCharge.push_back(LepCharge);
    } else if (ptype == KLFitter::Particles::kPhoton) {
      fPhotonIndex.push_back(measuredindex);
      fPhotonDetEta.push_back(DetEta);
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
int KLFitter::Particles::AddParticle(const TLorentzVector* const particle, double DetEta, float LepCharge, KLFitter::Particles::ParticleType ptype, std::string name, int measuredindex) {
  return AddParticle(*particle, DetEta, LepCharge, ptype, name, measuredindex);
}

// ---------------------------------------------------------
int KLFitter::Particles::AddParticle(const TLorentzVector& particle, double DetEta, KLFitter::Particles::ParticleType ptype, std::string name, int measuredindex, bool isBtagged, double bTagEff, double bTagRej, TrueFlavorType trueflav, double btagweight) {
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
  KLFitter::Particles::ParticleType temptype = kParton;

  // check if particle with name exists already
  if (!FindParticle(name, vect, &index, &temptype)) {
    // add particle
    // create pointer copy of particle content which is owend by Particles
    std::unique_ptr<TLorentzVector> cparticle{new TLorentzVector{particle.Px(), particle.Py(), particle.Pz(), particle.E()}};
    container->emplace_back(std::move(cparticle));
    ParticleNameContainer(ptype)->push_back(name);
    if (ptype == KLFitter::Particles::kParton) {
      fTrueFlavor.push_back(trueflav);
      fIsBTagged.push_back(isBtagged);
      fBTaggingEfficiency.push_back(bTagEff);
      fBTaggingRejection.push_back(bTagRej);
      fJetIndex.push_back(measuredindex);
      fJetDetEta.push_back(DetEta);
      fBTagWeight.push_back(btagweight);
      if (btagweight != 999) {
        fBTagWeightSet.push_back(true);
      } else {
        fBTagWeightSet.push_back(false);
      }
    } else if (ptype == KLFitter::Particles::kElectron) {
      fElectronIndex.push_back(measuredindex);
      fElectronDetEta.push_back(DetEta);
    } else if (ptype == KLFitter::Particles::kMuon) {
      fMuonIndex.push_back(measuredindex);
      fMuonDetEta.push_back(DetEta);
    } else if (ptype == KLFitter::Particles::kPhoton) {
      fPhotonIndex.push_back(measuredindex);
      fPhotonDetEta.push_back(DetEta);
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
int KLFitter::Particles::AddParticle(const TLorentzVector* const particle, double DetEta, KLFitter::Particles::ParticleType ptype, std::string name, int measuredindex, bool isBtagged, double bTagEff, double bTagRej, TrueFlavorType trueflav, double btagweight) {
  return AddParticle(*particle, DetEta, ptype, name, measuredindex, isBtagged, bTagEff, bTagRej, trueflav, btagweight);
}

// ---------------------------------------------------------
int KLFitter::Particles::AddParticle(const TLorentzVector& particle, KLFitter::Particles::ParticleType ptype, std::string name, int measuredindex, bool isBtagged, double bTagEff, double bTagRej, TrueFlavorType trueflav, double btagweight) {
  return AddParticle(particle, -999, ptype, name, measuredindex, isBtagged, bTagEff, bTagRej, trueflav, btagweight);
}

// --------------------------------------------------------- // THIS FUNCTION IS TO BE REMOVED IN THE NEXT MAJOR RELEASE
int KLFitter::Particles::AddParticle(const TLorentzVector* const particle, KLFitter::Particles::ParticleType ptype, std::string name, int measuredindex, bool isBtagged, double bTagEff, double bTagRej, TrueFlavorType trueflav, double btagweight) {
  return AddParticle(*particle, ptype, name, measuredindex, isBtagged, bTagEff, bTagRej, trueflav, btagweight);
}

// ---------------------------------------------------------
int KLFitter::Particles::AddParticle(const TLorentzVector& particle, KLFitter::Particles::ParticleType ptype, std::string name, int measuredindex, TrueFlavorType trueflav, double btagweight) {
  return AddParticle(particle, -999, ptype, name, measuredindex, false, -1., -1., trueflav, btagweight);
}

// --------------------------------------------------------- // THIS FUNCTION IS TO BE REMOVED IN THE NEXT MAJOR RELEASE
int KLFitter::Particles::AddParticle(const TLorentzVector* const particle, KLFitter::Particles::ParticleType ptype, std::string name, int measuredindex, TrueFlavorType trueflav, double btagweight) {
  return AddParticle(*particle, ptype, name, measuredindex, trueflav, btagweight);
}

// ---------------------------------------------------------
int KLFitter::Particles::RemoveParticle(int index, KLFitter::Particles::ParticleType ptype) {
  // check container and index
  if (!CheckIndex(ParticleContainer(ptype), index))
    return 0;

  // remove particle
  ParticleContainer(ptype)->erase(ParticleContainer(ptype)->begin() + index);
  ParticleNameContainer(ptype)->erase(ParticleNameContainer(ptype)->begin() + index);

  // no error
  return 1;
}

// ---------------------------------------------------------
int KLFitter::Particles::RemoveParticle(std::string name) {
  // get index and type
  TLorentzVector* vect = 0;
  int index = 0;
  KLFitter::Particles::ParticleType ptype = kParton;

  // remove particle
  if (FindParticle(name, vect, &index, &ptype)) {
    return RemoveParticle(index, ptype);
  } else {
    std::cout << "KLFitter::Particles::RemoveParticles(). Could not find particle with name " << name << "." << std::endl;
    return 0;
  }
}

// ---------------------------------------------------------
TLorentzVector* KLFitter::Particles::Particle(std::string name) {
  TLorentzVector* particle = 0;
  int index = 0;
  KLFitter::Particles::ParticleType ptype = kParton;

  // find particle
  if (!FindParticle(name, particle, &index, &ptype)) {
    std::cout << "KLFitter::Particles::Particle(). Could not find particles." << std::endl;
    return 0;
  }

  // return 4-vector
  return particle;
}

// ---------------------------------------------------------
TLorentzVector* KLFitter::Particles::Particle(int index, KLFitter::Particles::ParticleType ptype) {
  // get particle container
  auto container = ParticleContainer(ptype);

  if (index < 0 || index > NParticles(ptype)) {
    std::cout << "KLFitter::Particles::Particle(). Index out of range." << std::endl;
    return 0;
  }

  // return pointer
  return (*container)[index].get();
}

// ---------------------------------------------------------
int KLFitter::Particles::FindParticle(std::string name, TLorentzVector* &particle, int *index, KLFitter::Particles::ParticleType *ptype) {
  // loop over all partons
  unsigned int npartons = m_parton_names.size();
  for (unsigned int i = 0; i < npartons; ++i)
    if (name == m_parton_names[i]) {
      particle = m_partons[i].get();
      *index = i;
      *ptype = KLFitter::Particles::kParton;
      return 1;
    }

  // loop over all electrons
  unsigned int nelectrons = m_electron_names.size();
  for (unsigned int i = 0; i < nelectrons; ++i)
    if (name == m_electron_names[i]) {
      particle = m_electrons[i].get();
      *index = i;
      *ptype = KLFitter::Particles::kElectron;
      return 1;
    }

  // loop over all muons
  unsigned int nmuons = m_muon_names.size();
  for (unsigned int i = 0; i < nmuons; ++i)
    if (name == m_muon_names[i]) {
      particle = m_muons[i].get();
      *index = i;
      *ptype = KLFitter::Particles::kMuon;
      return 1;
    }

  // loop over all taus
  unsigned int ntaus = m_tau_names.size();
  for (unsigned int i = 0; i < ntaus; ++i)
    if (name == m_tau_names[i]) {
      particle = m_taus[i].get();
      *index = i;
      *ptype = KLFitter::Particles::kTau;
      return 1;
    }

  // loop over all neutrinos
  unsigned int nneutrinos = m_neutrino_names.size();
  for (unsigned int i = 0; i < nneutrinos; ++i)
    if (name == m_neutrino_names[i]) {
      particle = m_neutrinos[i].get();
      *index = i;
      *ptype = KLFitter::Particles::kNeutrino;
      return 1;
    }

  // loop over all bosons
  unsigned int nbosons = m_boson_names.size();
  for (unsigned int i = 0; i < nbosons; ++i)
    if (name == m_boson_names[i]) {
      particle = m_bosons[i].get();
      *index = i;
      *ptype = KLFitter::Particles::kBoson;
      return 1;
    }

  // loop over all photons
  unsigned int nphotons = m_photon_names.size();
  for (unsigned int i = 0; i < nphotons; ++i)
    if (name == m_photon_names[i]) {
      particle = m_photons[i].get();
      *index = i;
      *ptype = KLFitter::Particles::kPhoton;
      return 1;
    }

  // particle not found
  return 0;
}

// ---------------------------------------------------------
TLorentzVector* KLFitter::Particles::Parton(int index) {
  // no check on index range for CPU-time reasons
  return m_partons[index].get();
}

// ---------------------------------------------------------
TLorentzVector* KLFitter::Particles::Electron(int index) {
  // no check on index range for CPU-time reasons
  return m_electrons[index].get();
}

// ---------------------------------------------------------
TLorentzVector* KLFitter::Particles::Muon(int index) {
  // no check on index range for CPU-time reasons
  return m_muons[index].get();
}

// ---------------------------------------------------------
TLorentzVector* KLFitter::Particles::Tau(int index) {
  // no check on index range for CPU-time reasons
  return m_taus[index].get();
}

// ---------------------------------------------------------
TLorentzVector* KLFitter::Particles::Boson(int index) {
  // no check on index range for CPU-time reasons
  return m_bosons[index].get();
}

// ---------------------------------------------------------
TLorentzVector* KLFitter::Particles::Neutrino(int index) {
  // no check on index range for CPU-time reasons
  return m_neutrinos[index].get();
}

// ---------------------------------------------------------
TLorentzVector* KLFitter::Particles::Photon(int index) {
  // no check on index range for CPU-time reasons
  return m_photons[index].get();
}

// ---------------------------------------------------------
int KLFitter::Particles::NParticles(KLFitter::Particles::ParticleType ptype) {
  return static_cast<int>(ParticleContainer(ptype)->size());
}

// ---------------------------------------------------------
std::string KLFitter::Particles::NameParticle(int index, KLFitter::Particles::ParticleType ptype) {
  // get particle container
  auto container = ParticleContainer(ptype);

  // check container and index
  if (!CheckIndex(container, index))
    return "";

  // return name
  return (*ParticleNameContainer(ptype))[index];
}

// ---------------------------------------------------------
int KLFitter::Particles::CheckIndex(std::vector<std::unique_ptr<TLorentzVector> >* container, int index) {
  // check container
  if (!container) {
    std::cout << "KLFitter::Particles::CheckIndex(). Container does not exist." << std::endl;
    return 0;
  }

  // check index
  if (index < 0 || index >= static_cast<int>(container->size())) {
    std::cout << "KLFitter::Particles::CheckIndex(). Index out of range." << std::endl;
    return 0;
  }

  // no error
  return 1;
}

// ---------------------------------------------------------
std::vector<std::unique_ptr<TLorentzVector> >* KLFitter::Particles::ParticleContainer(KLFitter::Particles::ParticleType ptype) {
  // return particle container
  switch (ptype) {
  case KLFitter::Particles::kParton:
    return &m_partons;
    break;
  case KLFitter::Particles::kElectron:
    return &m_electrons;
    break;
  case KLFitter::Particles::kMuon:
    return &m_muons;
    break;
  case KLFitter::Particles::kTau:
    return &m_taus;
    break;
  case KLFitter::Particles::kNeutrino:
    return &m_neutrinos;
    break;
  case KLFitter::Particles::kBoson:
    return &m_bosons;
    break;
  case KLFitter::Particles::kPhoton:
    return &m_photons;
    break;
  }

  // or null pointer
  std::cout << "KLFitter::Particles::ParticleContainer(). Particle type unknown." << std::endl;
  return 0;
}

// ---------------------------------------------------------
std::vector <std::string>* KLFitter::Particles::ParticleNameContainer(KLFitter::Particles::ParticleType ptype) {
  // return container
  if (ptype == KLFitter::Particles::kParton) {
    return &m_parton_names;
  } else if (ptype == KLFitter::Particles::kElectron) {
    return &m_electron_names;
  } else if (ptype == KLFitter::Particles::kMuon) {
    return &m_muon_names;
  } else if (ptype == KLFitter::Particles::kTau) {
    return &m_tau_names;
  } else if (ptype == KLFitter::Particles::kBoson) {
    return &m_boson_names;
  } else if (ptype == KLFitter::Particles::kNeutrino) {
    return &m_neutrino_names;
  } else if (ptype == KLFitter::Particles::kPhoton) {
    return &m_photon_names;
  } else {
    // or null pointer
    std::cout << "KLFitter::Particles::ParticleNameContainer(). Particle type not known." << std::endl;
    return 0;
  }
}

// ---------------------------------------------------------
double KLFitter::Particles::DetEta(int index, KLFitter::Particles::ParticleType ptype) {
  if (index < 0 || index > NParticles(ptype)) {
    std::cout << "KLFitter::Particles::DetEta(). Index out of range." << std::endl;
    return 0;
  }

  if (ptype == KLFitter::Particles::kParton) {
    return fJetDetEta[index];
  } else if (ptype == KLFitter::Particles::kElectron) {
    return fElectronDetEta[index];
  } else if (ptype == KLFitter::Particles::kMuon) {
    return fMuonDetEta[index];
  } else if (ptype == KLFitter::Particles::kPhoton) {
    return fPhotonDetEta[index];
  }

  // return error value
  return -100;
}

// ---------------------------------------------------------
float KLFitter::Particles::LeptonCharge(int index, KLFitter::Particles::ParticleType ptype) {
  if (index < 0 || index > NParticles(ptype)) {
    std::cout << "KLFitter::Particles::LepCharge(). Index out of range." << std::endl;
    return 0;
  }

  if (ptype == KLFitter::Particles::kElectron) {
    if (fElectronCharge.size()== 0) {
      return -9;
    } else {
      return fElectronCharge[index];
    }
  } else if (ptype == KLFitter::Particles::kMuon) {
    if (fMuonCharge.size()== 0) {
      return -9;
    } else {
      return fMuonCharge[index];
    }
  } else {
    std::cout << "KLFitter::Particles::LepCharge NO LEPTON TYPE!" << std::endl;
  }

  // return error value
  return -9;
}

// ---------------------------------------------------------
int KLFitter::Particles::JetIndex(int index) {
  // no check on index range for CPU-time reasons
  return fJetIndex[index];
}

// ---------------------------------------------------------
int KLFitter::Particles::ElectronIndex(int index) {
  // no check on index range for CPU-time reasons
  return fElectronIndex[index];
}

// ---------------------------------------------------------
int KLFitter::Particles::MuonIndex(int index) {
  // no check on index range for CPU-time reasons
  return fMuonIndex[index];
}

// ---------------------------------------------------------
int KLFitter::Particles::PhotonIndex(int index) {
  // no check on index range for CPU-time reasons
  return fPhotonIndex[index];
}

// ---------------------------------------------------------
int KLFitter::Particles::SetIsBTagged(int index, bool isBTagged) {
  // check index
  if (index < 0 || index >= static_cast<int>(fIsBTagged.size())) {
    std::cout << " KLFitter::SetIsBTagged(). Index out of range." << std::endl;
    return 0;
  }

  fIsBTagged[index] = isBTagged;

  return 1;
}

// ---------------------------------------------------------
int KLFitter::Particles::SetBTagWeight(int index, double btagweight) {
  // check index
  if (index < 0 || index >= static_cast<int>(fBTagWeight.size())) {
    std::cout << " KLFitter::SetBTagWeight(). Index out of range." << std::endl;
    return 0;
  }

  fBTagWeight[index] = btagweight;
  SetBTagWeightSet(index, true);

  return 1;
}

// ---------------------------------------------------------
int KLFitter::Particles::SetBTagWeightSet(int index, bool btagweightset) {
  // check index
  if (index < 0 || index >= static_cast<int>(fBTagWeightSet.size())) {
    std::cout << " KLFitter::SetBTagWeightSet(). Index out of range." << std::endl;
    return 0;
  }

  fBTagWeightSet[index] = btagweightset;

  return 1;
}

// ---------------------------------------------------------
int KLFitter::Particles::SetBTaggingEfficiency(int index, double btagEff) {
  // check index
  if (index < 0 || index >= static_cast<int>(fBTaggingEfficiency.size())) {
    std::cout << " KLFitter::SetBTaggingEfficiency(). Index out of range." << std::endl;
    return 0;
  }

  fBTaggingEfficiency[index] = btagEff;

  return 1;
}

// ---------------------------------------------------------
int KLFitter::Particles::SetBTaggingRejection(int index, double btagRej) {
  // check index
  if (index < 0 || index >= static_cast<int>(fBTaggingRejection.size())) {
    std::cout << " KLFitter::SetBTaggingRejection(). Index out of range." << std::endl;
    return 0;
  }

  fBTaggingRejection[index] = btagRej;

  return 1;
}

// ---------------------------------------------------------
int KLFitter::Particles::NBTags() {
  unsigned int n = fIsBTagged.size();
  int sum = 0;

  for (unsigned int i = 0; i < n; ++i) {
    if (fIsBTagged[i]) {
      sum++;
    }
  }

  return sum;
}
