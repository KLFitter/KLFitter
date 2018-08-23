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

#ifndef KLFITTER_PARTICLES_H_
#define KLFITTER_PARTICLES_H_

#include <string>
#include <vector>

#include "TLorentzVector.h"
#include "KLFitter/Particle/Boson.h"
#include "KLFitter/Particle/Electron.h"
#include "KLFitter/Particle/Jet.h"
#include "KLFitter/Particle/Muon.h"
#include "KLFitter/Particle/Neutrino.h"
#include "KLFitter/Particle/Photon.h"
#include "KLFitter/Particle/Tau.h"
#include "KLFitter/Particle/Track.h"

// ---------------------------------------------------------

namespace KLFitter {
/**
 * A class to hold collections of various particle types. The
 * possible particle types are defined in
 * Particles::ParticleType, and new particles are added via one
 * of the overloaded AddParticle() functions. It provides
 * functionality to retrieve information about the contained
 * particles via various accessor methods.
 */
class Particles final {
 public:
  /// Enumerator for particle types.
  enum ParticleType { kParton,    ///< Particle type parton
                      kElectron,  ///< Particle type electron
                      kMuon,      ///< Particle type muon
                      kTau,       ///< Particle type tau
                      kNeutrino,  ///< Particle type neutrino
                      kBoson,     ///< Particle type boson
                      kPhoton,    ///< Particle type photon
                      kTrack      ///< Particle type track
  };

  /// An enumerator for the true jet flavor.
  enum TrueFlavorType { kLight,     ///< Light quark
                        kB,         ///< B-quark
                        kLightUp,   ///< Up-type light quark
                        kLightDown, ///< Down-type light quark
                        kNone       ///< Not specified
  };

  /// The default constructor.
  Particles();

  /// The (defaulted) copy constructor.
  explicit Particles(const Particles& o);

  /// The (defaulted) destructor.
  ~Particles();

  /// The (defaulted) assignment operator.
  Particles& operator=(const Particles& o);

  /** \name Member functions (Get)  */
  /** @{ */

  /**
   * Return the number of partons.
   * @return The number of partons.
   */
  int NPartons() const { return static_cast<int>(m_jets.size()); }

  /**
   * Return the number of electrons.
   * @return The number of electrons.
   */
  int NElectrons() const { return static_cast<int>(m_electrons.size()); }

  /**
   * Return the number of muons.
   * @return The number of muons.
   */
  int NMuons() const { return int (m_muons.size()); }

  /**
   * Return the number of taus.
   * @return The number of taus.
   */
  int NTaus() const { return int (m_taus.size()); }

  /**
   * Return the number of neutrinos.
   * @return The number of neutrinos.
   */
  int NNeutrinos() const { return int (m_neutrinos.size()); }

  /**
   * Return the number of bosons.
   * @return The number of bosons.
   */
  int NBosons() const { return int (m_bosons.size()); }

  /**
   * Return the number of photons.
   * @return The number of photons.
   */
  int NPhotons() const { return int (m_photons.size()); }

  /**
   * Return the number of tracks.
   * @return The number of tracks.
   */
  int NTracks() const { return int (m_tracks.size()); }

  /**
   * Return the particle with a certain name
   * @param name The name of the particle.
   * @return A pointer to the TLorentzVector of the particle.
   */
  TLorentzVector* Particle(const std::string& name);

  /**
   * Return a particle with some index and type.
   * @param index The index of the particle.
   * @param ptype The type of the particle.
   * @return A pointer to the TLorentzVector of the particle.
   */
  TLorentzVector* Particle(int index, KLFitter::Particles::ParticleType ptype);

  /**
   * Find a particle by name.
   * @param name The name of the particle.
   * @param particle A reference to the pointer to the 4-vector.
   * @param index The pointer to the index.
   * @param ptype The pointer to the particle type.
   * @return A flag (1: found, 0: not found).
   */
  int FindParticle(const std::string& name, TLorentzVector* &particle, int* index, KLFitter::Particles::ParticleType* ptype);

  /**
   * Return the parton at some index.
   * @param index The parton index
   * @return A pointer to the TLorentzVector of the parton.
   */
  TLorentzVector* Parton(int index);

  /**
   * Return the electron at some index.
   * @param index The electron index
   * @return A pointer to the TLorentzVector of the electron.
   */
  TLorentzVector* Electron(int index);

  /**
   * Return the muon at some index.
   * @param index The muon index
   * @return A pointer to the TLorentzVector of the muon.
   */
  TLorentzVector* Muon(int index);

  /**
   * Return the tau at some index.
   * @param index The tau index
   * @return A pointer to the TLorentzVector of the tau.
   */
  TLorentzVector* Tau(int index);

  /**
   * Return the boson at some index.
   * @param index The boson index
   * @return A pointer to the TLorentzVector of the boson.
   */
  TLorentzVector* Boson(int index);

  /**
   * Return the neutrino at some index.
   * @param index The neutrino index
   * @return A pointer to the TLorentzVector of the neutrino.
   */
  TLorentzVector* Neutrino(int index);

  /**
   * Return the photon at some index.
   * @param index The photon index
   * @return A pointer to the TLorentzVector of the photon.
   */
  TLorentzVector* Photon(int index);

  /**
   * Return the track at some index.
   * @param index The track index
   * @return A pointer to the TLorentzVector of the track.
   */
  TLorentzVector* Track(int index);

  /**
   * Return the number of particles.
   * @return The number of particles.
   */
  int NParticles() const { return static_cast<int>(m_jets.size() + m_electrons.size() + m_muons.size() + m_taus.size() + m_neutrinos.size() + m_bosons.size() + m_photons.size()) + m_tracks.size(); }

  /**
   * Return the number of particles of a certain type.
   * @param ptype The particle type.
   * @return The number of particles.
   */
  int NParticles(KLFitter::Particles::ParticleType ptype) const;

  /**
   * Return the name of a particle.
   * @param index The index of the particle.
   * @param ptype The type of the particle.
   * @return The name of the particle.
   */
  std::string NameParticle(int index, KLFitter::Particles::ParticleType ptype) const;

  /**
   * Return the name of a electron.
   * @param index The index of the electron.
   * @return The name of the electron.
   */
  std::string NameElectron(int index) const;

  /**
   * Return the name of a muon.
   * @param index The index of the muon.
   * @return The name of the muon.
   */
  std::string NameMuon(int index) const;

  /**
   * Return the name of a tau.
   * @param index The index of the tau.
   * @return The name of the tau.
   */
  std::string NameTau(int index) const;

  /**
   * Return the name of a boson.
   * @param index The index of the boson.
   * @return The name of the boson.
   */
  std::string NameBoson(int index) const;

  /**
   * Return the name of a neutrino.
   * @param index The index of the neutrino.
   * @return The name of the neutrino.
   */
  std::string NameNeutrino(int index) const;

  /**
   * Return the name of a photon.
   * @param index The index of the photon.
   * @return The name of the photon.
   */
  std::string NamePhoton(int index) const;

  /**
   * Return the index of the measured particle.
   * @param index The index of the model particle.
   * @return The index of the measured particle.
   */
  int JetIndex(int index) const;

  /**
   * Return the index of the measured particle.
   * @param index The index of the model particle.
   * @return The index of the measured particle.
   */
  int ElectronIndex(int index) const;

  /**
   * Return the index of the measured particle.
   * @param index The index of the model particle.
   * @return The index of the measured particle.
   */
  int MuonIndex(int index) const;

  /**
   * Return the index of the measured particle.
   * @param index The index of the model particle.
   * @return The index of the measured particle.
   */
  int PhotonIndex(int index) const;

  /**
   * Return the index of the measured particle.
   * @param index The index of the model particle.
   * @return The index of the measured particle.
   */
  int TrackIndex(int index) const;

  /**
   * Return the true flavor of a parton.
   * @param index The parton index
   * @return The parton true flavor.
   */
  TrueFlavorType TrueFlavor(int index) const { return static_cast<Particles::TrueFlavorType>(m_jets.at(index).GetTrueFlavor()); }

  /**
   * Return has the jet been b-tagged?
   * @param index The parton index
   * @return The parton b-tagging boolean.
   */
  bool IsBTagged(int index) const { return m_jets.at(index).GetIsBTagged(); }

  /**
   * Return the jet b-tagging efficiency.
   * @param index The parton index
   * @return The jet b-tagging efficiency.
   */
  double BTaggingEfficiency(int index) const { return m_jets.at(index).GetBTagEfficiency(); }

  /**
   * Return the jet b-tagging rejection.
   * @param index The parton index
   * @return The jet b-tagging rejection.
   */
  double BTaggingRejection(int index) const { return m_jets.at(index).GetBTagRejection(); }

  /**
   * Return the jet b-tagging weight.
   * @param index The parton index
   * @return The jet b-tagging weight.
   */
  double BTagWeight(int index) const { return m_jets.at(index).GetBTagWeight(); }

  /**
   * Return the bool of a set tagging weight.
   * @param index The parton index
   * @return The bool of a set tagging weight
   */
  bool BTagWeightSet(int index) const { return m_jets.at(index).GetBTagWeightIsSet(); }

  /**
   * Return the uncertainty of a particle with some index and type.
   * @param index The index of the particle
   * @param ptype The particle type.
   * @return The uncertaintie of the particle
   */
  const std::vector<double>* Uncertainties(int index, KLFitter::Particles::ParticleType ptype) const;

  /**
   * Return the detector eta of a particle with some index and type.
   * @param index The index of the particle
   * @param ptype The particle type.
   * @return The detector eta of the particle
   */
  double DetEta(int index, KLFitter::Particles::ParticleType ptype) const;

  /**
   * Return the charge of the lepton with some index and type.
   * @param index The index of the particle
   * @param ptype The particle type.
   * @return The charge of the lepton.
   */
  float LeptonCharge(int index, KLFitter::Particles::ParticleType ptype) const;

  /// Return the number of b-tags.
  int NBTags() const;

  /** @} */
  /** \name Member functions (Set)  */
  /** @{ */

  /**
   * Set has the jet been b-tagged?
   * @param index The parton index
   * @param isBTagged The parton b-tagging boolean.
   * @return An error flag.
   */
  int SetIsBTagged(int index, bool isBTagged);

  /**
   * Set the jet b-tagging efficiency.
   * @param index The parton index
   * @param btagEff The b-tagging efficiency.
   * @return An error flag.
   */
  int SetBTaggingEfficiency(int index, double btagEff);

  /**
   * Set the jet b-tagging rejection.
   * @param index The parton index
   * @param btagRej The b-tagging probability.
   * @return The jet b-tagging rejection.
   */
  int SetBTaggingRejection(int index, double btagRej);

  /**
   * Set the jet b-tagging weight.
   * @param index The parton index
   * @param btagweight The b-tagging weight.
   * @return An error flag.
   */
  int SetBTagWeight(int index, double btagweight);

  /**
   * Set bool for set jet b-tagging weight.
   * @param index The parton index
   * @param btagweightset The b-tagging probability.
   * @return An error flag.
   */
  int SetBTagWeightSet(int index, bool btagweightset);

  /** @} */
  /** \name Member functions (misc)  */
  /** @{ */

  /**
   * Add a particle to a list of particles.
   * @param particle A pointer to the particle.
   * @param DetEta The Detector Eta of the particle.
   * @param LepCharge The Charge of the particle.
   * @param ptype The type of particle.
   * @param name The name of the particle.
   * @param measuredindex The index of the associated measured particle.
   * @return An error code.
   */
  int AddParticle(const TLorentzVector& particle, double DetEta, float LepCharge, KLFitter::Particles::ParticleType ptype, std::string name = "", int measuredindex = -1);

  /**
   * DEPRECATED FUNCTION. This is an overloaded implementation of
   * the previous function using TLorentzVector pointers instead
   * of const ref. The usage of this function is deprecated and
   * it will be removed in the next major release. Please switch
   * to the above implementation.
   */
  int AddParticle(const TLorentzVector* const particle, double DetEta, float LepCharge, KLFitter::Particles::ParticleType ptype, std::string name = "", int measuredindex = -1);

  /**
   * Add a particle to a list of particles.
   * @param particle A pointer to the particle.
   * @param DetEta The Detector Eta of the particle.
   * @param ptype The type of particle.
   * @param name The name of the particle.
   * @param measuredindex The index of the associated measured particle.
   * @param isBtagged Has the particle been b-tagged?
   * @param bTagEff B-tagging efficiency of the particle.
   * @param bTagRej B-tagging rejection of the particle.
   * @param trueflav The true flavor (only for model particles).
   * @param btagweight The b tagger weight).
   * @return An error code.
   */
  int AddParticle(const TLorentzVector& particle, double DetEta, KLFitter::Particles::ParticleType ptype, std::string name = "", int measuredindex = -1, bool isBtagged = false, double bTagEff = -1., double bTagRej = -1., TrueFlavorType trueflav = kNone, double btagweight = 999);

  /**
   * DEPRECATED FUNCTION. This is an overloaded implementation of
   * the previous function using TLorentzVector pointers instead
   * of const ref. The usage of this function is deprecated and
   * it will be removed in the next major release. Please switch
   * to the above implementation.
   */
  int AddParticle(const TLorentzVector* const particle, double DetEta, KLFitter::Particles::ParticleType ptype, std::string name = "", int measuredindex = -1, bool isBtagged = false, double bTagEff = -1., double bTagRej = -1., TrueFlavorType trueflav = kNone, double btagweight = 999);

  /**
   * Add a particle to a list of particles.
   * @param particle A pointer to the particle.
   * @param ptype The type of particle.
   * @param name The name of the particle.
   * @param isBtagged Has the particle been b-tagged?
   * @param bTagEff B-tagging efficiency of the particle.
   * @param bTagRej B-tagging rejection of the particle.
   * @param measuredindex The index of the associated measured particle.
   * @param trueflav The true flavor (only for model particles).
   * @param btagweight The b tagger weight).
   * @return An error code.
   */
  int AddParticle(const TLorentzVector& particle, KLFitter::Particles::ParticleType ptype, std::string name = "", int measuredindex = -1, bool isBtagged = false, double bTagEff = -1., double bTagRej = -1., TrueFlavorType trueflav = kNone, double btagweight = 999);

  /**
   * DEPRECATED FUNCTION. This is an overloaded implementation of
   * the previous function using TLorentzVector pointers instead
   * of const ref. The usage of this function is deprecated and
   * it will be removed in the next major release. Please switch
   * to the above implementation.
   */
  int AddParticle(const TLorentzVector* const particle, KLFitter::Particles::ParticleType ptype, std::string name = "", int measuredindex = -1, bool isBtagged = false, double bTagEff = -1., double bTagRej = -1., TrueFlavorType trueflav = kNone, double btagweight = 999);

  /**
   * Add a particle to a list of particles (especially for model particles).
   * @param particle A pointer to the particle.
   * @param ptype The type of particle.
   * @param name The name of the particle.
   * @param measuredindex The index of the associated measured particle.
   * @param trueflav The true flavor (only for model particles).
   * @param btagweight The b tagger weight).
   * @return An error code.
   */
  int AddParticle(const TLorentzVector& particle, KLFitter::Particles::ParticleType ptype, std::string name, int measuredindex, TrueFlavorType trueflav, double btagweight = 999);

  /**
   * DEPRECATED FUNCTION. This is an overloaded implementation of
   * the previous function using TLorentzVector pointers instead
   * of const ref. The usage of this function is deprecated and
   * it will be removed in the next major release. Please switch
   * to the above implementation.
   */
  int AddParticle(const TLorentzVector* const particle, KLFitter::Particles::ParticleType ptype, std::string name, int measuredindex, TrueFlavorType trueflav, double btagweight = 999);

  /**
   * Add a particle to a list of particles.
   * @param particle A pointer to the particle.
   * @param ptype The type of particle.
   * @param name The name of the particle.
   * @param uncertainies The associated uncertainties.
   */
  int AddParticle(const TLorentzVector& particle, KLFitter::Particles::ParticleType ptype, std::string name = "", int measuredindex = -1, const std::vector<double>& uncertainies = std::vector<double>());

  /**
   * Removes a particle from a list of particles.
   * @param index The index of the particle.
   * @param ptype The type of the particle.
   * @return An error code.
   */
  int RemoveParticle(int index, KLFitter::Particles::ParticleType ptype);

  /**
   * Removes a particle from a list of particles.
   * @param name The name of the particle.
   * @return An error code.
   */
  int RemoveParticle(const std::string& name);

  /**
   * Return the const particle container of a type of particles
   * @param ptype The type of the particle.
   * @return The particle container.
   */
  const std::vector<TLorentzVector>* ParticleContainer(KLFitter::Particles::ParticleType ptype) const;

  /**
   * Return the particle container of a type of particles
   * @param ptype The type of the particle.
   * @return The particle container.
   */
  std::vector<TLorentzVector>* ParticleContainer(KLFitter::Particles::ParticleType ptype);

  /**
   * Return the (const) particle name container of a type of particles
   * @param ptype The type of the particle.
   * @return The particle name container.
   */
  const std::vector<std::string>* ParticleNameContainer(KLFitter::Particles::ParticleType ptype) const;

  /**
   * Return the particle name container of a type of particles
   * @param ptype The type of the particle.
   * @return The particle name container.
   */
  std::vector<std::string>* ParticleNameContainer(KLFitter::Particles::ParticleType ptype);

  /**
   * Checks if the index is within range.
   * @param container The particle container.
   * @param index The index of particle.
   * @return An error flag.
   */
  int CheckIndex(const std::vector<TLorentzVector>& container, int index) const;

  /** @} */

 private:
  /// Vector of all Particle::Jet objects.
  std::vector<Particle::Jet> m_jets;

  /// vector of all Particle::Electron objects.
  std::vector<Particle::Electron> m_electrons;

  /// Vector of all Particle::Muon objects.
  std::vector<Particle::Muon> m_muons;

  /// Vector of all Particle::Tau objects.
  std::vector<Particle::Tau> m_taus;

  /// Vector of all Particle::Neutrino objects.
  std::vector<Particle::Neutrino> m_neutrinos;

  /// Vector of all Particle::Boson objects.
  std::vector<Particle::Boson> m_bosons;

  /// Vector of all Particle::Photon objects.
  std::vector<Particle::Photon> m_photons;

  /// Vector of all Particle::Track objects.
  std::vector<Particle::Track> m_tracks;
};
}  // namespace KLFitter

inline KLFitter::Particles::ParticleType &operator++(KLFitter::Particles::ParticleType &ptype) {
  return ptype = KLFitter::Particles::ParticleType(ptype + 1);
}

#endif  // KLFITTER_PARTICLES_H_
