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
 * Particles::Type, and new particles are added via one
 * of the overloaded AddParticle() functions. It provides
 * functionality to retrieve information about the contained
 * particles via various accessor methods.
 */
class ParticleCollection final {
 public:
  /// The default constructor.
  ParticleCollection();

  /// The copy constructor.
  explicit ParticleCollection(const ParticleCollection& o);

  /// The (defaulted) destructor.
  ~ParticleCollection();

  /// The assignment operator.
  ParticleCollection& operator=(const ParticleCollection& o);

  /** \name Member functions (Get)  */
  /** @{ */

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
  TLorentzVector* Particle(int index, Particle::Type ptype);

  /**
   * Find a particle by name.
   * @param name The name of the particle.
   * @param particle A reference to the pointer to the 4-vector.
   * @param index The pointer to the index.
   * @param ptype The pointer to the particle type.
   * @return A flag (1: found, 0: not found).
   */
  int FindParticle(const std::string& name, TLorentzVector* &particle, int* index, Particle::Type* ptype);

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
  int NParticles() const { return static_cast<int>(jets.size() + electrons.size() + muons.size() + taus.size() + neutrinos.size() + bosons.size() + photons.size()) + tracks.size(); }

  /**
   * Return the number of particles of a certain type.
   * @param ptype The particle type.
   * @return The number of particles.
   */
  int NParticles(Particle::Type ptype) const;

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
  Particle::JetTrueFlavor TrueFlavor(int index) const { return jets.at(index).GetTrueFlavor(); }

  /**
   * Return has the jet been b-tagged?
   * @param index The parton index
   * @return The parton b-tagging boolean.
   */
  bool IsBTagged(int index) const { return jets.at(index).GetIsBTagged(); }

  /**
   * Return the jet b-tagging efficiency.
   * @param index The parton index
   * @return The jet b-tagging efficiency.
   */
  double BTaggingEfficiency(int index) const { return jets.at(index).GetBTagEfficiency(); }

  /**
   * Return the jet b-tagging rejection.
   * @param index The parton index
   * @return The jet b-tagging rejection.
   */
  double BTaggingRejection(int index) const { return jets.at(index).GetBTagRejection(); }

  /**
   * Return the jet b-tagging weight.
   * @param index The parton index
   * @return The jet b-tagging weight.
   */
  double BTagWeight(int index) const { return jets.at(index).GetBTagWeight(); }

  /**
   * Return the bool of a set tagging weight.
   * @param index The parton index
   * @return The bool of a set tagging weight
   */
  bool BTagWeightSet(int index) const { return jets.at(index).GetBTagWeightIsSet(); }

  /**
   * Return the uncertainty of a particle with some index and type.
   * @param index The index of the particle
   * @param ptype The particle type.
   * @return The uncertaintie of the particle
   */
  const std::vector<double>* Uncertainties(int index, Particle::Type ptype) const;

  /**
   * Return the detector eta of a particle with some index and type.
   * @param index The index of the particle
   * @param ptype The particle type.
   * @return The detector eta of the particle
   */
  double DetEta(int index, Particle::Type ptype) const;

  /**
   * Return the charge of the lepton with some index and type.
   * @param index The index of the particle
   * @param ptype The particle type.
   * @return The charge of the lepton.
   */
  float LeptonCharge(int index, Particle::Type ptype) const;

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
  int AddParticle(const TLorentzVector& particle, double DetEta, float LepCharge, Particle::Type ptype, std::string name = "", int measuredindex = -1);

  /**
   * DEPRECATED FUNCTION. This is an overloaded implementation of
   * the previous function using TLorentzVector pointers instead
   * of const ref. The usage of this function is deprecated and
   * it will be removed in the next major release. Please switch
   * to the above implementation.
   */
  int AddParticle(const TLorentzVector* const particle, double DetEta, float LepCharge, Particle::Type ptype, std::string name = "", int measuredindex = -1);

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
  int AddParticle(const TLorentzVector& particle, double DetEta, Particle::Type ptype, std::string name = "", int measuredindex = -1, bool isBtagged = false, double bTagEff = -1., double bTagRej = -1., Particle::JetTrueFlavor trueflav = Particle::JetTrueFlavor::kNone, double btagweight = 999);

  /**
   * DEPRECATED FUNCTION. This is an overloaded implementation of
   * the previous function using TLorentzVector pointers instead
   * of const ref. The usage of this function is deprecated and
   * it will be removed in the next major release. Please switch
   * to the above implementation.
   */
  int AddParticle(const TLorentzVector* const particle, double DetEta, Particle::Type ptype, std::string name = "", int measuredindex = -1, bool isBtagged = false, double bTagEff = -1., double bTagRej = -1., Particle::JetTrueFlavor trueflav = Particle::JetTrueFlavor::kNone, double btagweight = 999);

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
  int AddParticle(const TLorentzVector& particle, Particle::Type ptype, std::string name = "", int measuredindex = -1, bool isBtagged = false, double bTagEff = -1., double bTagRej = -1., Particle::JetTrueFlavor trueflav = Particle::JetTrueFlavor::kNone, double btagweight = 999);

  /**
   * DEPRECATED FUNCTION. This is an overloaded implementation of
   * the previous function using TLorentzVector pointers instead
   * of const ref. The usage of this function is deprecated and
   * it will be removed in the next major release. Please switch
   * to the above implementation.
   */
  int AddParticle(const TLorentzVector* const particle, Particle::Type ptype, std::string name = "", int measuredindex = -1, bool isBtagged = false, double bTagEff = -1., double bTagRej = -1., Particle::JetTrueFlavor trueflav = Particle::JetTrueFlavor::kNone, double btagweight = 999);

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
  int AddParticle(const TLorentzVector& particle, Particle::Type ptype, std::string name, int measuredindex, Particle::JetTrueFlavor trueflav, double btagweight = 999);

  /**
   * DEPRECATED FUNCTION. This is an overloaded implementation of
   * the previous function using TLorentzVector pointers instead
   * of const ref. The usage of this function is deprecated and
   * it will be removed in the next major release. Please switch
   * to the above implementation.
   */
  int AddParticle(const TLorentzVector* const particle, Particle::Type ptype, std::string name, int measuredindex, Particle::JetTrueFlavor trueflav, double btagweight = 999);

  /**
   * Add a particle to a list of particles.
   * @param particle A pointer to the particle.
   * @param ptype The type of particle.
   * @param name The name of the particle.
   * @param uncertainies The associated uncertainties.
   */
  int AddParticle(const TLorentzVector& particle, Particle::Type ptype, std::string name = "", int measuredindex = -1, const std::vector<double>& uncertainies = std::vector<double>());

  /**
   * Removes a particle from a list of particles.
   * @param index The index of the particle.
   * @param ptype The type of the particle.
   * @return An error code.
   */
  int RemoveParticle(int index, Particle::Type ptype);

  /**
   * Removes a particle from a list of particles.
   * @param name The name of the particle.
   * @return An error code.
   */
  int RemoveParticle(const std::string& name);

  /** @} */

  /** @{ */
  /** \name Particle containers */

  /// Vector of all Particle::Jet objects.
  std::vector<Particle::Jet> jets;

  /// vector of all Particle::Electron objects.
  std::vector<Particle::Electron> electrons;

  /// Vector of all Particle::Muon objects.
  std::vector<Particle::Muon> muons;

  /// Vector of all Particle::Tau objects.
  std::vector<Particle::Tau> taus;

  /// Vector of all Particle::Neutrino objects.
  std::vector<Particle::Neutrino> neutrinos;

  /// Vector of all Particle::Boson objects.
  std::vector<Particle::Boson> bosons;

  /// Vector of all Particle::Photon objects.
  std::vector<Particle::Photon> photons;

  /// Vector of all Particle::Track objects.
  std::vector<Particle::Track> tracks;

  /** @} */
};
}  // namespace KLFitter

inline KLFitter::Particle::Type &operator++(KLFitter::Particle::Type &ptype) {
  return ptype = KLFitter::Particle::Type(static_cast<int>(ptype) + 1);
}

#endif  // KLFITTER_PARTICLES_H_
