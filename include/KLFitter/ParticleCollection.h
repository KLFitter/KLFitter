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
#include "KLFitter/Particles/Boson.h"
#include "KLFitter/Particles/Electron.h"
#include "KLFitter/Particles/Parton.h"
#include "KLFitter/Particles/Muon.h"
#include "KLFitter/Particles/Neutrino.h"
#include "KLFitter/Particles/Photon.h"
#include "KLFitter/Particles/Tau.h"
#include "KLFitter/Particles/Track.h"

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

  /* @} */
  /** \name Add and remove particles  */
  /* @{ */

  /**
   * Add particle of type parton to the collection.
   * @param p Const reference to the particle object.
   */
  void AddParticle(const Particles::Parton& p);

  /**
   * Add particle of type electron to the collection.
   * @param p Const reference to the particle object.
   */
  void AddParticle(const Particles::Electron& p);

  /**
   * Add particle of type muon to the collection.
   * @param p Const reference to the particle object.
   */
  void AddParticle(const Particles::Muon& p);

  /**
   * Add particle of type photon to the collection.
   * @param p Const reference to the particle object.
   */
  void AddParticle(const Particles::Photon& p);

  /**
   * Add particle of type tau to the collection.
   * @param p Const reference to the particle object.
   */
  void AddParticle(const Particles::Tau& p);

  /**
   * Add particle of type neutrino to the collection.
   * @param p Const reference to the particle object.
   */
  void AddParticle(const Particles::Neutrino& p);

  /**
   * Add particle of type boson to the collection.
   * @param p Const reference to the particle object.
   */
  void AddParticle(const Particles::Boson& p);

  /**
   * Add particle of type track to the collection.
   * @param p Const reference to the particle object.
   */
  void AddParticle(const Particles::Track& p);

  /**
   * Removes a particle from a list of particles.
   * @param ptype The type of the particle.
   * @param index The index of the particle.
   */
  void RemoveParticle(Particles::Type ptype, size_t index);

  /** @} */
  /** \name Get particles and four-momenta */
  /** @{ */

  /**
   * Find a particle of any type by its name. Returns either the
   * pointer to the particle or a null pointer.
   * @param name The name of the particle.
   * @return Const pointer to the particle.
   */
  const Particles::Base* FindParticle(const std::string& name) const;

  /**
   * Find a particle of a given type by its name. Returns either
   * the pointer to the particle or a null pointer.
   * @param ptype The particle type.
   * @param name The name of the particle.
   * @return Const pointer to the particle.
   */
  const Particles::Base* FindParticle(Particles::Type ptype, const std::string& name) const;

  /**
   * Get the four-vector of a particle. This needs the type of
   * the particle and its index in the particle collection.
   * @param ptype The type of the particle.
   * @param index The index of the particle.
   * @return A (const) pointer to the particle four-momentum.
   */
  const TLorentzVector* GetP4(Particles::Type ptype, size_t index) const;

  /**
   * Get the four-vector of a particle. This needs the type of
   * the particle and its index in the particle collection.
   * @param ptype The type of the particle.
   * @param index The index of the particle.
   * @return A pointer to the particle four-momentum.
   */
  TLorentzVector* GetP4(Particles::Type ptype, size_t index);

  /** @} */
  /** \name Miscellaneous functions */
  /** @{ */

  /**
   * Return the total number of particles in the collection.
   * @return The number of particles.
   */
  size_t NParticles() const;

  /**
   * Return the number of particles of a certain type.
   * @param ptype The particle type.
   * @return The number of particles.
   */
  size_t NParticles(Particles::Type ptype) const;

  /**
   * Return the number of b-tagged Particles::Parton objects in
   * this collection (stored in #partons).
   */
  size_t NBTags() const;

  /** @} */
  /** \name Particle containers */
  /** @{ */

  /// Vector of all Particles::Parton objects.
  std::vector<Particles::Parton> partons;

  /// vector of all Particles::Electron objects.
  std::vector<Particles::Electron> electrons;

  /// Vector of all Particles::Muon objects.
  std::vector<Particles::Muon> muons;

  /// Vector of all Particles::Tau objects.
  std::vector<Particles::Tau> taus;

  /// Vector of all Particles::Neutrino objects.
  std::vector<Particles::Neutrino> neutrinos;

  /// Vector of all Particles::Boson objects.
  std::vector<Particles::Boson> bosons;

  /// Vector of all Particles::Photon objects.
  std::vector<Particles::Photon> photons;

  /// Vector of all Particles::Track objects.
  std::vector<Particles::Track> tracks;

  /** @} */
};
}  // namespace KLFitter

inline KLFitter::Particles::Type &operator++(KLFitter::Particles::Type &ptype) {
  return ptype = KLFitter::Particles::Type(static_cast<int>(ptype) + 1);
}

#endif  // KLFITTER_PARTICLES_H_
