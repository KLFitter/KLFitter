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

#ifndef KLFITTER_PARTICLES_BASE_H_
#define KLFITTER_PARTICLES_BASE_H_

#include <string>

#include "TLorentzVector.h"

namespace KLFitter {
/**
 * Namespace within the KLFitter namespace that holds classes for
 * various types of particles. The purely virtual class
 * Particles::Base is the parent class of all others.
 */
namespace Particles {
/// An enumerator class for the particle type.
enum class Type {kParton,    ///< Particle type parton
                 kElectron,  ///< Particle type electron
                 kMuon,      ///< Particle type muon
                 kTau,       ///< Particle type tau
                 kNeutrino,  ///< Particle type neutrino
                 kBoson,     ///< Particle type boson
                 kPhoton,    ///< Particle type photon
                 kTrack      ///< Particle type track
};

/**
 * The base class for all other particles. This only provides a
 * very generic structure, that other particles can inherit from.
 * Note that this is a purely virtual class, so it cannot be
 * instantiated.
 */
class Base {
 public:
  /// Default constructor.
  Base() = default;

  /// The (defaulted) destructor.
  virtual ~Base() = default;

  /// The (defaulted) copy constructor.
  Base(const Base&) = default;

  /// The (defaulted) assignment operator.
  Base& operator=(const Base&) = default;

  /** @} */
  /** \name Get methods */
  /** @{ */

  /// Const accessor to the particle's name.
  virtual const std::string& GetName() const = 0;

  /// Const accessor to the assigned identifier of the particle.
  virtual int GetIdentifier() const = 0;

  /// Const accessor to the particle's four momentum.
  virtual const TLorentzVector& GetP4() const = 0;

  /// Non-const accessor to the particle's four momentum.
  virtual TLorentzVector& GetP4() = 0;

  /** @} */
  /** \name Set methods */
  /** @{ */

  /// Set the identifier of the particle.
  virtual void SetIdentifier(int) = 0;

  /// Set the four momentum.
  virtual void SetP4(const TLorentzVector&) = 0;

  /** @} */
};
}  // namespace Particles
}  // namespace KLFitter

#endif  // KLFITTER_PARTICLES_BASE_H_
