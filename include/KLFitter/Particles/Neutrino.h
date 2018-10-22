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

#ifndef KLFITTER_PARTICLES_NEUTRINO_H_
#define KLFITTER_PARTICLES_NEUTRINO_H_

#include "KLFitter/Particles/Base.h"

namespace KLFitter {
namespace Particles {
/**
 * Implementation of a neutrino class. This inherits a generic
 * structure from Particles::Base, reimplements all purely virtual
 * functions and adds functionality specific to neutrinos.
 */
class Neutrino : public Base {
 public:
  /**
   * Neutrino constructor taking name and four momentum.
   * @param name The name of the particle
   * @param p4 The four momentum
   */
  Neutrino(const std::string& name, const TLorentzVector& p4)
    : m_name(name)
    , m_p4(p4) {
    // empty
  }

  /// The (defaulted) destructor.
  ~Neutrino() = default;

  /// The (defaulted) copy constructor.
  Neutrino(const Neutrino&) = default;

  /// The (defaulted) assignment operator.
  Neutrino& operator=(const Neutrino&) = default;

  /** \name Inherited, reimplemented methods */
  /** @{ */

  /// Const accessor for the particle's name.
  const std::string& GetName() const override { return m_name; }

  /// Const accessor to the assigned identifier of the particle.
  int GetIdentifier() const override { return m_identifier; }

  /// Const accessor for the particle's four momentum.
  const TLorentzVector& GetP4() const override { return m_p4; }

  /// Non-const accessor for the particle's four momentum.
  TLorentzVector& GetP4() override { return m_p4; }

  /// Set the identifier of the particle.
  void SetIdentifier(int id) override { m_identifier = id; }

  /// Set the four momentum.
  void SetP4(const TLorentzVector& p4) override { m_p4 = p4; }

  /** @} */

 private:
  /// String with the name of the particle.
  std::string m_name;

  /// Assigned identifier of the particle.
  int m_identifier{-1};

  /// Four vector of the particle.
  TLorentzVector m_p4;
};
}  // namespace Particles
}  // namespace KLFitter

#endif  // KLFITTER_PARTICLES_NEUTRINO_H_
