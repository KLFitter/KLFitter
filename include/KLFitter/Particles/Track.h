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

#ifndef KLFITTER_PARTICLES_TRACK_H_
#define KLFITTER_PARTICLES_TRACK_H_

#include "KLFitter/Particles/Base.h"

namespace KLFitter {
namespace Particles {
/**
 * Implementation of a track class. This inherits a generic
 * structure from Particles::Base, reimplements all purely virtual
 * functions and adds functionality specific to tracks.
 */
class Track : public Base {
 public:
  /**
   * Track constructor taking name and four momentum.
   * @param name The name of the particle
   * @param p4 The four momentum
   */
  Track(const std::string& name, const TLorentzVector& p4)
    : m_name(name)
    , m_p4(p4) {
    // empty
  }

  /// The (defaulted) destructor.
  ~Track() = default;

  /// The (defaulted) copy constructor.
  Track(const Track&) = default;

  /// The (defaulted) assignment operator.
  Track& operator=(const Track&) = default;

  /** \name Inherited, reimplemented methods */
  /** @{ */

  /// Const accessor for the particle's name.
  const std::string& GetName() const override { return m_name; }

  /// Const accessor to the assigned identifier of the particle.
  unsigned int GetIdentifier() const override { return m_identifier; }

  /// Const accessor for the particle's four momentum.
  const TLorentzVector& GetP4() const override { return m_p4; }

  /// Non-const accessor for the particle's four momentum.
  TLorentzVector& GetP4() override { return m_p4; }

  /// Set the identifier of the particle.
  void SetIdentifier(unsigned int id) override { m_identifier = id; }

  /// Set the four momentum.
  void SetP4(const TLorentzVector& p4) override { m_p4 = p4; }

  /** @} */
  /** \name Get methods */
  /** @{ */

  /// Const accessor to the uncertainties.
  const std::vector<double>& GetUncertainties() const { return m_uncertainties; }

  /// Non-const accessor to the uncertainties.
  std::vector<double>& GetUncertainties() { return m_uncertainties; }

  /** @} */
  /** \name Set methods */
  /** @{ */

  /// Set the uncertainties.
  void SetUncertainties(const std::vector<double>& unc) { m_uncertainties = unc; }

  /** @} */

 private:
  /// String with the name of the particle.
  std::string m_name;

  /// Assigned identifier of the particle.
  unsigned int m_identifier{0};

  /// Four vector of the particle.
  TLorentzVector m_p4;

  std::vector<double> m_uncertainties;
};
}  // namespace Particles
}  // namespace KLFitter

#endif  // KLFITTER_PARTICLES_TRACK_H_
