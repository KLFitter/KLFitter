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

#ifndef KLFITTER_PARTICLES_PARTON_H_
#define KLFITTER_PARTICLES_PARTON_H_

#include "KLFitter/Particles/Base.h"

namespace KLFitter {
namespace Particles {
/// An enumerator class for the true parton flavor.
enum class PartonTrueFlavor {kLight,     ///< Light quark
                             kB,         ///< B-quark
                             kLightUp,   ///< Up-type light quark
                             kLightDown, ///< Down-type light quark
                             kNone       ///< Not specified
};

/**
 * Implementation of a parton class. This inherits a generic
 * structure from Particles::Base, reimplements all purely virtual
 * functions and adds functionality specific to partons.
 */
class Parton : public Base {
 public:
  /**
   * Parton constructor taking name and four momentum.
   * @param name The name of the particle
   * @param p4 The four momentum
   */
  Parton(const std::string& name, const TLorentzVector& p4)
      : m_name(name)
      , m_p4(p4) {
    // empty
  }

  /**
   * Parton constructor taking name, four momentum, and whether the
   * parton is b-tagged or not.
   * @param name The name of the particle
   * @param p4 The four momentum
   * @param is_tagged Boolean whether the parton is b-tagged
   */
  Parton(const std::string& name, const TLorentzVector& p4, bool is_tagged)
      : m_name(name)
      , m_p4(p4)
      , m_is_btagged(is_tagged) {
    // empty
  }

  /// The (defaulted) destructor.
  ~Parton() = default;

  /// The (defaulted) copy constructor.
  Parton(const Parton&) = default;

  /// The (defaulted) assignment operator.
  Parton& operator=(const Parton&) = default;

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
  /** \name Get methods */
  /** @{ */

  /// Get the detector eta of the parton.
  double GetDetEta() const { return m_det_eta; }

  /// Get the b-tagging efficiency.
  double GetBTagEfficiency() const { return m_btag_eff; }

  /// Get the b-tagging rejection.
  double GetBTagRejection() const { return m_btag_rej; }

  /// Get the b-tagging weight.
  double GetBTagWeight() const { return m_btag_weight; }

  /// Get whether the b-tag weight is set.
  bool GetBTagWeightIsSet() const { return m_btag_weight_set; }

  /// Get whether the parton is b-tagged.
  bool GetIsBTagged() const { return m_is_btagged; }

  /// Get the true parton flavor.
  PartonTrueFlavor GetTrueFlavor() const { return m_true_flavor; }

  /** @} */
  /** \name Set methods */
  /** @{ */

  /// Set the detector eta.
  void SetDetEta(double eta) { m_det_eta = eta; }

  /// Set the b-tagging efficiency.
  void SetBTagEfficiency(double val) { m_btag_eff = val; }

  /// Set the b-tagging rejection.
  void SetBTagRejection(double val) { m_btag_rej = val; }

  /// Set the b-tagging weight.
  void SetBTagWeight(double val) { m_btag_weight = val; m_btag_weight_set = true; }

  /// Set whether the b-tag weight is set.
  void SetBTagWeightIsSet(bool set) { m_btag_weight_set = set; }

  /// Set whether the parton is b-tagged.
  void SetIsBTagged(bool is_tagged) { m_is_btagged = is_tagged; }

  /// Set the true parton flavor.
  void SetTrueFlavor(PartonTrueFlavor flavor) { m_true_flavor = flavor; }

  /** @} */

 private:
  /// String with the name of the particle.
  std::string m_name;

  /// Assigned identifier of the particle.
  int m_identifier{-1};

  /// Four vector of the particle.
  TLorentzVector m_p4{};

  double m_det_eta{0.};
  double m_btag_eff{0.};
  double m_btag_rej{0.};
  double m_btag_weight{0.};
  bool m_btag_weight_set{false};
  bool m_is_btagged{false};
  PartonTrueFlavor m_true_flavor{PartonTrueFlavor::kNone};
};
}  // namespace Particles
}  // namespace KLFitter

#endif  // KLFITTER_PARTICLES_PARTON_H_
