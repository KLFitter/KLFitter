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

#ifndef KLFITTER_DETECTORBASE_H_
#define KLFITTER_DETECTORBASE_H_

#include <string>

// ---------------------------------------------------------

namespace KLFitter {
class ResolutionBase;

/**
 * A base class for describing detectors. This base class
 * contains the energy resolution of different objects. More
 * information (angular resolutions, acceptance, correections,
 * etc.) can be added here.
 */
class DetectorBase {
 public:
  /**
   * The default constructor.
   * @param folder The folder with transfer function parameters.
   */
  explicit DetectorBase(std::string folder = "");

  /// The (defaulted) destructor.
  virtual ~DetectorBase();

  /** \name Member functions (Get)  */
  /* @{ */

  /**
   * Return the energy resolution of light jets.
   * @param eta The eta of the particle.
   * @return A pointer to the energy resolution object.
   */
  virtual ResolutionBase* ResEnergyLightJet(double /*eta*/) { return nullptr; }

  /**
   * Return the energy resolution of b jets.
   * @param eta The eta of the particle.
   * @return A pointer to the energy resolution object.
   */
  virtual ResolutionBase* ResEnergyBJet(double /*eta*/) { return nullptr; }

  /**
   * Return the energy resolution of gluon jets.
   * @param eta The eta of the particle.
   * @return A pointer to the energy resolution object.
   */
  virtual ResolutionBase* ResEnergyGluonJet(double /*eta*/) { return nullptr; }

  /**
   * Return the energy resolution of electrons.
   * @param eta The eta of the particle.
   * @return A pointer to the energy resolution object.
   */
  virtual ResolutionBase* ResEnergyElectron(double /*eta*/) { return nullptr; }

  /**
   * Return the energy resolution of muons.
   * @param eta The eta of the particle.
   * @return A pointer to the energy resolution object.
   */
  virtual ResolutionBase* ResEnergyMuon(double /*eta*/) { return nullptr; }

  /**
   * Return the energy resolution of photons.
   * @param eta The eta of the particle.
   * @return A pointer to the energy resolution object.
   */
  virtual ResolutionBase* ResEnergyPhoton(double /*eta*/) { return nullptr; }

  /**
   * Return the missing ET resolution.
   * @return A pointer to the missing ET resolution.
   */
  virtual ResolutionBase* ResMissingET() { return nullptr; }

  /**
   * Return the eta resolution of light jets.
   * @param eta The eta of the particle.
   * @return A pointer to the eta resolution object.
   */
  virtual ResolutionBase* ResEtaLightJet(double /*eta*/) { return nullptr; }

  /**
   * Return the eta resolution of b jets.
   * @param eta The eta of the particle.
   * @return A pointer to the eta resolution object.
   */
  virtual ResolutionBase* ResEtaBJet(double /*eta*/) { return nullptr; }

  /**
   * Return the phi resolution of light jets.
   * @param eta The phi of the particle.
   * @return A pointer to the phi resolution object.
   */
  virtual ResolutionBase* ResPhiLightJet(double /*eta*/) { return nullptr; }

  /**
   * Return the phi resolution of b jets.
   * @param eta The phi of the particle.
   * @return A pointer to the phi resolution object.
   */
  virtual ResolutionBase* ResPhiBJet(double /*eta*/) { return nullptr; }

  /* @} */
  /** \name Member functions (misc)  */
  /* @{ */

  int Status();

  /* @} */
};
}  // namespace KLFitter

#endif  // KLFITTER_DETECTORBASE_H_
