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

#include <set>
#include <string>

// ---------------------------------------------------------

namespace KLFitter {
class ResolutionBase;

/**
 * Enum class to describe all possible resolution types that a
 * detector can hold. Resolution functions can be implemented for
 * all these types. If a detector does not implement a type, it
 * should return DetectorBase::ResolutionUndefined().
 */
enum class ResolutionType {
  EnergyLightJet,  ///< Energy resolution of light jets
  EnergyBJet,      ///< Energy resolution of b-jets
  EnergyGluonJet,  ///< Energy resolution of gluon jets
  EnergyElectron,  ///< Energy resolution of electrons
  EnergyMuon,      ///< Energy resolution of muons
  EnergyPhoton,    ///< Energy resolution of photons
  MissingET,       ///< Missing ET resolution
  EtaLightJet,     ///< Eta resolution of light jets
  EtaBJet,         ///< Eta resolution of b-jets
  PhiLightJet,     ///< Phi resolution of light jets
  PhiBJet,         ///< Phi resolution of b-jets
};

/**
 * A base class for describing detectors. This base class
 * contains accessors for the energy resolution functions for all
 * members of ResolutionType.
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
  /** @{ */

  /**
   * Return the energy resolution of light jets.
   * @return A pointer to the energy resolution object.
   */
  virtual ResolutionBase* ResEnergyLightJet(double /*eta*/) = 0;

  /**
   * Return the energy resolution of b jets.
   * @return A pointer to the energy resolution object.
   */
  virtual ResolutionBase* ResEnergyBJet(double /*eta*/) = 0;

  /**
   * Return the energy resolution of gluon jets.
   * @return A pointer to the energy resolution object.
   */
  virtual ResolutionBase* ResEnergyGluonJet(double /*eta*/) = 0;

  /**
   * Return the energy resolution of electrons.
   * @return A pointer to the energy resolution object.
   */
  virtual ResolutionBase* ResEnergyElectron(double /*eta*/) = 0;

  /**
   * Return the energy resolution of muons.
   * @return A pointer to the energy resolution object.
   */
  virtual ResolutionBase* ResEnergyMuon(double /*eta*/) = 0;

  /**
   * Return the energy resolution of photons.
   * @return A pointer to the energy resolution object.
   */
  virtual ResolutionBase* ResEnergyPhoton(double /*eta*/) = 0;

  /**
   * Return the missing ET resolution.
   * @return A pointer to the missing ET resolution.
   */
  virtual ResolutionBase* ResMissingET() = 0;

  /**
   * Return the eta resolution of light jets.
   * @return A pointer to the eta resolution object.
   */
  virtual ResolutionBase* ResEtaLightJet(double /*eta*/) = 0;

  /**
   * Return the eta resolution of b jets.
   * @return A pointer to the eta resolution object.
   */
  virtual ResolutionBase* ResEtaBJet(double /*eta*/) = 0;

  /**
   * Return the phi resolution of light jets.
   * @return A pointer to the phi resolution object.
   */
  virtual ResolutionBase* ResPhiLightJet(double /*eta*/) = 0;

  /**
   * Return the phi resolution of b jets.
   * @return A pointer to the phi resolution object.
   */
  virtual ResolutionBase* ResPhiBJet(double /*eta*/) = 0;

  /** @} */
  /** \name Member functions (misc)  */
  /** @{ */

  int Status();

  /**
   * Request a resolution type from the detector. The likelihoods
   * should define which resolution types are needed and make the
   * appropriate requests to the detector class. Requested
   * resolution types are then checked in the Status() function.
   * @param type The ResolutionType object
   */
  void RequestResolutionType(const ResolutionType& type);

  /** @} */

 protected:
  /**
   * Handle the call to a resolution object that is not defined,
   * and throw an exception. The derived detector classes are
   * expected to reimplement *all* resolution functions and, if
   * not supported by the detector, return this dummy object that
   * throws the exception accordingly.
   * @param type The type of resolution function
   * @return A nullptr, because an exception is thrown
   */
  ResolutionBase* ResolutionUndefined(const std::string& type);

  /**
   * Handle cases where the parameters for a resolution type are
   * not available. This prints an error message and throws an
   * exception of type std::invalid_argument.
   * @param type The name of the resolution function
   */
   void ResolutionParametersUnavailable(const std::string& type);

  /// Requested resolutions that the detector must provide.
  std::set<ResolutionType> res_type_requested;
};
}  // namespace KLFitter

#endif  // KLFITTER_DETECTORBASE_H_
