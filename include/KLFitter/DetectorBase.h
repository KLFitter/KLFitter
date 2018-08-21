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

  /** \name Enumerators */
  /* @{ */

  /// Enumerate for beam centre-of-mass energy
  enum BeamCMEnergy {k7TeV, k8TeV, k10TeV};

  /* @} */
  /** \name Member functions (Get)  */
  /* @{ */

  /**
   * Return the energy resolution of light jets.
   * @param eta The eta of the particle.
   * @return A pointer to the energy resolution object.
   */
  virtual ResolutionBase* ResEnergyLightJet(double /*eta*/) { return fResEnergyLightJet; }

  /**
   * Return the energy resolution of b jets.
   * @param eta The eta of the particle.
   * @return A pointer to the energy resolution object.
   */
  virtual ResolutionBase* ResEnergyBJet(double /*eta*/) { return fResEnergyBJet; }

  /**
   * Return the energy resolution of gluon jets.
   * @param eta The eta of the particle.
   * @return A pointer to the energy resolution object.
   */
  virtual ResolutionBase* ResEnergyGluonJet(double /*eta*/) { return fResEnergyGluonJet; }

  /**
   * Return the energy resolution of electrons.
   * @param eta The eta of the particle.
   * @return A pointer to the energy resolution object.
   */
  virtual ResolutionBase* ResEnergyElectron(double /*eta*/) { return fResEnergyElectron; }

  /**
   * Return the energy resolution of muons.
   * @param eta The eta of the particle.
   * @return A pointer to the energy resolution object.
   */
  virtual ResolutionBase* ResEnergyMuon(double /*eta*/) { return fResEnergyMuon; }

  /**
   * Return the energy resolution of photons.
   * @param eta The eta of the particle.
   * @return A pointer to the energy resolution object.
   */
  virtual ResolutionBase* ResEnergyPhoton(double /*eta*/) { return fResEnergyPhoton; }

  /**
   * Return the missing ET resolution.
   * @return A pointer to the missing ET resolution.
   */
  virtual ResolutionBase* ResMissingET() { return fResMissingET; }

  /**
   * Return the eta resolution of light jets.
   * @param eta The eta of the particle.
   * @return A pointer to the eta resolution object.
   */
  virtual ResolutionBase* ResEtaLightJet(double /*eta*/) { return fResEtaLightJet; }

  /**
   * Return the eta resolution of b jets.
   * @param eta The eta of the particle.
   * @return A pointer to the eta resolution object.
   */
  virtual ResolutionBase* ResEtaBJet(double /*eta*/) { return fResEtaBJet; }

  /**
   * Return the phi resolution of light jets.
   * @param eta The phi of the particle.
   * @return A pointer to the phi resolution object.
   */
  virtual ResolutionBase* ResPhiLightJet(double /*eta*/) { return fResPhiLightJet; }

  /**
   * Return the phi resolution of b jets.
   * @param eta The phi of the particle.
   * @return A pointer to the phi resolution object.
   */
  virtual ResolutionBase* ResPhiBJet(double /*eta*/) { return fResPhiBJet; }

  /* @} */
  /** \name Member functions (Set)  */
  /* @{ */

  /**
    * Set the energy resolution parameterization of b jets.
    * @param res A pointer to the resolution object.
    * @return An error code.
    */
  int SetResEnergyBJet(ResolutionBase* res);

  /**
    * Set the energy resolution parameterization of light jets.
    * @param res A pointer to the resolution object.
    * @return An error code.
    */
  int SetResEnergyLightJet(ResolutionBase* res);

  /**
    * Set the energy resolution parameterization of gluon jets.
    * @param res A pointer to the resolution object.
    * @return An error code.
    */
  int SetResEnergyGluonJet(ResolutionBase* res);

  /**
    * Set the energy resolution parameterization of electrons.
    * @param res A pointer to the resolution object.
    * @return An error code.
    */
  int SetResEnergyElectron(ResolutionBase* res);

  /**
    * Set the energy resolution parameterization of muons.
    * @param res A pointer to the resolution object.
    * @return An error code.
    */
  int SetResEnergyMuon(ResolutionBase* res);

  /**
    * Set the energy resolution parameterization of photons.
    * @param res A pointer to the resolution object.
    * @return An error code.
    */
  int SetResEnergyPhoton(ResolutionBase* res);

  /**
    * Set the missing ET resolution parameterization.
    * @param res A pointer to the resolution object.
    * @return An error code.
    */
  int SetResMissingET(ResolutionBase* res);

  /**
    * Set the beam centre-of-mass energy in the current detector.
    * @param beamenergy The beam energy.
    * @return An error code.
    */
  int SetBeamCMEnergy(BeamCMEnergy beamenergy) { fBeamCMEnergy = beamenergy; return 1; }

  /* @} */
  /** \name Member functions (misc)  */
  /* @{ */

  int Status();

  /**
    * Get the beam centre-of-mass energy in the current detector.
    * @return An error code.
    */
  BeamCMEnergy GetBeamCMEnergy() { return fBeamCMEnergy; }
  /* @} */

 protected:
  /**
   * The energy resolution of light jets.
   */
  ResolutionBase* fResEnergyLightJet;

  /**
   * The energy resolution of b jets.
   */
  ResolutionBase* fResEnergyBJet;

  /**
   * The energy resolution of gluon jets.
   */
  ResolutionBase* fResEnergyGluonJet;

  /**
   * The energy resolution of electrons.
   */
  ResolutionBase* fResEnergyElectron;

  /**
   * The energy resolution of muons.
   */
  ResolutionBase* fResEnergyMuon;

  /**
   * The energy resolution of photons.
   */
  ResolutionBase* fResEnergyPhoton;

  /**
   * The missing ET resolution.
   */
  ResolutionBase* fResMissingET;

  /**
   * The eta resolution of light jets.
   */
  ResolutionBase* fResEtaLightJet;

  /**
   * The eta resolution of b jets.
   */
  ResolutionBase* fResEtaBJet;

  /**
   * The phi resolution of light jets.
   */
  ResolutionBase* fResPhiLightJet;

  /**
   * The phi resolution of b jets.
   */
  ResolutionBase* fResPhiBJet;

  /**
   * The current beam centre-of-mass energy in the detector
   */
  BeamCMEnergy fBeamCMEnergy;
};
}  // namespace KLFitter

#endif  // KLFITTER_DETECTORBASE_H_
