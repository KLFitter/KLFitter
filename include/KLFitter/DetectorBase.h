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

#ifndef DETECTORBASE
#define DETECTORBASE

// ---------------------------------------------------------

#include <string>

// ---------------------------------------------------------

/**
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter {
  class ResolutionBase;

  /**
   * \class KLFitter::DetectorBase
   * \brief A base class for describing detectors.
   * \author Kevin Kr&ouml;ninger
   *
   * This base class contains the energy resolution of different
   * objects. More information (angular resolutions, acceptance,
   * correections, etc.) can be added here.
   */
  class DetectorBase {
  public:
    /** \name Enumerators */
    /* @{ */

    /**
     * Enumerate for beam centre-of-mass energy
     */
    enum BeamCMEnergy {k7TeV,k8TeV, k10TeV};

    /* @} */

    /** \name Constructors and destructors */
    /* @{ */

    /**
     * The default constructor.
     * @param folder The folder with transfer function parameters.
     */
    DetectorBase(std::string folder = "");

    /**
     * The default destructor.
     */
    virtual ~DetectorBase();

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /**
     * Return the energy resolution of light jets.
     * @param eta The eta of the particle.
     * @return A pointer to the energy resolution object.
     */
    virtual KLFitter::ResolutionBase * ResEnergyLightJet(double eta = 0.) { return fResEnergyLightJet; };

    /**
     * Return the energy resolution of b jets.
     * @param eta The eta of the particle.
     * @return A pointer to the energy resolution object.
     */
    virtual KLFitter::ResolutionBase * ResEnergyBJet(double eta = 0.) { return fResEnergyBJet; };

    /**
     * Return the energy resolution of gluon jets.
     * @param eta The eta of the particle.
     * @return A pointer to the energy resolution object.
     */
    virtual KLFitter::ResolutionBase * ResEnergyGluonJet(double eta = 0.) { return fResEnergyGluonJet; };

    /**
     * Return the energy resolution of electrons.
     * @param eta The eta of the particle.
     * @return A pointer to the energy resolution object.
     */
    virtual KLFitter::ResolutionBase * ResEnergyElectron(double eta = 0.) { return fResEnergyElectron; };

    /**
     * Return the energy resolution of muons.
     * @param eta The eta of the particle.
     * @return A pointer to the energy resolution object.
     */
    virtual KLFitter::ResolutionBase * ResEnergyMuon(double eta = 0.) { return fResEnergyMuon; };

    /**
     * Return the energy resolution of photons.
     * @param eta The eta of the particle.
     * @return A pointer to the energy resolution object.
     */
    virtual KLFitter::ResolutionBase * ResEnergyPhoton(double eta = 0.) { return fResEnergyPhoton; };

    /**
     * Return the missing ET resolution.
     * @return A pointer to the missing ET resolution.
     */
    virtual KLFitter::ResolutionBase * ResMissingET() { return fResMissingET; };

    /**
     * Return the eta resolution of light jets.
     * @param eta The eta of the particle.
     * @return A pointer to the eta resolution object.
     */
    virtual KLFitter::ResolutionBase * ResEtaLightJet(double eta = 0.) { return fResEtaLightJet; };

    /**
     * Return the eta resolution of b jets.
     * @param eta The eta of the particle.
     * @return A pointer to the eta resolution object.
     */
    virtual KLFitter::ResolutionBase * ResEtaBJet(double eta = 0.) { return fResEtaBJet; };

    /**
     * Return the phi resolution of light jets.
     * @param eta The phi of the particle.
     * @return A pointer to the phi resolution object.
     */
    virtual KLFitter::ResolutionBase * ResPhiLightJet(double eta = 0.) { return fResPhiLightJet; };

    /**
     * Return the phi resolution of b jets.
     * @param eta The phi of the particle.
     * @return A pointer to the phi resolution object.
     */
    virtual KLFitter::ResolutionBase * ResPhiBJet(double eta = 0.) { return fResPhiBJet; };

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */

    /**
     * Set the energy resolution parameterization of b jets.
     * @param res A pointer to the resolution object.
     * @return An error code.
     */
    int SetResEnergyBJet(KLFitter::ResolutionBase * res);

    /**
     * Set the energy resolution parameterization of light jets.
     * @param res A pointer to the resolution object.
     * @return An error code.
     */
    int SetResEnergyLightJet(KLFitter::ResolutionBase * res);

    /**
     * Set the energy resolution parameterization of gluon jets.
     * @param res A pointer to the resolution object.
     * @return An error code.
     */
    int SetResEnergyGluonJet(KLFitter::ResolutionBase * res);

    /**
     * Set the energy resolution parameterization of electrons.
     * @param res A pointer to the resolution object.
     * @return An error code.
     */
    int SetResEnergyElectron(KLFitter::ResolutionBase * res);

    /**
     * Set the energy resolution parameterization of muons.
     * @param res A pointer to the resolution object.
     * @return An error code.
     */
    int SetResEnergyMuon(KLFitter::ResolutionBase * res);

    /**
     * Set the energy resolution parameterization of photons.
     * @param res A pointer to the resolution object.
     * @return An error code.
     */
    int SetResEnergyPhoton(KLFitter::ResolutionBase * res);

    /**
     * Set the missing ET resolution parameterization.
     * @param res A pointer to the resolution object.
     * @return An error code.
     */
    int SetResMissingET(KLFitter::ResolutionBase * res);

    /**
     * Set the beam centre-of-mass energy in the current detector.
     * @param beamenergy The beam energy.
     * @return An error code.
     */
    int SetBeamCMEnergy(KLFitter::DetectorBase::BeamCMEnergy beamenergy) {fBeamCMEnergy = beamenergy; return 1;};

    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */

    int Status();

    /**
     * Get the beam centre-of-mass energy in the current detector.
     * @return An error code.
     */
    KLFitter::DetectorBase::BeamCMEnergy GetBeamCMEnergy() {return fBeamCMEnergy;};
    /* @} */

  protected:
    /**
     * The energy resolution of light jets.
     */
    KLFitter::ResolutionBase * fResEnergyLightJet;

    /**
     * The energy resolution of b jets.
     */
    KLFitter::ResolutionBase * fResEnergyBJet;

    /**
     * The energy resolution of gluon jets.
     */
    KLFitter::ResolutionBase * fResEnergyGluonJet;

    /**
     * The energy resolution of electrons.
     */
    KLFitter::ResolutionBase * fResEnergyElectron;

    /**
     * The energy resolution of muons.
     */
    KLFitter::ResolutionBase * fResEnergyMuon;

    /**
     * The energy resolution of photons.
     */
    KLFitter::ResolutionBase * fResEnergyPhoton;

    /**
     * The missing ET resolution.
     */
    KLFitter::ResolutionBase * fResMissingET;

    /**
     * The eta resolution of light jets.
     */
    KLFitter::ResolutionBase * fResEtaLightJet;

    /**
     * The eta resolution of b jets.
     */
    KLFitter::ResolutionBase * fResEtaBJet;

    /**
     * The phi resolution of light jets.
     */
    KLFitter::ResolutionBase * fResPhiLightJet;

    /**
     * The phi resolution of b jets.
     */
    KLFitter::ResolutionBase * fResPhiBJet;

    /**
     * The current beam centre-of-mass energy in the detector
     */
    KLFitter::DetectorBase::BeamCMEnergy fBeamCMEnergy;
  };
}  // namespace KLFitter

// ---------------------------------------------------------

#endif
